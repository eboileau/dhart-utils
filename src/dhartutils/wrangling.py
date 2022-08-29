#! /usr/bin/env python3

"""Data wrangling wrapper to prepare and reformat GEO 
   supplementary files into h5ad for input into DHART. 
"""

import string
import sys
import logging
import argparse
import re
import json
import tarfile

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad

from pathlib import Path

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

import dhartutils.utils.utils as utils
import dhartutils.utils.data_utils as data_utils

logger = logging.getLogger(__name__)


# defined as in accession.py
data_choices = ['single-cell RNA-Seq', 
                'bulk RNA-Seq', 
                'microarray', 
                'ChIP-Seq', 
                'ATAC-Seq']


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Data wrangling wrapper.""")

    parser.add_argument('-d', '--dest', help="""Path to directory where GEO files reside, 
                        including output from accession script and/or pre-processed bulk 
                        RNA-seq count matrix.""", type=Path, default=Path(__file__).parent)
    
    parser.add_argument('-n', '--name', help="""Name of H5AD output.""", type=str, 
                        default='adata')

    parser.add_argument('-fmt', '--input-format', help="""Input format: either MEX (e.g. TSV
                        and MTX), or TXT (e.g. TXT, TSV, or CSV). Market Exchange (MEX) format
                        is for gene-barcode matrix output from Cell Ranger; TXT format is 
                        for read count matrices (single file, genes x barcodes/wells/samples).""", 
                        type=str, choices=['MEX', 'TXT'], default=None)
    
    parser.add_argument('--mex-gene', help="""File name pattern, typically "features" or
                        "genes", used with [--input-format MEX]. Silently ignored if
                        [--input-format TXT].""", type=str, default='feature')
    
    parser.add_argument('--mex-barcode', help="""File name pattern, typically "barcodes",
                        used with [--input-format MEX]. Silently ignored if
                        [--input-format TXT].""", type=str, default='barcode')
    
    parser.add_argument('--mex-matrix', help="""File name pattern, typically "matrix",
                        used with [--input-format MEX]. Silently ignored if
                        [--input-format TXT].""", type=str, default='matrix')

    parser.add_argument('--mex-gene-ncols', help="""Number of feature file columns, 
                        used with [--input-format MEX]. If [--mex-gene-ncols 1], it should
                        be "gene_symbol", otherwise they must be in the following 
                        order: ["gene_id", "gene_symbol", "feature_type"].
                        Silently ignored if [--input-format TXT].""", type=int, 
                        choices=[1, 2, 3], default=2)
    
    parser.add_argument('--txt-ext', help="""Extension that indicates the file type. 
                        If None, uses extension of filename, see the anndata read_txt 
                        function. Silently ignored if [--input-format MEX].""", 
                        type=str, default=None)
    
    parser.add_argument('--txt-delimiter', help="""If None, will split at arbitrary 
                        number of white spaces, see the anndata read_txt function.
                        Silently ignored if [--input-format MEX].""", type=str,
                        default=None)
    
    parser.add_argument('--txt-first-column-names', help="""Assume the first column 
                        stores row names. This is only necessary if these are not 
                        strings: strings in the first column are automatically assumed 
                        to be row names. See anndata read_txt function.
                        Silently ignored if [--input-format MEX].""", 
                        action='store_true')
    
    parser.add_argument('--txt-first-row-names', help="""Assume the first row stores 
                        column names (barcodes/wells). The behaviour of this flag is 
                        mostly dependent on the anndata read_txt function. In general, 
                        it might wrongly be read into X (and var), e.g. if column names 
                        are well numbers. If True, then column names are read separately 
                        and re-assigned, and this depends on [--txt-delimiter]. If False, 
                        fall back entirely to read_txt.. Silently ignored 
                        if [--input-format MEX].""", action='store_true')
    
    parser.add_argument('--txt-gene-cols', help="""Content of first column(s), 
                        used with [--input-format TXT]. [--txt-gene-cols] can either be
                        "gene_id", "gene_symbol", or both, but then they must be in the 
                        following order: ["gene_id", "gene_symbol"]. Any additional column
                        will cause an error. Silently ignored if [--input-format MEX].""", 
                        type=str, choices=["gene_id", "gene_symbol", "both"], 
                        default="gene_id")
    
    parser.add_argument('--gene-var', help="""Column to use as "var_names". Can be used 
                        with either [--input-format] options.""", type=str, 
                        choices=['gene_symbol', 'gene_id'], default='gene_symbol')
    
    parser.add_argument('-p', '--pattern', help="""A space separated list of patterns 
                        used to glob files. Can be used with either [--input-format] options.
                        Search is performed independently (and before) file name pattern 
                        search for [--input-format MEX].""", type=str, default='', nargs='*')
    
    parser.add_argument('--normalised', help="""If this flag is present, then 
                        files are treated as normalised input. This has currently no effect:
                        for single-cell (presumably TXT fmt) wrangling is the same, for bulk,
                        this overrides [--do-normalise], if set.""", action='store_true')
    
    parser.add_argument('--do-normalise', help="""If this flag is present, then 
                        output is written as normalised (only for bulk, otherwise ignored). 
                        Input files must be raw.""", action='store_true')

    parser.add_argument('--normalization-method', help="""Set the normalization method
                        for bulk RNA-seq only.""", required='--do-normalise' in sys.argv, 
                        type=str, choices=['TMM', 'DESeq', 'TPM'], default='TMM')
    
    parser.add_argument('--gene-lengths', help="""If [--do-normalise] is set, and 
                        [--normalization-method] is TPM, gene lengths required for 
                        normalization. A single column (one length per row), wtihout 
                        header, where genes are in the same order as the read count
                        matrix.""", type=str, default=None)
    
    parser.add_argument('--mygene-scopes', help="""Type of identifiers as defined 
                        by MyGene.py. Refer to official MyGene.info docs for full list of 
                        fields. This is used if [--mex-gene-ncols] is 1, or if 
                        [--txt-gene-cols] is not both, to fill in either gene ids or 
                        symbols. If ids or symbols have one or more underscores, e.g.
                        "ZYX__chr7", then this is automatically converted to "ZXY". Any
                        other non-standard format needs to be pre-formatted accordingly.""", 
                        type=str, default='ensembl.gene')
    
    parser.add_argument('--replicate-column', help="""Column name from "observations.tab.gz"
                        to use as replicate, e.g. "characteristics_ch1.2.biological replicate".
                        By default, "replicate" is used to find a suitable column (substring
                        match).""", type=str, default="replicate")
    
    parser.add_argument('--cluster-column', help="""Column name from "observations.tab.gz"
                        to use for cell type or cluster (for scRNA-seq), e.g. 
                        "characteristics_ch1.1.cell type". By default, "cluster" is used 
                        to find a suitable column (substring match).""", type=str, 
                        default="cluster")
        
    utils.add_logging_options(parser)
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    utils.update_logging(args)
    
    msg = f'[wrangling]: {" ".join(sys.argv)}'
    logger.info(msg)
    
    # load json config and observations dataframe (from accession.py)
    fl = list(filter(Path.is_file, args.dest.glob('*.json')))
    if len(fl) != 1:
        msg = "[--dest] must include a single configuration file (*.json). " \
              f"{len(fl)} json file(s) were found! Check accession.py logs. " \
              "Remove non-essential files or specify a different directory. " \
              "Terminating program execution!"
        logger.critical(msg)
        return
    with open(fl[0], "r") as f:
        config = json.load(f)

    filen = Path(args.dest, "observations.tab.gz")
    if not filen.is_file():
        msg = "Missing input file observations.tab.gz under [--dest]. " \
              "Check accession.py logs or specify a different directory. " \
              "Terminating program execution!"
        logger.critical(msg)
        return
    observations = pd.read_csv(filen, sep="\t")
    
    # output directory
    h5ad_dir = Path(args.dest, 'h5ad')
    h5ad_dir.mkdir(parents=True, exist_ok=False)
        
    msg = f"Processing {config['dataset_type']} {config['geo_accession']}. "
    logger.info(msg)
    
    # data wrangling specific to single-cell 
    if "single" in config['dataset_type'].lower():
        filen = Path(args.dest, f"{config['geo_accession']}_filelist.txt")
        if not filen.is_file():
            msg = f"Missing input file {config['geo_accession']}_filelist.txt under [--dest]. " \
                  "Terminating program execution!"
            logger.critical(msg)
            return
        
        filelist = pd.read_csv(
            filen,
            header=None,
            names=['File', 'Name', 'Time', 'Size', 'Type'],
            comment='#',
            sep='\t'
        )
        tarname = Path(args.dest, filelist[filelist.File=='Archive'].Name[0])
        filelist = filelist[~(filelist.File=='Archive')].copy()
        filelist.reset_index(drop=True, inplace=True)
        filelist['Sample'] = filelist['Name'].str.split('_', n=1, expand=True)[0].values
        
        # check which files to extract
        pattern = fr'.*({"|".join(args.pattern)}).*'
        filelist = filelist[filelist.Name.str.match(pattern, flags=re.IGNORECASE)]
        filelist.reset_index(drop=True, inplace=True)
        observations = observations[observations.geo_accession.isin(filelist.Sample)].copy()
        observations.reset_index(drop=True, inplace=True)
        
        # NOTE: we don't check observations (supplementary_file_) for consistency... 
        if args.input_format == 'MEX':
            mex_files = "|".join([args.mex_gene, args.mex_barcode, args.mex_matrix])
            pattern = fr'.*({mex_files}).*'
            filelist = filelist[filelist.Name.str.match(pattern, flags=re.IGNORECASE)]
            filelist.reset_index(drop=True, inplace=True)
            observations = observations[observations.geo_accession.isin(filelist.Sample)].copy()
            observations.reset_index(drop=True, inplace=True)
            
            parse_fmt = data_utils.parse_mex_fmt
            kwargs = {
                'mex_gene': args.mex_gene,
                'mex_barcode': args.mex_barcode,
                'mex_matrix': args.mex_matrix,
                'n_vars': args.mex_gene_ncols,
                'var_names': args.gene_var,
                'feature': 'gene_symbol',
                'rep_col': args.replicate_column,
                'clust_col': args.cluster_column
            }
        elif args.input_format == 'TXT':
            parse_fmt = data_utils.parse_txt_fmt
            kwargs = {
                'ext': args.txt_ext,
                'delimiter': args.txt_delimiter,
                'first_column_names': args.txt_first_column_names,
                'first_row_names': args.txt_first_row_names,
                'n_vars': 2 if args.txt_gene_cols == 'both' else 1,
                'var_names': args.gene_var,
                'feature': args.txt_gene_cols,
                'rep_col': args.replicate_column,
                'clust_col': args.cluster_column
            }
        else:
            msg = "[--input-format] undefined! Cannot process files." \
                  "Terminating program execution!"
            logger.critical(msg)
            return
                
        # extract
        extract_dir = Path(args.dest, 'extract')
        # make sure there are no filenames starting with "/"
        # if these are inconsistent, tarfile will throw an error...
        to_extract = [Path(n).name for n in filelist.Name]
        try:
            extract_dir.mkdir(parents=True, exist_ok=False)
            try:
                msg = f'Extracting {tarname}...'
                logger.info(msg)
                msg = f'... the following files will be extracted:'
                logger.info(msg)
                for sample, group in filelist.groupby('Sample'):
                    msg = f'{sample} {observations[observations.geo_accession==sample].title.values[0]}:\n' \
                        f'\t{", ".join(group.Name)}'
                    logger.info(msg)
                    
                tar = tarfile.open(tarname)
                
                def _members(tar, names):
                    return (m for m in tar.getmembers() if m.name in names)
            
                tar.extractall(members=_members(tar, to_extract), path=extract_dir)
                tar.close()
            except BaseException:
                raise
        except FileExistsError:
            msg = f'{extract_dir} already exists! Checking if files were extracted...'
            logging.warning(msg)
            all_files = [p.name for p in list(filter(Path.is_file, extract_dir.glob('*.*')))]
            if not (all(m in all_files for m in to_extract)):
                msg = f'... some (or all) files are missing! Existing files will not be overwritten. ' \
                    f'Terminating!'
                logger.critical(msg)
                return
            else:
                if len(all_files) > len(to_extract):
                    msg = f'There are more files ({len(all_files)}) under {extract_dir} ' \
                        f'than the number of files to extract ({len(to_extract)}). Terminating!'
                    logger.critical(msg)
                    return
                msg = f'... all files were already extracted! Continuing...'
                logger.warning(msg)
                
        # files are ready to be read in...
        # currently same processing whether data is raw or normalised (this shouln't 
        # really occur, given input formats)
        # if not args.normalised:
        # else: # data is already normalised
        #    pass
        
        h5ad_name = f"{args.name}_raw.h5ad"
        raw = parse_fmt(
            extract_dir,
            observations,
            filelist,
            **kwargs
        )
        adata = ad.concat(
            raw, 
            join='outer', 
            index_unique="-", 
            label='batch', 
            merge="same"
        )
            
        if args.do_normalise:
            msg = "[--do-normalise] is set, but this option has currently no " \
                  "effect on single-cell RNA-seq. Continuing!"
            logger.warning(msg)
        
        
    elif "bulk" in config['dataset_type'].lower():
        # force text format - pre-formatted input
        args.input_format = "TXT"
        h5ad_name = f"{args.name}_normalized.h5ad"
        if args.normalised:
            args.do_normalise = False
            msg = "Input treated as normalized data. Forcing [--do-normalise] to False." \
                  "Continuing!"
            logger.warning(msg)
        if not args.normalised and not args.do_normalise:
            h5ad_name = f"{args.name}_raw.h5ad"
            msg = "Input treated as un-normalized data, but [--do-normalise] is unset." \
                  "Continuing, but bulk RNA-seq data will NOT be normalized!"
            logger.warning(msg)
        all_files = [p.name for p in list(filter(Path.is_file, args.dest.glob('*.*')))]
        # any TXT format, compressed or not
        filenames = list(data_utils.glob_re(f"{config['geo_accession']}_read_count", all_files))
        if len(filenames) != 1:
            msg = "[--dest] must include a pre-formatted read count matrix  " \
                  f"with the following name: {config['geo_accession']}_read_count. " \
                  "Terminating program execution!"
            logger.critical(msg)
            return
        filelist = pd.DataFrame({"Sample": [config['geo_accession']], "Name": filenames})
        parse_fmt = data_utils.parse_txt_fmt
        kwargs = {
            'ext': args.txt_ext,
            'delimiter': args.txt_delimiter,
            'first_column_names': args.txt_first_column_names,
            'first_row_names': args.txt_first_row_names,
            'n_vars': 2 if args.txt_gene_cols == 'both' else 1,
            'var_names': args.gene_var,
            'feature': args.txt_gene_cols,
            'rep_col': args.replicate_column,
            'bulk': True
        }
        adata = parse_fmt(
            args.dest,
            observations,
            filelist,
            **kwargs
        )[config['geo_accession']]
        
        # bulk RNA-seq-specific normalization
        if args.do_normalise:
            msg = f"Normalizing using {args.normalization_method}."
            logger.info(msg)
            
            # defining R script and loading instance in Python
            r = robjects.r
            rootpath = Path(__file__).parent
            filepath = str(Path(rootpath, 'R-normalization.R'))
            r['source'](filepath)
            
            # get matrix
            X = adata.X.todense().T
            
            def convert(X, func, gene_lengths=None):
                # convert it into an R object to pass into the R function
                with localconverter(robjects.default_converter + numpy2ri.converter):
                    input_data_r = robjects.conversion.py2rpy(X)
                # invoke the R function and get the result
                if gene_lengths is None:
                    result_r = func(input_data_r)
                else:
                    with localconverter(robjects.default_converter + numpy2ri.converter):
                        gene_lengths_r = robjects.conversion.py2rpy(gene_lengths)
                    result_r = func(input_data_r, gene_lengths_r)
                # convert it back to a pandas dataframe (or numpy array)
                with localconverter(robjects.default_converter + numpy2ri.converter):
                    Y = robjects.conversion.rpy2py(result_r)
                return Y
                    
            normalisation_methods = {
                "TMM": robjects.globalenv['normalized_tmm_data'],
                "DESeq": robjects.globalenv['normalized_deseq2_data'],
                "TPM": robjects.globalenv['normalized_tpm_data'],
            }
            # check missing arguments for TPM
            if args.normalization_method == 'TPM':
                try:
                    gene_lengths = np.loadtxt(args.gene_lengths)
                except Exception as e:
                    msg = f"{e}\n Missing or cannot load [--gene-lengths] required " \
                        "for TPM normalization. Terminating!"
                    logger.critical(msg)
                X_normalised = convert(X, 
                                       normalisation_methods[args.normalization_method],
                                       gene_lengths=gene_lengths)
            else:
                X_normalised = convert(X, 
                                       normalisation_methods[args.normalization_method])
            
            X_normalised = csr_matrix(X_normalised.astype(np.float32)).T
            adata.layers["raw"] = adata.X
            adata.X = X_normalised
    else:
        msg = f"{config['dataset_type']} not supported!"
        logger.critical(msg)
        return
        
    # add missing features (id or symbol) to the final object
    if kwargs['n_vars'] == 1:
        qterm = "gene_id" if kwargs['feature'] == "gene_symbol" else "gene_symbol"
        scopes = args.mygene_scopes if qterm == "gene_symbol" else "symbol"
        adata = data_utils.add_missing_var_cols(adata,
                                                qterm,
                                                scopes,
                                                args.mygene_scopes,
                                                config['sample_taxid'],
                                                kwargs['var_names']
        )

    # finally, write to disk
    adata.write_h5ad(Path(h5ad_dir, h5ad_name))
    msg = f"Writing data to {Path(h5ad_dir, h5ad_name).as_posix()}. Completed!"
    logger.info(msg)
    

if __name__ == '__main__':
    main()
    
