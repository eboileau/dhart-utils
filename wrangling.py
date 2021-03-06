#! /usr/bin/env python3

"""Data wrangling wrapper to prepare and reformat GEO 
   supplementary files into h5ad for input into DHART. 
"""

import sys
import logging
import argparse
import re

import tarfile
from pathlib import Path
from paths import PROJECT_ROOT
UTILS_PATH = PROJECT_ROOT / 'utils'

import pandas as pd
import anndata as ad
import utils.utils as utils
import utils.data_utils as data_utils

logger = logging.getLogger(__name__)


required_files = ['_filelist.txt', 
                  '_RAW.tar', 
                  'observations.tab.gz']


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Data wrangling wrapper.""")

    parser.add_argument('-d', '--dest', help="""Path to directory where GEO files reside.""",
                        type=Path, default=Path(__file__).parent)
    
    parser.add_argument('-n', '--name', help="""Name of H5AD output.""", type=str, 
                        default='adata')
    
    parser.add_argument('-fmt', '--input-format', help="""Input format: either MEX (e.g. TSV
                        and MTX), or TXT (e.g. TXT, TSV, or CSV). Market Exchange (MEX) format
                        is for gene-barcode matrix output from Cell Ranger; TXT format is 
                        for read count matrices (single file, genes x barcodes/wells).""", 
                        type=str, choices=['MEX', 'TXT'], required=True)
    
    parser.add_argument('--mex-gene', help="""File name pattern, typically "features" or
                        "genes", used with [--input-format MEX]. Silently ignored if
                        [--input-format TXT].""", default='feature', type=str)
    
    parser.add_argument('--mex-barcode', help="""File name pattern, typically "barcodes",
                        used with [--input-format MEX]. Silently ignored if
                        [--input-format TXT].""", default='barcode', type=str)
    
    parser.add_argument('--mex-matrix', help="""File name pattern, typically "matrix",
                        used with [--input-format MEX]. Silently ignored if
                        [--input-format TXT].""", default='matrix', type=str)

    parser.add_argument('--mex-gene-ncols', help="""Number of feature file columns, 
                        used with [--input-format MEX]. If [--mex-gene-ncols 1], it should
                        be "gene_symbol", otherwise they must be in the following 
                        order: ["gene_id", "gene_symbol", "feature_type"].
                        Silently ignored if [--input-format TXT].""", choices=[1, 2, 3], 
                        default=2, type=int)
                        
    parser.add_argument('--mex-gene-var', help="""Column to use as "var_names". If
                        [--mex-gene-ncols 1], uses whatever is there. Silently ignored if
                        [--input-format TXT].""", choices=['gene_symbol', 'gene_id'], 
                        default='gene_symbol', type=str)
    
    parser.add_argument('--txt-ext', help="""Extension that indicates the file type. 
                        If None, uses extension of filename, see the anndata read_txt 
                        function. Silently ignored if [--input-format MEX].""", 
                        default=None, type=str)
    
    parser.add_argument('--txt-delimiter', help="""If None, will split at arbitrary 
                        number of white spaces, see the anndata read_txt function.
                        Silently ignored if [--input-format MEX].""", default=None,
                        type=str)
    
    parser.add_argument('--txt-first-column-names', help="""Assume the first column 
                        stores row names. This is only necessary if these are not 
                        strings: strings in the first column are automatically assumed 
                        to be row names. See anndata read_txt function.
                        Silently ignored if [--input-format MEX].""", default=False,
                        type=bool)
    
    parser.add_argument('--txt-first-row-names', help="""Assume the first row stores 
                        column names (barcodes/wells). The behaviour of this flag is 
                        mostly dependent on the anndata read_txt function. If the 
                        first line starts with "#", this is ignored. If the first line 
                        does not start with "#", it might wrongly be read into X 
                        (and var), e.g. if column names are well numbers. If True, 
                        then column names are read separately and re-assigned, and this 
                        depends on [--txt-delimiter]. If False, only the first row is 
                        removed, and obs_names are numbered by default. Silently ignored 
                        if [--input-format MEX].""", default=False, type=bool)
    
    parser.add_argument('-p', '--pattern', help="""A space separated list of patterns 
                        used to glob files. Can be used with either [--input-format] options.
                        Search is performed independently (and before) file name pattern 
                        search for [--input-format MEX].""", default='', type=str, nargs='*')
    
    parser.add_argument('--normalised', help="""If this flag is present, then 
                        files are treated as normalised input.""", action='store_true')
    # TODO: if input is normalised, what do we do?
    
    parser.add_argument('--do-normalise', help="""If this flag is present, then 
                        output is also written as normalised, in addition to raw.
                        Input files must be raw.""", action='store_true')
    
    utils.add_logging_options(parser)
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    utils.update_logging(args)
    
    msg = f'[wrangling]: {" ".join(sys.argv)}'
    logger.info(msg)
    
    # first check that all required files are available
    all_files = [p.name for p in list(filter(Path.is_file, args.dest.glob('*.*')))]
    pattern = fr'.*({"|".join(required_files)})'
    filenames = list(data_utils.glob_re(pattern, all_files))
    if len(filenames) != len(required_files):
        msg = f'Missing input files! Required files are: {", ".join(required_files)}. ' \
              f'Only these files were found: {", ".join(filenames)}'
        logger.critical(msg)
        return
    # must be consistent with the order of required_files
    filen = [f for f in filenames if required_files[0] in f][0]
    filelist = pd.read_csv(
        Path(args.dest, filen),
        header=None,
        names=['File', 'Name', 'Time', 'Size', 'Type'],
        comment='#',
        sep='\t'
    )
    filen = [f for f in filenames if required_files[2] in f][0]
    observations = pd.read_csv(
        Path(args.dest, filen),
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
            'var_names': args.mex_gene_var
        }
    else: # TXT
        parse_fmt = data_utils.parse_txt_fmt
        kwargs = {
            'ext': args.txt_ext,
            'delimiter': args.txt_delimiter,
            'first_column_names': args.txt_first_column_names,
            'first_row_names': args.txt_first_row_names
        }
            
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
    h5ad_dir = Path(args.dest, 'h5ad')
    h5ad_dir.mkdir(parents=True, exist_ok=False)
    if not args.normalised:
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
            merge="same")
        adata.write_h5ad(Path(h5ad_dir, f'{args.name}_raw.h5ad'))
        
        if args.do_normalise:
            pass
        

    else: # data is already normalised
        pass
        
    
    

    
if __name__ == '__main__':
    main()
    
