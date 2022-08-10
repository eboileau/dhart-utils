#! /usr/bin/env python3

"""Data wrangling wrapper to prepare and reformat GEO 
   supplementary files into h5ad for input into DHART. 
"""

from ast import arg
import sys
import logging
import argparse
import re
from importlib_metadata import metadata
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import scipy as sc
import anndata as ad
import utils.utils as utils
import utils.data_utils as data_utils
import mygene

import tarfile
from pathlib import Path
from paths import PROJECT_ROOT
UTILS_PATH = PROJECT_ROOT / 'utils'

logger = logging.getLogger(__name__)


required_files = []

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Data wrangling wrapper.""")

    parser.add_argument('-d', '--dest', help="""Path to directory where GEO files reside.""",
                        type=Path, default=Path(__file__).parent)
    
    parser.add_argument('-n', '--name', help="""Name of H5AD output.""", type=str, 
                        default='adata')

    parser.add_argument('-b', '--bulk', help="""Specifies whether input file is bulk or single-cell data.""", 
                        action='store_true')

    parser.add_argument('--species', help="""Specifies the species to be used for missing gene ID or gene name info in bulk data. 
                        If both gene ID and name are provided, this is silently ignored. 
                        This is also silently ignored, if the -b flag isn't set. The default value is human. """, 
                        type=str, default='human')

    parser.add_argument('--gene-info', help="""Indicate presence of gene ID and/or gene name.""", type=str, 
                        default='ID', choices=['ID', 'name', 'ID+name'])

    parser.add_argument('--gene-info-id', help="""Provide input gene types as defined by mygene. 
                        The default value is ensembl.gene, hence ensembl is assumed as the standard input format.""", 
                        type=str, default='ensembl.gene')
    
    parser.add_argument('-fmt', '--input-format', help="""Input format: either MEX (e.g. TSV
                        and MTX), or TXT (e.g. TXT, TSV, or CSV). Market Exchange (MEX) format
                        is for gene-barcode matrix output from Cell Ranger; TXT format is 
                        for read count matrices (single file, genes x barcodes/wells).""", 
                        type=str, choices=['MEX', 'TXT']) # This may be reintroduced for scRNA processing only: , required=True
    
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

    parser.add_argument('--normalization-method', help="""Set the normalization method. 
                        Flag is silently ignored if --do-normalize isn't set""", 
                        default='edgeR-TMM', choices=['edgeR-TMM', 'deseq2', 'TPM'])
    
    utils.add_logging_options(parser)
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    utils.update_logging(args)
    
    msg = f'[wrangling]: {" ".join(sys.argv)}'
    logger.info(msg)
    
    # Set required files according to what data was provided
    if args.bulk: 
        required_files = ['bulk_input.csv', 'observations.tab.gz']
    else:
        required_files = ['_filelist.txt', '_RAW.tar', 'observations.tab.gz']
        
    # first check that all required files are available
    all_files = [p.name for p in list(filter(Path.is_file, args.dest.glob('*.*')))]
    pattern = fr'.*({"|".join(required_files)})'
    filenames = list(data_utils.glob_re(pattern, all_files))
    if len(filenames) != len(required_files):
        msg = f'Missing input files! Required files are: {", ".join(required_files)}. ' \
              f'Only these files were found: {", ".join(filenames)}'
        logger.critical(msg)
        return

    if not args.bulk:
        # This section applies to scRNA data

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
    else: 
        # This section applies to bulk RNA data

        # must be consistent with the order of required_files
        filen = [f for f in filenames if required_files[0] in f][0]
        bulkdata = pd.read_csv(
            Path(args.dest, filen),
            header=0,
            comment='#',
        )
        filen = [f for f in filenames if required_files[1] in f][0]
        observations = pd.read_csv(
            Path(args.dest, filen),
            sep='\t'
        )
        
        








        # Read data from bulk_input.csv

        # TODO: CSV reader for both input files (observations.tab.gz)

        if args.gene_info == 'ID':
            # Only the gene ID is present
            logger.debug('gene ID is present')

            gene_metadata = bulkdata.iloc[:,[0]]
            gene_metadata_series = bulkdata.iloc[:,0]
            sampledata = bulkdata.iloc[:,1:]
            print("input data shape: " + str(sampledata.shape))

            gene_format = args.gene_info_id

            mg = mygene.MyGeneInfo()
            gene_metadata_list = gene_metadata_series.tolist()
            species_from_args = args.species

            mygene_query_result = mg.querymany(gene_metadata_list, scopes=gene_format, fields='symbol,'+gene_format, species=species_from_args, as_dataframe=True)

            matching_table = mygene_query_result[[gene_format, "symbol"]]
            
            #df["Subjects"] = df["first_name"].map(Subjects)

            print(gene_metadata)

            gene_metadata["name"] = gene_metadata["ID"].map(matching_table[gene_format])

            print(mygene_query_result)
            print(gene_metadata)



        if args.gene_info == 'name': 
            # Only the gene name is present
            logger.debug('gene name is present')

            gene_metadata = bulkdata.iloc[:,[0]]
            gene_metadata_series = bulkdata.iloc[:,0]
            sampledata = bulkdata.iloc[:,1:]
            print("input data shape: " + str(sampledata.shape))

            mg = mygene.MyGeneInfo()
            gene_metadata_list = gene_metadata_series.tolist()
            species_from_args = args.species

            mygene_query_result = mg.querymany(gene_metadata_list, scopes='symbol', fields='symbol,ensembl.gene', species=species_from_args, as_dataframe=True)

            matching_table = mygene_query_result[["ensembl.gene", "symbol"]]
            
            #df["Subjects"] = df["first_name"].map(Subjects)

            gene_metadata["ID"] = gene_metadata["name"].map(matching_table["symbol"])

            print(mygene_query_result)
            print(gene_metadata)





        if args.gene_info == 'ID+name':
            # Both the gene ID and name are present
            logger.debug('gene ID and name are present')

            gene_metadata = bulkdata.iloc[:,0:2]
            sampledata = bulkdata.iloc[:,2:]
            print("input data shape: " + sampledata.shape)

        else: 
            # Unrecognized gene_info parameter
            msg = f'Unrecognized gene_info parameter. '
            logger.critical(msg)
        
        

        if args.do_normalise:
            # Data needs to be normalized, check for normalization method and normalize accordingly
            # Defining the R script and loading the instance in Python
            r = robjects.r
            rootpath = Path(__file__).parent
            filepath = str(Path(rootpath, 'R-normalization.R'))
            print(filepath)

            if args.normalization_method == 'edgeR-TMM':
                # edgeR's TMM normalization selected
                r['source'](filepath)
                
                # Load the function defined in R.
                edgeR_TMM = robjects.globalenv['normalized_tmm_data']

                # Convert it into an R object to pass into the R function
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    input_data_r = robjects.conversion.py2rpy(sampledata)
                
                #Invoke the R function and get the result
                df_result_r = edgeR_TMM(input_data_r)
                
                # Convert it back to a pandas dataframe.
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    df_result = robjects.conversion.rpy2py(df_result_r)
                
                logger.debug(df_result)
                print(df_result)
                print(df_result.shape)

            if args.normalization_method == 'deseq2':
                # deseq2 normalization selected

                # TODO: This is experimental, needs to be checked over

                # Defining the R script and loading the instance in Python
                r = robjects.r
                r['source'](filepath)

                # Load the function defined in R.
                deseq2 = robjects.globalenv['normalized_deseq2_data']

                # Convert it into an R object to pass into the R function
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    input_data_r = robjects.conversion.py2rpy(sampledata)

                #Invoke the R function and get the result
                df_result_r = deseq2(input_data_r)

                # Convert it back to a pandas dataframe.
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    df_result = robjects.conversion.rpy2py(df_result_r)

                logger.debug(df_result)
                print(df_result)
                print(df_result.shape)

            if args.normalization_method == 'TPM':
                # TPM normalization selected

                # Check if length information is present
                if 'length' in sampledata.columns:
                    lengthdata = sampledata["length"]
                else: 
                    logger.critical("Column named \"length\" must be present and must contain gene lenght information when using TPM. ")
                
                # Remove length data from sample for clean sample data
                sampledata = sampledata.drop('length', axis=1)

                r['source'](filepath)
                
                # Load the function defined in R.
                TPM = robjects.globalenv['normalized_tpm_data']

                # Convert it into an R object to pass into the R function
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    input_data_r = robjects.conversion.py2rpy(sampledata)
                
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    length_data_r = robjects.conversion.py2rpy(lengthdata)

                #Invoke the R function and get the result
                df_result_r = TPM(input_data_r, length_data_r)
                
                # Convert it back to a pandas dataframe.
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    df_result = robjects.conversion.rpy2py(df_result_r)

                logger.debug(df_result)
                print(df_result)
                print(df_result.shape)

            else: 
                # Unsupported normalization method selected
                msg = f'An unsupported normalization method was selected'
                logger.critical(msg)
                




            # TODO: Write normalized data to the H5AD file. 
            # TODO: Customize the next few lines of code to properly write the data


            kwargs = {
                'mex_gene': args.mex_gene,
                'mex_barcode': args.mex_barcode,
                'mex_matrix': args.mex_matrix,
                'n_vars': args.mex_gene_ncols,
                'var_names': args.mex_gene_var
            }


            # Write the h5ad file
            h5ad_dir = Path(args.dest, 'h5ad')
            h5ad_dir.mkdir(parents=True, exist_ok=False)
            #raw = parse_fmt(
            #        extract_dir,
            #        observations,
            #        filelist,
            #        **kwargs
            #    )
            adata = ad.concat(
                #raw, 
                join='outer', 
                index_unique="-", 
                label='batch', 
                merge="same")
            adata.write_h5ad(Path(h5ad_dir, f'{args.name}_raw.h5ad'))
            # Gene names go to adata.var_names and ids to a column "gene_id" in adata.var
            # Actually, I think we also need to add (redundantly) a "gene_symbol" column ( i.e. adata.var_names).
            # The actual matrix goes to adata.X (compressed to CSC format e.g. X = sp.sparse.csc_matrix(mat.T)). 
            # The sample names should match those of observations.tab.gz - title, and will be the adata.obs_names.
            # adata.obs can be populated from observations.tab.gz using minimally title, geo_accession, source_name_ch1, characteristics_ch1. 
            # The field characteristics may contain multiple entries e.g. characteristics_ch1.1.cell type characteristics_ch1.2.age, etc.



            # Remove comment: New Stuff
            X = sc.sparse.csc_matrix(df_result) # mat is your matrix of counts normalized, etc.
            adata = ad.AnnData(X)
            adata.var_names = gene_names # here features are e.g. gene names/symbols
            adata.var_names_make_unique('_')
            adata.var['gene_symbol'] = gene_names
            adata.var['gene_id'] = gene_ids


            # Now for the observations:
            adata.obs_names = obs_names # these are the sample names, they should match the column names from the count matrix, and I think the `title` from observations.tab.gz
            adata.obs = obs # where obs is a dataframe with observations from observations.tab.gz, where you selected the relevant entries

        else: 
            # Data needs to be normalized, alert the user to the missing flag
            msg = f'Data needs to be normalized when processing bulk data. Set the flag accordingly. '
            logger.critical(msg)
            return
    

if __name__ == '__main__':
    main()
    
