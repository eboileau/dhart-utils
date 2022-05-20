#! /usr/bin/env python3

"""Accession script to download GEO SERIES using geo_utils and
   prepare metadata for DHART. 
"""

import sys
import logging
import argparse

from pathlib import Path
from paths import PROJECT_ROOT
UTILS_PATH = PROJECT_ROOT / 'utils'

import pandas as pd
import utils.utils as utils
import utils.geo_utils as geo_utils

logger = logging.getLogger(__name__)


data_choices = ['single-cell RNA-Seq', 
                'bulk RNA-Seq', 
                'microarray', 
                'ChIP-Seq', 
                'ATAC-Seq']
annotation_choices = ['Ensembl', 
                      'Genbank']

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Accession script for geo_utils.""")

    parser.add_argument('-d', '--dest', help="""Path to directory where files will be written.
                        For each GEO identifier, a new directory will be created.""",
                        type=Path, default=Path(__file__).parent)

    parser.add_argument('-g', '--geo', help="""A GEO (GSE) identifier, or a space separated 
                        list of identifiers. Use this option only if all datasets have the
                        same [dataset_type] and [annotation_source]. Tags cannot be specified
                        with this option.""", nargs="+", type=str)
    
    parser.add_argument('-f', '--gfile', help="""Path to a tab-delimited file with 
                        GEO (GSE) identifier, data type, annotation source, and tags, 
                        one entry per line, without header. Tags can be empty. If both
                        [--geo] and [--file] are given, the latter is silently ignored.""",
                        type=Path)
    
    parser.add_argument('-t', '--dtype', help="""Allowed data types. If type include a space, 
                        it must be passed with quotes.""", type=str, choices=data_choices, 
                        default='single-cell RNA-Seq')
    
    parser.add_argument('-a', '--annotation', help="""Annotation source.""", type=str,
                        choices=annotation_choices, default='Ensembl')
    
    
    utils.add_logging_options(parser)
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    utils.update_logging(args)
    
    msg = f'[accession]: {" ".join(sys.argv)}'
    logger.info(msg)
    
    if args.geo:
        args.geo = [geo.upper() for geo in args.geo]
        msg = f'Using data_type = "{args.dtype}" and ' \
              f'annotation_source = "{args.annotation}" ' \
              f'for {", ".join(args.geo)}'
        logger.warning(msg)
    
        for geo in args.geo:
            geotype = geo[:3]
            if geotype != 'GSE':
                msg = f'GEO type "{geotype}" not supported. Skipping!'
                logger.warning(msg)
                continue
            gse = geo_utils.get_GEO(geo, dest=args.dest)
            gse.get_supp_files(dest=args.dest)
            gse.get_input_files(
                dest=args.dest,
                dataset_type=args.dtype,
                annotation_source=args.annotation
            )
            
    elif args.gfile:
        geo = pd.read_csv(
            args.gfile, 
            sep='\t',
            header=None,
            names=['geo_accession', 'dataset_type', 'annotation_source', 'tags']
        )
        geo['geo_accession'] = geo['geo_accession'].str.upper()
        geo = geo[geo.geo_accession.str.startswith('GSE')].copy()
        geo = geo[geo['dataset_type'].isin(data_choices)].copy()
        geo = geo[geo['annotation_source'].isin(annotation_choices)].copy()
        geo['tags'] = geo['tags'].fillna('')
        
        def _get_geo(row):
            gse = geo_utils.get_GEO(row.geo_accession, dest=args.dest)
            gse.get_supp_files(dest=args.dest)
            gse.get_input_files(
                dest=args.dest,
                dataset_type=row.dataset_type,
                annotation_source=row.annotation_source,
                tags=row.tags
            )
        geo.apply(_get_geo, axis=1)
        
    else:
        # only if another argument is passed, otherwise we call [--help]
        logger.critical('Input [--geo] or [--file] must be specified.')
    
    logger.info('Completed successfully.')

    
if __name__ == '__main__':
    main()
    
