""" Python library to help in data preparation
"""

import re
import gzip
import logging

from pathlib import Path
from typing import Union, Optional, List, Literal, Dict

from anndata import AnnData
from anndata.utils import make_index_unique
from scipy.io import mmread
from scipy.sparse import csr_matrix, issparse
from pandas import DataFrame

logger = logging.getLogger(__name__)


def glob_re(pattern, strings):
    return filter(re.compile(pattern).match, strings)
    
    
def read_mex(    
    path: Union[Path, str],
    prefix: List[Union[Path, str]],
    var_names: Literal['gene_symbol', 'gene_id'] = 'gene_symbol',
    n_vars: int = 2,
    make_unique: bool = True,
    dtype: str = "float32"
) -> AnnData:
    
    """\
    Read 10x-Genomics-formatted mtx files. A simplified/modified
    version of scanpy read_10x_mtx.
    
    Parameters
    ----------
    path
        Path to directory where files reside.
    prefix
        List of MEX files (typically feature, barcode, matrix).
    var_names
        The variables index.
    n_vars
        Number of columns for feature file.
    make_unique
        Whether to make the variables index unique by appending '-1',
        '-2' etc. or not.
    dtype
        Numpy data type.
        
    Returns
    -------
    :class:`~AnnData` object
    """
    
    import pandas as pd
    
    genefile = Path(path, prefix[0])
    barcodefile = Path(path, prefix[1])
    matrixfile = Path(path, prefix[2])
    
    mat = mmread(matrixfile)
    mat = csr_matrix(mat).T
    adata = AnnData(mat, dtype=dtype)
    
    cols = ['gene_id', 'gene_symbol', 'feature_type']
    col = cols[0]
    if n_vars == 1:
        # should be symbols to work with the portal, but currently we don't check...
        cols = [cols[1]]
        var_names = cols[0]
    else:
        cols = cols[:n_vars]
        if var_names == col:
            col = 'gene_symbol'
            
    genes = pd.read_csv(
        genefile, 
        header=None, 
        names=cols,
        usecols=cols,
        sep='\t')
    var_names = genes[var_names].values
    if make_unique:
        var_names = make_index_unique(pd.Index(var_names), join='_')
    adata.var_names = var_names
    if n_vars > 1:
        adata.var[col] = genes[col].values
        if col != 'gene_symbol':
            adata.var['gene_symbol'] = var_names
    if n_vars > 2:
        adata.var['feature_type'] = genes['feature_type'].values
        
    adata.obs_names = pd.read_csv(
        barcodefile, 
        header=None)[0].values
    
    return adata


def read_txt(    
    path: Union[Path, str],
    ext: Optional[str] = None,
    delimiter: Optional[str] = None,
    first_column_names: bool = False,
    first_row_names: bool = False,
    var_names: Literal['gene_symbol', 'gene_id'] = 'gene_symbol',
    n_vars: int = 1
) -> AnnData:
    
    """\
    Read `.txt`, `.tab`, `.data` (text) file.. A modified
    version of scanpy read_text (:func:`~anndata.read_csv` but 
    with default delimiter `None`).
    
    Parameters
    ----------
    path
        Path to input file.
    ext
        Extension that indicates the file type. If ``None``, uses extension of
        filename.
    delimiter
        Delimiter that separates data within text file. If ``None``, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at any single white space ``' '``. Note that this is also passed
        to pandas.read_csv if n_vars == 2.
    first_column_names
        Assume the first column stores row names. This is only necessary if
        these are not strings: strings in the first column are automatically
        assumed to be row names.
    first_row_names
        Assume the first row stores column names. The behaviour of this flag is 
        mostly dependent on the anndata read_txt function. The first line might 
        wrongly be read into X (and var), e.g. if column names are well 
        numbers. If True, then column names are read separately and re-assigned, 
        and this depends on delimiter. If False, fall back entirely to read_txt.
    var_names
        The variables index.
    n_vars
        Number of columns (features). 
        
    Returns
    -------
    :class:`~AnnData` object
    """
    
    import scanpy as sc 
    import pandas as pd
    
    path = Path(path)
    
    # TXT file contains 2 columns for gene ids (1st col.) 
    # and symbols (2nd col.)
    if n_vars == 2:
        cols = ['gene_id', 'gene_symbol']
        idx = 0
        if var_names == cols[0]:
            idx = 1
        col = cols[idx]
        # if sep is None, fallback to Python parsing engine
        mat = pd.read_csv(path, sep=delimiter)
        var_col = mat.iloc[:,idx].values
        mat.drop(mat.columns[idx], axis=1, inplace=True)
        
        path = Path(str(path).removesuffix("".join(path.suffixes)))
        path = Path(path.parent, f"{path.stem}_tmp.csv.gz")
        mat.to_csv(path, sep=',', compression="gzip", index=False)

    if path.suffix == ".gz":
        with gzip.open(path, mode="rt") as f:
            first_line = f.readline().strip()
    else:
        with path.open() as f:
            first_line = f.readline().strip()
    
    adata = sc.read(
        path,
        ext=ext,
        delimiter=delimiter,
        first_column_names=first_column_names).T

    if first_row_names:
        adata = adata[:,1:].copy()
        # assume first entry is actually the first column name (e.g. gene_ids)...
        obs_names = first_line.split(delimiter)[1:]
        if len(obs_names) == len(adata.obs_names):
            adata.obs_names = obs_names
        else:
            msg = f'Cannot re-assign column names as ' \
                f'len(adata.obs_names) = {len(adata.obs_names)} != ' \
                f'len(obs_names) = {len(obs_names)}. Skipping!'
            logger.warning(msg)
            
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
        
    # adjust var_names
    if n_vars == 2:
        adata.var[col] = var_col
        if col != 'gene_symbol':
            adata.var['gene_symbol'] = adata.var_names
        path.unlink() # remove tmp file

    return adata


def add_obs(
    geo: str,
    obs: DataFrame,
    adata: AnnData,
    bulk: bool = False
) -> AnnData:
    
    pretty_cols = dict()
    cols = [c for c in obs.columns if c.startswith('source_name_')]
    if cols:
        pretty_cols[cols[0]] = 'source_name' # there should only be one such entry...
    cols.extend([c for c in obs.columns if c.startswith('characteristics_')])
    pretty_cols = dict(pretty_cols, **{c:c.rsplit('.', 1)[1] for c in cols if not c.startswith('source_name_')})
    
    if not cols:
        return adata
    if bulk:
        obs.rename(columns=pretty_cols, inplace=True)
        cols = ['title', 'geo_accession']
        cols.extend(list(pretty_cols.values()))
        obs = obs[cols]
        obs_names = adata.obs_names
        adata.obs = obs
        adata.obs_names = obs_names
        return adata
    
    info = obs[obs.geo_accession==geo]
    for c in cols:
        adata.obs[pretty_cols[c]] = info[c].values[0]
    
    return adata


def parse_mex_fmt(
    path: Union[Path, str],
    observations: DataFrame,
    files: DataFrame,
    **kwargs
) -> Dict[str, AnnData]:
    
    file_order = [kwargs['mex_gene'], kwargs['mex_barcode'], kwargs['mex_matrix']]
             
    adatad = {}
    grouped = files.groupby('Sample')
    for geo, info in grouped:
        files = info.Name.values
        files = [list(glob_re(fr'.*({p})', files))[0] for p in file_order]
        adata = read_mex(
            path=path,
            prefix=files,
            var_names=kwargs['var_names'],
            n_vars=kwargs['n_vars']
        )
        adata = add_obs(geo, observations, adata)
        adatad[geo] = adata
    
    return adatad


def parse_txt_fmt(
    path: Union[Path, str],
    observations: DataFrame,
    files: DataFrame,
    **kwargs
) -> Dict[str, AnnData]:
    
    def _parse(
        sample,
        filename,
        observations,
        **kwargs
    ) -> AnnData:
        
        adata = read_txt(
            path=filename,
            ext=kwargs['ext'],
            delimiter=kwargs['delimiter'],
            first_column_names=kwargs['first_column_names'],
            first_row_names=kwargs['first_row_names'],
            var_names=kwargs['var_names'],
            n_vars=kwargs['n_vars']
        )
        bulk = kwargs.get("bulk", False)
        adata = add_obs(sample, observations, adata, bulk=bulk)
        return adata
      
    adatad = {sample: _parse(sample, Path(path, filename), observations, **kwargs) for 
                  sample, filename in zip(files['Sample'], files['Name'])} 

    return adatad


def add_missing_var_cols(
    adata: AnnData,
    qterm: str, 
    scopes: str,
    fields: str,
    species: int,
    var_name: str
) -> AnnData:
    
    import pandas as pd
    import mygene
    
    mg = mygene.MyGeneInfo()
    mg_query = mg.querymany(adata.var_names, 
                            scopes=scopes, 
                            fields='symbol,' + fields, 
                            species=species, 
                            as_dataframe=True)
    mg_query.reset_index(inplace=True)
    mg_query.drop_duplicates(subset=['query'], inplace=True)
        
    feature_present = "gene_symbol" if qterm == "gene_id" else "gene_id"
    col = "symbol"
    if feature_present == "gene_symbol":
        col = fields
    
    mg_query.rename(columns={"query": feature_present, col: qterm}, inplace=True)
    mg_query = mg_query[[feature_present, qterm]]
    adata.var[feature_present] = adata.var_names
    adata.var = pd.merge(adata.var, mg_query, on=feature_present, how='left')
    
    adata.var[qterm].fillna(adata.var[feature_present], inplace=True)
    adata.var_names = adata.var[var_name].values
    
    return adata

