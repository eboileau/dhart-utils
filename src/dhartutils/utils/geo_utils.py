""" Python library to query and retrieve data from Gene Expression Omnibus (GEO)
    and interact with https://github.com/saketkc/pysradb.
    
    Classes and functions are based on https://github.com/guma44/GEOparse.
    
    * Download GEO SOFT file for GSE (not GSM, GPL, GDS)
    * Download associated SRA metadata
    * Download GEO supplementary data
    * Prepare metadata template for DHART (json) and observations file
"""


import re
import abc
import gzip
import json
import logging

import numpy as np

from pathlib import Path, PurePath
from typing import Union, Optional, List, Literal, Dict
from collections import defaultdict
from itertools import groupby
from pandas import DataFrame, read_csv
from contextlib import contextmanager
from six import StringIO, iteritems

from pysradb.geoweb import GEOweb
from pysradb import SRAweb

logger = logging.getLogger(__name__)


########################################################################
# Classes from pysradb to interact with GEO/SRA online
########################################################################


geoweb = GEOweb()
sraweb = SRAweb()
                
# 20.05.22
# sraweb.close() is a dummy method
# and 'GEOweb' object has no attribute 'db'


########################################################################
# Classes to store GEO metadata
########################################################################


class DataIncompatibilityException(Exception):
    pass


class NoMetadataException(Exception):
    pass


class BaseGEO(object):
    __metaclass__ = abc.ABCMeta

    geotype = None
    
    def __init__(
        self, 
        name: str = None, 
        metadata: Dict[str, List[str]] = None
    ):
        
        """\
        Initialize base GEO object.

        Parameters
        ----------
        name
            Name of the object.
        metadata
            Metadata information.

        Raises
        ------
        TypeError: Metadata should be a dict.
        """

        if not isinstance(metadata, dict):
            raise TypeError(
                "Metadata should be a dictionary not a %s" % str(type(metadata))
            )

        self.name = name
        self.metadata = metadata
        

    def get_metadata_attribute(self, metaname: str):
        
        """\
        Get metadata attribute by name.

        Parameters
        ----------
        metaname
            Name of the attribute

        Returns
        -------
        :obj:`list` or :obj:`str`: Value(s) of the requested metadata attribute

        Raises
        ------
        NoMetadataException: Attribute error
        TypeError: Metadata should be a list
        """
        
        metadata_value = self.metadata.get(metaname, None)
        if metadata_value is None:
            raise NoMetadataException(f'{metaname} is not a metadata attribute.')
        if not isinstance(metadata_value, list):
            raise TypeError('Metadata is not a list!')

        if len(metadata_value) > 1:
            return metadata_value
        else:
            return metadata_value[0]


    def get_accession(self):
        
        """\
        Return accession ID of the sample.

        Returns
        -------
        :obj:`str`: GEO accession ID
        """
        
        return self.get_metadata_attribute('geo_accession')


class SimpleGEO(BaseGEO):
    __metaclass__ = abc.ABCMeta

    def __init__(
        self, 
        name, 
        metadata, 
        table, 
        columns
    ):
        
        """\
        Initialize simple GEO object.

        Parameters
        ----------
        name
            Name of the object
        metadata
            Metadata information
        table (:obj:`pandas.DataFrame`)
            Table with the data from SOFT file. May be absent.
        columns (:obj:`pandas.DataFrame`)
            Description of the columns, number of columns, order and names represented as index
            in this DataFrame has to be the same as table.columns. May be absent.

        Raises
        ------
        :obj:`ValueError`: Table should be a DataFrame
        :obj:`ValueError`: Columns' description should be a DataFrame
        :obj:`DataIncompatibilityException`: Columns are wrong
        :obj:`ValueError`: Description has to be present in columns
        """
        
        if not isinstance(table, DataFrame):
            raise ValueError(f'Table data should be an instance ' \
                             f'of "pandas.DataFrame" not {str(type(table))}')
        if not isinstance(columns, DataFrame):
            raise ValueError(f'Columns description should be an instance ' \
                             f'of "pandas.DataFrame" not {str(type(table))}')

        BaseGEO.__init__(self, name=name, metadata=metadata)

        self.table = table
        self.columns = columns
        columns_are_correct = False
        if self.columns.index.tolist() != self.table.columns.tolist():
            if not self.columns.index.is_unique:
                # try to correct duplicate index in the same way pandas is
                # doing the columns
                logger.warning(f'Detected duplicated columns ' \
                               f'in {self.geotype}, {self.name}. Correcting.\n')
                indices = {}
                new_index = []
                for idx in self.columns.index:
                    if idx not in indices:
                        indices[idx] = 0
                        new_index.append(idx)
                    else:
                        indices[idx] += 1
                        new_index.append("%s.%i" % (idx, indices[idx]))
                self.columns.index = new_index
                if self.columns.index.tolist() == self.table.columns.tolist():
                    columns_are_correct = True
            if not columns_are_correct:
                # if the columns are still not correct check the order.
                if sorted(self.columns.index.tolist()) == sorted(
                    self.table.columns.tolist()
                ):
                    logger.warning(f'Data columns in d{self.geotype}, {self.name} ' \
                                   f'are not in order. Reordering.\n')
                    self.columns = self.columns.loc[self.table.columns]
                    if self.columns.index.tolist() == self.table.columns.tolist():
                        columns_are_correct = True
        else:
            columns_are_correct = True

        if not columns_are_correct:
            rows_in_columns = ", ".join(self.columns.index.tolist())
            columns_in_table = ", ".join(self.table.columns.tolist())
            raise DataIncompatibilityException(
                f'Data columns do not match columns description index in {self.name}\n' \
                f'Columns in table are: {columns_in_table}\n' \
                f'Index in columns are: {rows_in_columns}\n'
            )
        if self.columns.columns[0] != 'description':
            raise ValueError(f'Columns table must contain a column named "description".' \
                             f'Here columns are: {", ".join(map(str, self.columns.columns))}')


class GSM(SimpleGEO):
    
    """Class that represents a sample from the GEO database."""

    geotype = "SAMPLE"


class GSE(BaseGEO):
    
    
    """Class representing a GEO series"""

    geotype = "SERIES"

    def __init__(
        self, 
        name, 
        gsms=None,
        metadata=None, 
        platform_metadata=None,
        database_metadata=None
    ):
        
        """Initialize GSE.

        Parameters
        ----------
        name
            Name of the object.
        gsms
            A dictionary with SAMPLE information (GSM).
        metadata
            SERIES information.
        platform_metadata
            PLATFORM information.
        database_metadata
            DATABASE information. 
        """
        
        gsms = {} if gsms is None else gsms
        if not isinstance(gsms, dict):
            raise ValueError(f'GSMs should be a dictionary not a {str(type(gsms))}')

        for gsm_name, gsm in iteritems(gsms):
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"

        BaseGEO.__init__(self, name=name, metadata=metadata)

        self.gsms = gsms
        self.platform_metadata = platform_metadata
        self.database_metadata = database_metadata
        self.observations = None

    @property
    def get_observations(self):
        
        """Get the phenotype data for each sample."""
        
        if self.observations is None:
            pheno_data = {}
            for gsm_name, gsm in iteritems(self.gsms):
                tmp = {}
                for key, value in iteritems(gsm.metadata):
                    if len(value) == 0:
                        tmp[key] = np.nan
                    elif key.startswith("characteristics_"):
                        for i, char in enumerate(value):
                            char = re.split(r":\s+", char)
                            char_type, char_value = [char[0], ": ".join(char[1:])]
                            tmp[key + "." + str(i) + "." + char_type] = char_value
                    else:
                        tmp[key] = ",".join(value)
                pheno_data[gsm_name] = tmp
            self.observations = DataFrame(pheno_data).T  
        return self.observations

   
    def get_supp_files(
        self,
        dest: Union[Path, str] = '.',
        download_sra_meta: bool = True
    ):
       
        """\
        Download supplementary data from GEO/SRA using pysradb.


        Parameters
        ----------
        dest
            Directory to download the data
        download_sra_meta 
            Download SRA metadata 
        """
        
        geo = self.get_accession()
        links, root_url = geoweb.get_download_links(geo)
        geoweb.download(
            links=links, 
            root_url=root_url, 
            gse=geo, 
            out_dir=dest
        )
        
        out_dir = Path(dest, geo, 'SraRunTable.txt')
        logger.info(f'Writing {out_dir}')
        try:
            SRP = sraweb.gse_to_srp(geo).study_accession.values[0]
            metadata = sraweb.sra_metadata(SRP)
            metadata.to_csv(
                out_dir, 
                index=False, 
                header=True
            )
        except:
            # we had some failed fetch BioProject... 
            msg = f'Failed to fetch (and download) SRA metadata for {geo}...' 
            logger.warning(msg)
        
    
    def get_input_files(
        self,
        dest: Union[Path, str] = '.',
        dataset_type: Literal['single-cell RNA-Seq', 
                              'bulk RNA-Seq', 
                              'microarray', 
                              'ChIP-Seq',
                              'ATAC-Seq'] = 'single-cell RNA-Seq',
        genome_assembly: str = None,
        annotation_source: str = None,
        annotation_release: int = None,
        tags: Union[str, List[str]] = None
    ):
       
        """\
        Write observations (data wrangling) and template (metadata) to disk.


        Parameters
        ----------
        dest
            Directory to download the data.
        dataset_type
            Allowed data types.
        genome_assembly
            Genome assembly. If None, try to guess from the SOFT file.
        annotation_source
            Source of annotation. If None, leave out.
        annotation_release
            Annotation release. If None, leave out.
        tags
            Comma-separated string or list of keywords.
        """
        
        geo = self.get_accession()
        
        out_dir = Path(dest, geo, 'observations.tab.gz')
        logger.info(f'Writing {out_dir}')
        metadata = self.get_observations
        metadata.to_csv(
            out_dir, 
            index=False, 
            header=True,
            sep='\t',
            compression='gzip'
        )
        
        # try to guess the assembly using first SAMPLE
        try:
            sample = self.observations.geo_accession.iloc[1]
            data_processing = self.gsms[sample].get_metadata_attribute('data_processing')
            genome_assembly = [l.split(': ')[1] for l in data_processing if l.split(': ')[0] == 'Genome_build'][0]
        except:
            genome_assembly = ''
        
        if isinstance(tags, list):
            tags = ','.join(tags)
        elif isinstance(tags, str):
            tags = ','.join(tags.split(','))
        else:
            tags = ''
            
        # Only fill essential fields, remaining fields are filled at upload
        
        # in some instances, contact email is missing...
        try:
            contact_email = self.get_metadata_attribute('contact_email')
        except NoMetadataException:
            msg = f'No contact email found for {geo}! Replacing by default admin@localhost.'
            logger.warning(msg)
            contact_email = 'admin@localhost'
        
        metadata = {
            'title': self.get_metadata_attribute('title'),
            'summary': self.get_metadata_attribute('summary'),
            'dataset_type': dataset_type,
            'assembly': genome_assembly,
            'annotation_source': annotation_source if annotation_source is not None else '',
            'annotation_release_number': annotation_release if annotation_release is not None else '',
            'geo_accession': self.get_accession(),
            'contact_email': contact_email,
            'contact_institute': self.get_metadata_attribute('contact_institute'),
            'contact_name': self.get_metadata_attribute('contact_name'),
            'sample_taxid': self.get_metadata_attribute('sample_taxid'), # check if consistent with platform_taxid?
            'tags': tags
        }
        # in particular for summary, if there are multiple "!Series_summary" entries...
        if type(metadata['summary']) == list:
            metadata['summary'] = " ".join(metadata['summary'])
        out_dir = Path(dest, geo, f'{geo}.json')
        logger.info(f'Writing {out_dir}')
        with open(out_dir, "w") as out_file:
            json.dump(metadata, out_file)
        

        
########################################################################
# GEO utils
########################################################################

        
def get_GEO(
        geo: Optional[str] = None,
        filepath: Optional[Union[Path, str]] = None,
        dest: Union[Path, str] = "."
    ):
    
    """\
    Get GEO (SERIES) entry, from GEO or from SOFT file.
    
    Parameters
    ----------
    geo 
        GEO identifier.
    filepath
        Path to local SOFT file. If both geo and filepath 
        are given, the latter is ignored. The SOFT file name is
        of the form e.g. "GSEnnnnnn_family.soft.gz". 
    dest
        Directory to download data.
        
    Returns
    -------
    :obj:`BaseGEO`: A GEO object.
    """
    
    if geo is not None:
        if filepath is not None:
            logger.info(f'Using GEO identifier {geo}. '\
                        f'Ignoring {filepath}.')
        filepath = _get_file(
            geo,
            dest=Path(dest)
        )
    else:
        if filepath is None:
            raise Exception('GEO identifier or SOFT file missing!')
        filepath = Path(filepath)
            
    logger.info(f'Parsing {filepath}')
    return parse_GSE(filepath)
    
    
def _get_file(geo, dest):
    
    """\
    Download SOFT file for a given GEO accession.
    
    Parameters
    ----------
    geo 
        GEO identifier.
    dest
        Directory to download data.
        
    Returns
    -------
    :obj:`str`: Path to downloaded file.
    """
    
    geo = geo.upper()
    range_subdir = re.sub(r'\d{1,3}$', 'nnn', geo)
    root_url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/' \
               f'{range_subdir}/{geo}/soft/'
    links = [f'{geo}_family.soft.gz']
    out_dir = dest # Path(dest, geo).mkdir(parents=True, exist_ok=True) 

    geoweb.download(
        links=links, 
        root_url=root_url, 
        gse=geo, 
        out_dir=out_dir
    )
    
    return Path(out_dir, geo, links[0])


@contextmanager
def _sopen(filepath):
    
    """I/O context manager
    
    Parameters
    ----------
    filepath 
        Path to the file.
        
    Returns
    -------
    Context manager for file handle.
    """
    
    kwargs = {}
    if filepath.suffix == ".gz":
        kwargs["mode"] = "rt"
        fopen = gzip.open
    else:
        kwargs["mode"] = "r"
        fopen = open
    fh = fopen(filepath, **kwargs)
    
    try:
        yield fh
    except IOError:
        fh.close()
    finally:
        fh.close()


def __parse_entry(entry_line: str):
    
    """\
    Parse the SOFT file entry name line that starts with '^', '!' or '#'.
    
    Parameters
    ----------
    entry_line
        Line to be parsed.
        
    Returns
    -------
    :obj:`2-tuple`: Type of entry, value of entry.
    """
    
    if entry_line.startswith("!"):
        entry_line = re.sub(r"!\w*?_", "", entry_line)
    else:
        entry_line = entry_line.strip()[1:]
    try:
        entry_type, entry_name = [i.strip() for i in entry_line.split("=", 1)]
    except ValueError:
        entry_type = [i.strip() for i in entry_line.split("=", 1)][0]
        entry_name = ""
    return entry_type, entry_name


def _parse_entry_name(nameline: str):
    
    """\
    Parse line that starts with ^ and assign the name to it.
    
    Parameters
    ----------
    nameline
        Line to parse.
        
    Returns
    -------
    :obj:`str`: Entry name.
    """
    
    entry_type, entry_name = __parse_entry(nameline)
    return entry_name


def _parse_metadata(lines: List[str]):
    
    """\
    Parse list of lines with metadata information from SOFT file.
    
    Parameters
    ----------
    lines
        List of lines.
    
    Returns
    -------
    :obj:`dict`: Metadata from SOFT file.
    """
    
    meta = defaultdict(list)
    for line in lines:
        line = line.rstrip()
        if line.startswith("!"):
            if "_table_begin" in line or "_table_end" in line:
                continue
            key, value = __parse_entry(line)
            meta[key].append(value)

    return dict(meta)


def _parse_columns(lines: List[str]):
    
    """\
    Parse list of lines with columns description from SOFT file.
    
    Parameters
    ----------
    lines 
        List of lines.
    
    Returns
    -------
    :obj:`pandas.DataFrame`: Columns description.
    """
    
    data = []
    index = []
    for line in lines:
        line = line.rstrip()
        if line.startswith("#"):
            tmp = __parse_entry(line)
            data.append(tmp[1])
            index.append(tmp[0])

    return DataFrame(data, index=index, columns=["description"])


def _parse_table_data(lines):
    
    """\
    Parse list of lines from SOFT file into DataFrame.
    
    Parameters
    ----------
    lines
        Iterator over the lines.
    
    Returns
    -------
    :obj:`pandas.DataFrame`: Table data.
    """
    
    # filter lines that do not start with symbols
    data = "\n".join(
        [i.rstrip() for i in lines if not i.startswith(("^", "!", "#")) and i.rstrip()]
    )
    if data:
        return read_csv(StringIO(data), index_col=None, sep="\t")
    else:
        return DataFrame()
    


def parse_GSM(
        filepath: Union[Path, str, List[str]], 
        entry_name: str = None):
    
    """\
    Parse GSM entry from SOFT file.
    
    Parameters
    ----------
    filepath
        Path to file with 1 GSM entry or list of lines representing GSM from GSE file.
    entry_name
        Name of the entry, inferred from the data.

    Returns
    -------
    :obj:`GEOparse.GSM`: A GSM object.
    """
    
    if isinstance(filepath, str) or isinstance(filepath, PurePath):
        with _sopen(filepath) as f:
            soft = []
            has_table = False
            for line in f:
                if "_table_begin" in line or (not line.startswith(("^", "!", "#"))):
                    has_table = True
                soft.append(line.rstrip())
    else:
        soft = []
        has_table = False
        for line in filepath:
            if "_table_begin" in line or (not line.startswith(("^", "!", "#"))):
                has_table = True
            soft.append(line.rstrip())

    if entry_name is None:
        sets = [i for i in soft if i.startswith("^")]
        if len(sets) == 0:
            raise NoEntriesException('No entries found!')
        entry_name = _parse_entry_name(sets[0])

    columns = _parse_columns(soft)
    metadata = _parse_metadata(soft)
    if has_table:
        table_data = _parse_table_data(soft)
    else:
        table_data = DataFrame()

    gsm = GSM(name=entry_name, 
              table=table_data, 
              columns=columns, 
              metadata=metadata)

    return gsm


def parse_GSE(filepath: Union[Path, str]):
    
    """\
    Parse GSE SOFT file.
    
    Parameters
    ----------
    filepath
        Path to GSE SOFT file.

    Returns
    -------
    :obj:`GSE`: A GSE object.
    """

    gsms = {}
    metadata = {}
    gse_name = None
    platform_metadata = {}
    database_metadata = None
    
    with _sopen(filepath) as soft:
        groupper = groupby(soft, lambda x: x.startswith("^"))
        for is_new_entry, group in groupper:
            if is_new_entry:
                entry_type, entry_name = __parse_entry(next(group))
                logger.debug(f'{entry_type.upper()}, {entry_name}')
                if entry_type == "SERIES":
                    gse_name = entry_name
                    is_data, data_group = next(groupper)
                    metadata = _parse_metadata(data_group)
                elif entry_type == "SAMPLE":
                    is_data, data_group = next(groupper)
                    gsms[entry_name] = parse_GSM(data_group, entry_name)
                elif entry_type == "PLATFORM":
                    is_data, data_group = next(groupper)
                    platform_metadata[entry_name] = _parse_metadata(data_group)
                elif entry_type == "DATABASE":
                    is_data, data_group = next(groupper)
                    database_metadata = _parse_metadata(data_group)
                else:
                    logger.error(f'Cannot recognize type {entry_type}')
    
    gse = GSE(name=gse_name, 
              gsms=gsms, 
              metadata=metadata, 
              platform_metadata=platform_metadata, 
              database_metadata=database_metadata)
    
    return gse

