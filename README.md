
# DHART utils

This repository contains modules and scripts associated with the DHART (Digital HeART) portal, _e.g._ data retrieval from GEO, semi-automated
metadata preparation, data wrangling, _etc._

The DHART is based on the [gEAR framework](https://github.com//dieterich-lab/gEAR), developped by the [Institute for Genome Sciences](https://github.com/IGS/gEAR).

!!! info

    This repository is currently under development.


## Quickstart

### Dependencies

When using this library without docker, the following dependencies need to be installed on your system: 

- R version 4.2.x
- Python 3.9.5
- pandas
- anndata
- rpy2
- dependencies from requirements.txt and R-dependencies.R

Dependency files can be found in the folder [dependencies](./dependencies/).

### accession.py

This script is used to download GEO Series and will prepare the data for ingestion into DHART.

The following parameters are supported:

`-d, --dest`: Path to directory where files will be written. For each GEO identifier, a new directory will be created.

`-g, --geo`: A GEO (GSE) identifier, or a space separated list of identifiers. Use this option only if all datasets have the same [dataset_type] and [annotation_source]. Tags cannot be specified.

`-f, --gfile`: Path to a tab-delimited file with GEO (GSE) identifier, data type, annotation source, and tags, one entry per line, without header. Tags can be empty. If both [--geo] and [--file] are given, the latter is silently ignored.

`-t, --dtype`: Allowed data types. If type include a space, it must be passed with quotes.

- default: 'single-cell RNA-Seq'

choices:

- 'single-cell RNA-Seq'
- 'bulk RNA-Seq'
- 'microarray'
- 'ChIP-Seq'
- 'ATAC-Seq'

`-a, --annotation`: Annotation source.

- default='Ensembl'

choices:

- 'Ensembl'
- 'Genbank'

### wrangling.py

This is a data wrangling wrapper to prepare and reformat GEO supplementary files into h5ad for input into DHART.

The following parameters are supported:

`-d, --dest`: Path to directory where GEO files reside.

- default is the parent directory.

`-n, --name`: Name of H5AD output.

- default='adata'

`-b, --bulk`: Specifies whether input file is bulk or single-cell data.

`--species`: Specifies the species to be used for missing gene ID or gene name info in bulk data. If both gene ID and name are provided, this is silently ignored. This is also silently ignored, if the -b flag isn't set.

- default='human'

`--gene-info`: Indicate presence of gene ID and/or gene name.

- default='ID'

choices:

- 'ID'
- 'name'
- 'ID+name'

`-fmt, --input-format`: Input format: either MEX (e.g. TSV and MTX), or TXT (e.g. TXT, TSV, or CSV). Market Exchange (MEX) format is for gene-barcode matrix output from Cell Ranger; TXT format is for read count matrices (single file, genes x barcodes/wells).

choices:
- 'MEX'
- 'TXT'

`--mex-gene`: File name pattern, typically "features" or "genes", used with [--input-format MEX]. Silently ignored if [--input-format TXT].

- default='feature'

`--mex-barcode`: File name pattern, typically "barcodes", used with [--input-format MEX]. Silently ignored if [--input-format TXT].

- default='barcode'

`--mex-matrix`: File name pattern, typically "matrix", used with [--input-format MEX]. Silently ignored if [--input-format TXT].

- default='matrix'

`--mex-gene-ncols`: Number of feature file columns, used with [--input-format MEX]. If [--mex-gene-ncols 1], it should be "gene_symbol", otherwise they must be in the following order: ["gene_id", "gene_symbol", "feature_type"]. Silently ignored if [--input-format TXT].

- default=2

choices:

- 1
- 2
- 3

`--mex-gene-var`: Column to use as "var_names". If [--mex-gene-ncols 1], uses whatever is there. Silently ignored if [--input-format TXT].

- default='gene_symbol'

choices:

- 'gene_symbol'
- 'gene_id'

`--txt-ext`: Extension that indicates the file type. If None, uses extension of filename, see the anndata read_txt function. Silently ignored if [--input-format MEX].

- default=None

`--txt-delimiter`: If None, will split at arbitrary number of white spaces, see the anndata read_txt function.
Silently ignored if [--input-format MEX].

- default=None

`--txt-first-column-names`: Assume the first column stores row names. This is only necessary if these are not strings: strings in the first column are automatically assumed to be row names. See anndata read_txt function.
Silently ignored if [--input-format MEX].

- default=False

`--txt-first-row-names`: Assume the first row stores column names (barcodes/wells). The behaviour of this flag is mostly dependent on the anndata read_txt function. If the first line starts with "#", this is ignored. If the first line does not start with "#", it might wrongly be read into X (and var), e.g. if column names are well numbers. If True, then column names are read separately and re-assigned, and this depends on [--txt-delimiter]. If False, only the first row is removed, and obs_names are numbered by default. Silently ignored if [--input-format MEX].

- default=False

`-p, --pattern`: A space separated list of patterns used to glob files. Can be used with either [--input-format] options. Search is performed independently (and before) file name pattern search for [--input-format MEX].

- default=''

`--normalised`: If this flag is present, then files are treated as normalised input.
    # TODO: if input is normalised, what do we do?

`--do-normalise`: If this flag is present, then output is also written as normalised, in addition to raw. Input files must be raw.

`--normalization-method`: Set the normalization method. Flag is silently ignored if --do-normalize isn't set.

- default='edgeR-TMM'

choices:

- 'edgeR-TMM'
- 'deseq2'
- 'TPM'

## Installation

This repository is currently not *installable* ( *i.e.* it is not a package ).


```bash
# clone dhart-utils
git clone https://github.com/eboileau/dhart-utils.git
# install dependencies (see below) in a virtual environment
module load python3/3.9.5_deb10
python3 -m venv ~/.virtualenvs/dhart_utils
source ~/.virtualenvs/dhart_utils/bin/activate
```

### Dependencies

Pinned version of selected dependencies are listed in the _requirements.txt_ file for reproducible installation.

## RNA bulk data download and processing

`accession.py` can download bulk RNA data. When processing bulk data, the following flags are important to set in `wrangling.py`:

- `-b, --bulk`
- `--do-normalise`
- `--gene-info`

The `-b` or `--bulk` flag needs to be provided along with the `--do-normalise` flag. Using the `--gene-info` flag, the data scheme needs to be set. Please find data scheme options below. Data can only be read if it conforms to one of the provided data schemes. Using the `--normalization-method` flag, a normalization method can be specified. `edgeR-TMM`, `deseq2` and `TPM` are all valid normalization methods, `edgeR-TMM` being the default. This flag is optional. A bundle of `observations.tab.gz` and the bulk data in a file named `bulk_input.csv` needs to be submitted for each dataset. `observations.tab.gz` is provided by when downloading data using `accession.py`.

### Data Schemes

Data schmeme choices include the following:

- 'ID': Inlcudes just the gene ID, this is later used to look up the gene name. Both fields are then written to the H5AD file.
- 'name': Includes just the gene name.
- 'ID+name': Both the gene name and gene ID are already present in the input scheme, both will be written to the H5AD file. 

Please conform input data to one of the following schemes:

'ID':
| ID | <sample#1> | ... | <sample#n> |
|----|------------|-----|------------|
|    | 1.343      | ... | 1.472      |
|    | 2.872      | ... | 3.971      |
|    | 1.341      | ... | 1.349      |

'name':
| name | <sample#1> | ... | <sample#n> |
|------|------------|-----|------------|
|      | 1.343      | ... | 1.472      |
|      | 2.872      | ... | 3.971      |
|      | 1.341      | ... | 1.349      |

'ID+name':
| ID | name | <sample#1> | ... | <sample#n> |
|----|------|------------|-----|------------|
|    |      |            |     |            |
|    |      | 1.343      | ... | 1.472      |
|    |      | 2.872      | ... | 3.971      |
|    |      | 1.341      | ... | 1.349      |

## Running the tests


## Contributing


## Versioning

This project currently does NOT adheres to [Semantic Versioning](http://semver.org/). See [CHANGELOG](CHANGELOG.md).

## Notes

GEO classes and parsing functions in `geo_utils` are based on [GEOparse](https://github.com/guma44/GEOparse), which does
not seem to be maintained anymore. These were modified to deal with GEO SERIES and SAMPLES metadata only.
Interacting with GEO/SRA online (data download) is done with [pysradb](https://github.com/saketkc/pysradb).
