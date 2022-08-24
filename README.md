
# DHART utils

This repository contains modules and scripts associated with the DHART (Digital HeART) portal, _e.g._ data retrieval from GEO, semi-automated
metadata preparation, data wrangling, _etc._

The DHART is based on the [gEAR framework](https://github.com//dieterich-lab/gEAR), developped by the [Institute for Genome Sciences](https://github.com/IGS/gEAR).

!!! info

    This repository is currently under development.

## Quickstart

### Docker
This software is supplied in a Docker image. Docker images enable the user to easily deploy software without any dependency tree issues. Since DHART utils isn't hosted in an image registry, the image needs to be built from the Dockerfile.
To do so, execute the following command from the root directory:

- accession: `docker build -t accession -f ./docker/accession/Dockerfile .`
- wrangling: `docker build -t wrangling -f ./docker/wrangling/Dockerfile .`

Once the image is built, it can then be launched by running the following command:

- accession: `docker run accession`
- wrangling: `docker run wrangling`

### Dependencies

When using this library without docker, the following dependencies need to be installed on your system:

- R version 4.2.x
- Python 3.9.5
- pandas
- anndata
- rpy2
- mygene
- dependencies from requirements.txt and R-dependencies.R

Dependency files can be found in the folder [dependencies](./dependencies/).

### accession.py

This script is used to download GEO Series and will prepare the data for ingestion into DHART.

### wrangling.py

This is a data wrangling wrapper to prepare and reformat GEO supplementary files into h5ad for input into DHART.

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

When processing bulk data, the following flags are important to set in `wrangling.py`:

- `--do-normalise`

Please find data scheme options below. Data can only be read if it conforms to one of the provided data schemes. Using the `--normalization-method` flag, a normalization method can be specified. `TMM`, `DESeq` and `TPM` are all valid normalization methods, `TMM` being the default. This flag is optional. A bundle of `observations.tab.gz` and the bulk data in a file named `GSEXXXXXX_read_count` needs to be submitted for each dataset. `observations.tab.gz` is provided by when downloading data using `accession.py`.

When processing data using TPM normalization, gene lengths must be provided separately. See Data Schemes for details.

### Data Schemes

Data schmeme choices include the following:

- 'gene_id': Includes just the gene ID, this is later used to look up the gene symbol. Both fields are then written to the H5AD file.
- 'gene_symbol': Includes just the gene symbol.
- 'both': Both the gene symbol and gene ID are already present in the input scheme, both will be written to the H5AD file. 

Please conform input data to one of the following schemes (note that column header for extra columns must be empty):

'gene_id':
|    | <sample#1> | ... | <sample#n> |
|----|------------|-----|------------|
|    | 1.343      | ... | 1.472      |
|    | 2.872      | ... | 3.971      |
|    | 1.341      | ... | 1.349      |

'gene_symbol':
|      | <sample#1> | ... | <sample#n> |
|------|------------|-----|------------|
|      | 1.343      | ... | 1.472      |
|      | 2.872      | ... | 3.971      |
|      | 1.341      | ... | 1.349      |

'both':
|    | name | <sample#1> | ... | <sample#n> |
|----|------|------------|-----|------------|
|    |      |            |     |            |
|    |      | 1.343      | ... | 1.472      |
|    |      | 2.872      | ... | 3.971      |
|    |      | 1.341      | ... | 1.349      |

## Running the tests

Currently no test are implemented. We're working on adding integration/regression tests baed on some sample data.
Test data for bulk under /prj/DHART/DATA/test_data/RNAseq/GSEXXXXX1.

## Contributing


## Versioning

This project currently does NOT adheres to [Semantic Versioning](http://semver.org/). See [CHANGELOG](CHANGELOG.md).

## Notes

GEO classes and parsing functions in `geo_utils` are based on [GEOparse](https://github.com/guma44/GEOparse), which does
not seem to be maintained anymore. These were modified to deal with GEO SERIES and SAMPLES metadata only.
Interacting with GEO/SRA online (data download) is done with [pysradb](https://github.com/saketkc/pysradb).
