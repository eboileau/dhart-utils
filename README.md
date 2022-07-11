
# DHART utils

This repository contains modules and scripts associated with the DHART (Digital HeART) portal, _e.g._ data retrieval from GEO, semi-automated
metadata preparation, data wrangling, _etc._

The DHART is based on the [gEAR framework](https://github.com//dieterich-lab/gEAR), developped by the [Institute for Genome Sciences](https://github.com/IGS/gEAR).


!!! info

    This repository is currently under development.


## Quickstart

There is currently no documentation available.

#### accession.py

This script is used to download GEO Series and will prepare the data for ingestion into DHART. 

The following parameters are supported: 

-d, --dest: Path to directory where files will be written. For each GEO identifier, a new directory will be created.

-g, --geo: A GEO (GSE) identifier, or a space separated list of identifiers. Use this option only if all datasets have the same [dataset_type] and [annotation_source]. Tags cannot be specified.

-f, --gfile: Path to a tab-delimited file with GEO (GSE) identifier, data type, annotation source, and tags, one entry per line, without header. Tags can be empty. If both [--geo] and [--file] are given, the latter is silently ignored.

-t, --dtype: Allowed data types. If type include a space, it must be passed with quotes. 

- default: 'single-cell RNA-Seq'

choices:
- 'single-cell RNA-Seq'
- 'bulk RNA-Seq'
- 'microarray'
- 'ChIP-Seq'
- 'ATAC-Seq'

-a, --annotation: Annotation source.

- default='Ensembl'

choices: 
- 'Ensembl'
- 'Genbank'

### Installation

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


## Running the tests


## Contributing


## Versioning

This project currently does NOT adheres to [Semantic Versioning](http://semver.org/). See [CHANGELOG](CHANGELOG.md).


## Notes

GEO classes and parsing functions in `geo_utils` are based on [GEOparse](https://github.com/guma44/GEOparse), which does
not seem to be maintained anymore. These were modified to deal with GEO SERIES and SAMPLES metadata only.
Interacting with GEO/SRA online (data download) is done with [pysradb](https://github.com/saketkc/pysradb).

