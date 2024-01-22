
# DHART-utils

This repository contains tools associated with the DHART (Digital HeART) portal, _e.g._ data retrieval from GEO, semi-automated
metadata preparation, data wrangling, _etc._

The DHART is based on the [gEAR framework](https://github.com//dieterich-lab/gEAR), developped by the [Institute for Genome Sciences](https://github.com/IGS/gEAR).


<a id="quickstart"></a>

## Quickstart

For Docker (Singularity) usage, see [Containers](#docker).

If installing the package, first see [Dependencies](#dependencies).
Create a virtual environment, clone the repository and install **dhart-utils**:

```bash
# create virtual environment (Python >=3.9)
python3 -m venv /path/to/virtual/environment

# activate the environment
source /path/to/virtual/environment/bin/activate

# clone the repository
git clone https://github.com/eboileau/dhart-utils.git
cd dhart-utils

# install
pip3 install --upgrade pip setuptools wheel
pip3 --verbose install . 2>&1 | tee install.log
```

### accession

This script is used to download GEO series, prepare the metadata for ingestion into DHART, and data files for wrangling.
Call the program with `--geo` for a single GEO identifier, or with `--gfile` for many identifiers *e.g.* [ACCESSION.txt](data/ACCESSION.txt).
To see the full list of parameters, call the program on the command line without any option: `accession`.

### wrangling

This script is a data wrangling wrapper to prepare and reformat GEO supplementary files (output of `accession`) into H5AD for input into the DHART portal.
The program currently only handles *single-cell RNA-Seq* and *bulk RNA-Seq*, and adheres to the [gEAR standard formatting guidelines](https://github.com/IGS/gEAR/blob/main/docs/Documentation/UploadingOverview.md). To enable some of the portal functionalities, use `--replicate-column` and `--cluster-column` if relevant. 

**Note** As DHART is being actively developed, formatting guidelines may change in the future!

The program handles two basic `--input-format`: *TXT (e.g. TXT, TSV, or CSV)* and *MEX*, and each has is own set of options. 
For the portal, single-cell RNA-Seq data is *NOT* normalized, but bulk RNA-Seq must be normalized (use `--do-normalise`).

For plate-based single-cell RNA-Seq ( *e.g.* with well numbers as observations), pass `--txt-first-row-names`, in particular if there is a 
header for gene ids/symbol column. In general, however, avoid using headers! To see the full list of parameters, call the program on the command line without any option: `wrangling`.


<a id="dependencies"></a>

## Dependencies

Unless you are using Docker (or Singularity), installing **dhart-utils** also requires a working R >4.1.x installation
with packages listed in [R-dependencies.r](dependencies/R-dependencies.r)

Pinned version of selected Python packages are listed in [requirements.txt](dependencies/requirements.txt) for reproducible installation.

**Note:** Depending on your setup, environment variables may need to be updated, *e.g.* `export LD_LIBRARY_PATH="/path/to/your/R/installation/lib/R/lib:$LD_LIBRARY_PATH"` after installing **rpy2**.

<a id="docker"></a>

## Containers

### Docker 

Docker images allow to deploy the package resources without any dependency tree issues.
However, **dhart-utils** images are not yet published (hosted in an image registry) or integrated to GitHub actions, and thus
needs to be built from the Dockerfile. 

Assuming that Docker is installed, execute the following command from the root directory:

- accession: `docker build --rm -t accession:latest -f docker/accession/Dockerfile .`
- wrangling: `docker build --rm -t wrangling:latest -f docker/wrangling/Dockerfile .`

where the tag *latest* can be changed.

Once the image is built, the Docker `run` command can be used to create a writeable container and start it:

- accession: `docker run accession:latest` (without options, shows help message and exits, equivalent to [--help/-h])

```
docker run -v "`pwd`/local":/out accession:latest -f /out/ACCESSION.txt -d /out --logging-level INFO --log-file /out/accession.log
```

Here, we use volume driver options, where the first field is an existing directory on the host machine (a directory named *local* under the 
current directory, where we have placed required input files *e.g.* ACCESSION.txt), and the second field is the path where the file or directory 
are mounted in the container. Options with input files or paths that are passed to the program must be relative to the mounted directory in the container.

- wrangling: `docker run wrangling:latest` (without options, shows help message and exits, equivalent to [--help/-h])

```
docker run -v /mnt/smb/prj/DHART/DATA/GSEXXXXX1:/out wrangling:latest --dest /out --name docker -fmt TXT --txt-gene-cols gene_id --gene-var gene_symbol --do-normalise --normalization-method TMM --logging-level INFO --log-file /out/wrangling.log
```

Here, the host machine's directory is /mnt/smb/prj/DHART/DATA/GSEXXXXX1, which *must* be the path to the directory where GEO files reside, including output from accession and/or pre-processed bulk RNA-seq count matrix (the path specified by [--dest]).


### Singularity

We strongly recommend to use Singularity. However, at present, we do not have recipes to build Singularity images, but they can be built from 
Docker. First create a tarball of the Docker image, *e.g.* `docker save 5fad7a46ecfd -o docker/accession/accession.tar`
Build the singularity image `singularity build accession.sif docker-archive:///prj/DHART/dhart-utils/local/accession.tar`
Use the `--sandbox` option as `singularity build --sandbox accession docker-archive:///prj/DHART/dhart-utils/local/accession.tar`
for development, by accessing the container `singularity shell --bind /prj/DHART/dhart-utils/local/:/local accession`.
We use the `--bind` to mount volumes, as for Docker above.


## RNA bulk data download and processing

Data can only be read if it conforms to one of the provided data schemes. Using the `--normalization-method` flag, a normalization method can be specified. `TMM`, `DESeq` and `TPM` are all valid normalization methods, `TMM` being the default. When processing data using TPM normalization, gene lengths must be provided separately. See below for details. 

The file *observations.tab.gz* (output form `accession`) and the bulk data in a file named `GSEXXXXXX_read_count` (whith the correct GEO identifier, and any *TXT* extension) needs to be submitted for each dataset. 

**Note:** The input read count matrix *e.g.* `GSEXXXXXX_read_count.csv.gz` *MUST* be prepared by you. This is NOT automatically done by `accession` or `wrangling`. In some cases, data downloaded from GEO by `accession` will be ready to input, but you must first ensure that it satisfies the input format.


### Data Schemes

Data schmeme choices include the following:

- 'gene_id': Includes just the gene ID, this is later used to look up the gene symbol. Both fields are then written to the H5AD file.
- 'gene_symbol': Includes just the gene symbol.
- 'both': Both the gene symbol and gene ID are already present in the input scheme, both will be written to the H5AD file. 

**Note:** The column header for extra columns must be empty!

* gene_id

|    | <sample#1> | ... | <sample#n> |
|----|------------|-----|------------|
| ENS...   | 1.343      | ... | 1.472      |
| ENS...   | 2.872      | ... | 3.971      |
| ENS...   | 1.341      | ... | 1.349      |

* gene_symbol

|      | <sample#1> | ... | <sample#n> |
|------|------------|-----|------------|
|  A...    | 1.343      | ... | 1.472      |
|  B...    | 2.872      | ... | 3.971      |
|  C...    | 1.341      | ... | 1.349      |

* both

|    |      | <sample#1> | ... | <sample#n> |
|----|------|------------|-----|------------|
| ENS...   |  A...    | 1.343      | ... | 1.472      |
| ENS...   |  B...    | 2.872      | ... | 3.971      |
| ENS...   |  C...    | 1.341      | ... | 1.349      |

## Running the tests

Currently no test are implemented. We're working on adding integration/regression tests baed on some sample data.
Some minimal test data for bulk under /prj/DHART/DATA/test_data/RNAseq/GSEXXXXX1.

## Contributing


## Versioning

This project adheres to [Semantic Versioning](http://semver.org/). See [CHANGELOG](CHANGELOG.md).

## Notes

GEO classes and parsing functions in `geo_utils` are based on [GEOparse](https://github.com/guma44/GEOparse), which does
not seem to be maintained anymore. These were modified to deal with GEO SERIES and SAMPLES metadata only.
Interacting with GEO/SRA online (data download) is done with [pysradb](https://github.com/saketkc/pysradb).
Bulk RNA-seq normalization relies on standard functions implemented in the R packages **edgeR** and **DESeq2**.

