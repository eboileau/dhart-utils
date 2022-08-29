# Change Log

Changes to dhart-utils will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and 
this project adheres to [Semantic Versioning](http://semver.org/).

## [0.1.0] - 2022-08-25

### Changed

- Package structure, setup install, first version (alpha)

### Fixed

- Containers

## [Unreleased] - 2022-05

### Changed

- `wrangling.py` and `data_utils.py` (bulk processing)
- Annotation source, input file format.

### Added

- Assembly, annotation release.

### Fixed

- Docker
- TPM normalization
- Column names in adata.var to match what is expected in gEAR. TODO: if n_vars == 1, ids or symbols, this actualy does make a difference... they should be symbols!


---

### TODOs
- This would probably not happen if we select our datasets, but if there are no supp data, SOFT file will be downloaded, and
an error will be thrown in `get_supp_files`. We could add a check and terminate gracefully.
