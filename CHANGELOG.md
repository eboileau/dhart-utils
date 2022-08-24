# Change Log

Changes to dhart-utils will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/).
This project currently does NOT adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] - 2022-05

### Changed
- `wrangling.py` and `data_utils.py` (bulk processing)
- Annotation source, input file format.

### Added
- Assembly, annotation release.

### Fixed
- TPM normalization
- Column names in adata.var to match what is expected in gEAR. TODO: if n_vars == 1, ids or symbols, this actualy does make a difference... they should be symbols!
