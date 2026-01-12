# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2025-01-12

### Added
- **Distance matrix and contact map calculations** (`geometry` feature)
  - `distance_matrix()` for all-atom pairwise distance calculations
  - `distance_matrix_ca()` for CA-only distance matrix
  - `contact_map(threshold)` for all-atom contact detection
  - `contact_map_ca(threshold)` for CA contact maps
  - Returns efficient numpy arrays in Python bindings
- **RMSD calculation and structure superposition** (`geometry` feature)
  - `rmsd_to()` for RMSD calculation between structures
  - `align_to()` for structure superposition using Kabsch algorithm
  - `per_residue_rmsd_to()` for per-residue RMSD analysis
  - `AlignmentResult` with RMSD, rotation matrix, and translation vector
  - `PerResidueRmsd` for flexibility analysis
- **Python bindings for geometry features**
  - Full numpy integration for distance matrices and contact maps
  - Structure alignment returns (aligned_structure, AlignmentResult) tuple

### Changed
- **Platform-specific Python builds**
  - Linux wheels built with `core` features (without RCSB) to avoid OpenSSL cross-compilation issues
  - macOS and Windows wheels built with `full` features (including RCSB)
  - All other features (parsing, filtering, descriptors, geometry, numpy) work on all platforms
- Updated dependencies

### Fixed
- Linux wheel builds now work reliably without OpenSSL dependency issues

### Notes
- v0.4.0 and v0.4.1 were skipped (accidentally published during CI debugging)

## [0.3.0] - 2025-01-04

### Added
- **Gzip decompression support** (`gzip` feature)
  - `parse_gzip_pdb_file()` for compressed PDB files (`.pdb.gz`, `.ent.gz`)
  - `parse_gzip_mmcif_file()` for compressed mmCIF files (`.cif.gz`)
  - `parse_gzip_structure_file()` with automatic format detection
  - Reader-based variants for streaming: `parse_gzip_pdb_reader()`, `parse_gzip_mmcif_reader()`
  - Handles PDB archive naming conventions (e.g., `pdb1ubq.ent.gz`)
- **Full PDB archive benchmark example** (`examples/full_pdb_benchmark.rs`)
  - Parallel processing of large datasets with Rayon
  - Detailed error categorization and failure logging
  - Comprehensive timing statistics (mean, median, P99, throughput)
  - Progress reporting with ETA for long-running operations

### Changed
- Updated `full` feature to include `gzip` support

### Validated
- **Full PDB Archive Benchmark**: Successfully parsed 230,655 structures with 100% success rate
  - 2,057,302,767 total atoms processed
  - 92 structures/second throughput (parallel processing)
  - Zero parsing failures across entire RCSB PDB mirror

## [0.2.0] - 2025-01-02

### Added
- **mmCIF format support** with automatic format detection
  - `parse_mmcif_file()` for explicit mmCIF parsing
  - `parse_structure_file()` with auto-detection based on content
  - Full support for `_atom_site`, `_entity_poly_seq`, `_struct_disulfid` categories
- **Filter module** (`filter` feature)
  - Fluent API for structure filtering: `remove_ligands()`, `keep_only_chain()`, `keep_only_ca()`, `keep_only_backbone()`
  - Cleaning operations: `normalize_chain_ids()`, `reindex_residues()`, `center_structure()`
  - Coordinate extraction: `get_ca_coords()`
- **Descriptors module** (`descriptors` feature)
  - Geometric descriptors: `radius_of_gyration()`, `max_ca_distance()`, `ca_density()`, `compactness_index()`
  - Composition analysis: `aa_composition()`, `hydrophobic_ratio()`, `glycine_ratio()`, `polar_ratio()`
  - Secondary structure estimation: `secondary_structure_ratio()`
  - Combined descriptor extraction: `structure_descriptors()`
- **Quality module** (`quality` feature)
  - Structure quality assessment: `quality_report()`
  - Quality checks: `has_altlocs()`, `has_multiple_models()`, `has_ca_only()`, `has_hydrogens()`
  - Resolution parsing from REMARK 2
  - Analysis readiness evaluation: `is_analysis_ready()`
- **Summary module** (`summary` feature)
  - Unified structure summaries combining quality and descriptors
  - Batch processing: `batch_summarize()`
  - CSV export: `summaries_to_csv()`
- **RCSB PDB integration** (`rcsb` feature)
  - Search API with query builder: `SearchQuery`, `rcsb_search()`
  - Structure download: `download_structure()`, `download_pdb_string()`
  - Multiple download: `download_multiple()`
  - Support for both PDB and mmCIF formats
- **Comprehensive documentation**
  - `guide` module with detailed usage examples
  - `docs/GETTING_STARTED.md` quick start guide
  - Example programs for common workflows

### Changed
- Reorganized module structure for feature-gated components
- Enhanced error types with more detailed context

## [0.1.0] - 2025-03-26

### Added
- Initial release
- Complete PDB file parsing support
  - ATOM/HETATM records
  - MODEL/ENDMDL records
  - SEQRES records
  - CONECT records
  - SSBOND records
  - REMARK records
  - HEADER and TITLE metadata
- Utility functions for structure analysis
  - Chain identification
  - Residue sequence extraction
  - Connectivity analysis
  - Model-based structure organization
  - Disulfide bond analysis
- Comprehensive error handling
- Property-based testing
- Performance benchmarks
- Optional parallel processing support
- Optional geometric analysis support

## [Unreleased]

### Added
- None

### Changed
- None

### Fixed
- None
