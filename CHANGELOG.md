# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **DockQ v2 interface quality assessment** (`dockq` feature)
  - Standard metric for CAPRI/CASP-multimer protein-protein interface evaluation
  - `dockq_to(native)` - Compute DockQ with automatic chain mapping
  - `dockq_to_with_options(native, options)` - Compute with custom options
  - `find_chain_mapping(model, native)` - Auto-detect chain correspondence via sequence alignment
  - `calculate_interface_dockq(model, native, chains, options)` - Score specific interfaces
  - Needleman-Wunsch sequence alignment for chain and residue matching
  - Interface contact detection (fnat, fnonnat, F1)
  - Interface RMSD (iRMSD) with Kabsch superposition of interface backbone atoms
  - Ligand RMSD (LRMSD) after receptor alignment
  - DockQ formula: (fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD/8.5)^2)) / 3
  - Quality classification: Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>=0.80)
  - Multi-interface support with contact-weighted averaging
  - Automatic chain permutation search (optimal for <=8 chains, greedy for larger)
  - `DockQResult`, `InterfaceResult`, `DockQOptions`, `DockQQuality` types
  - 21 integration tests + unit tests in each submodule
  - `dockq_demo.rs` example
- **Ligand pose quality validation** (`ligand-quality` feature)
  - PoseBusters-style geometry checks for protein-ligand complexes
  - `ligand_pose_quality(ligand_name)` - Validate a specific ligand's geometry
  - `all_ligand_pose_quality()` - Validate all ligands in structure
  - `get_ligand_names()` - List all ligand residue names
  - VDW radii-based clash detection with 0.75 × sum(vdW radii) threshold
  - Grid-based volume overlap calculation with 7.5% threshold
  - Cofactor clash detection with metal coordination support (more permissive thresholds)
  - CONECT record handling for covalent ligands (excludes bonded pairs from clash detection)
  - Water molecule exclusion (HOH, WAT, H2O, DOD)
  - `LigandPoseReport` with comprehensive validation results:
    - `min_protein_ligand_distance` - Closest approach in Angstroms
    - `clashes` - List of `AtomClash` with severity scoring
    - `protein_volume_overlap_pct` - Percentage of ligand volume overlapping protein
    - `passes_distance_check`, `passes_overlap_check`, `is_geometry_valid` - Pass/fail verdicts
  - `AtomClash` struct with protein/ligand atom details, distance, expected minimum, severity
  - Van der Waals radii table (Bondi/Alvarez radii for ~30 elements)
  - Covalent radii table (Cordero radii for metal coordination)
  - Helper functions: `vdw_radius()`, `covalent_radius()`, `is_metal()`
  - Full Python bindings: `PyLigandPoseReport`, `PyAtomClash` classes
  - New feature flag: `ligand-quality` (included in `analysis` feature group)

- **LDDT (Local Distance Difference Test) calculation** (`geometry` feature)
  - `lddt_to()` - Calculate LDDT score to a reference structure (superposition-free)
  - `lddt_to_with_options()` - LDDT with custom atom selection and options
  - `per_residue_lddt_to()` - Per-residue LDDT scores for quality analysis
  - `per_residue_lddt_to_with_options()` - Per-residue LDDT with custom options
  - `LddtOptions` struct with configurable `inclusion_radius` (default: 15.0 Å) and `thresholds` (default: [0.5, 1.0, 2.0, 4.0] Å)
  - `LddtResult` struct with global score (0.0-1.0), per-threshold scores, pair counts
  - `PerResidueLddt` struct for identifying poorly modeled regions
  - Same metric used by AlphaFold (pLDDT) and CASP evaluations
  - Translation and rotation invariant (superposition-free)
  - Full Python bindings: `LddtOptions`, `LddtResult`, `PerResidueLddt` classes

- **AlphaFold/pLDDT confidence score support** (`descriptors` feature)
  - `is_predicted()` - Detect AlphaFold/ESMFold predicted structures from metadata or B-factor patterns
  - `plddt_mean()` - Mean pLDDT confidence score across the structure
  - `plddt_scores()` - Raw pLDDT values (B-factors interpreted as confidence)
  - `per_residue_plddt()` - Per-residue pLDDT with confidence categories (VeryHigh >90, Confident 70-90, Low 50-70, VeryLow <50)
  - `low_confidence_regions(threshold)` - Identify disordered/uncertain regions
  - `high_confidence_regions(threshold)` - Identify well-predicted regions
  - `plddt_distribution()` - Fraction of residues in each confidence category
  - `ResiduePlddt` struct with `is_confident()`, `is_disordered()` helper methods
  - `ConfidenceCategory` enum with `is_reliable()`, `needs_caution()` methods
  - Full Python bindings: `PyConfidenceCategory`, `PyResiduePlddt` classes

- **Phi/Psi dihedral angles and Ramachandran analysis** (`descriptors` + `dssp` features)
  - `phi_psi_angles()` - Backbone dihedral angles (φ, ψ, ω) for all residues
  - `ramachandran_outliers()` - Identify residues with strained conformations
  - `ramachandran_statistics()` - Favored/allowed/outlier fractions for structure validation
  - `cis_peptide_bonds()` - Detect cis peptide bonds (ω ≈ 0°)
  - `ResidueDihedrals` struct with phi, psi, omega angles and Ramachandran classification
  - `RamachandranRegion` enum: Core, Allowed, Generous, Outlier, Glycine, Proline, PrePro
  - `RamachandranStats` with favored_fraction, allowed_fraction, outlier_fraction, cis counts
  - Proper IUPAC sign convention for dihedral angles
  - Full Python bindings: `PyRamachandranRegion`, `PyResidueDihedrals`, `PyRamachandranStats`

- **Hydrogen bond network API** (`descriptors` + `dssp` features)
  - `mainchain_hbonds()` - All backbone hydrogen bonds with energy and geometry
  - `hbonds_for_residue(chain, resid)` - H-bonds donated/accepted by specific residue
  - `hbond_statistics()` - Network analysis (total, intra-helical, beta-sheet, turns, long-range)
  - `MainchainHBond` struct with donor/acceptor info, energy (kcal/mol), N-O distance
  - `HBondType` classification: IntraHelical, BetaSheet, Turn, LongRange, InterChain
  - `ResidueHBonds` with donated and accepted H-bond lists
  - `HBondStats` with counts by type and mean energy
  - H-bond classification based on sequence separation and chain identity
  - Full Python bindings: `PyHBondType`, `PyMainchainHBond`, `PyResidueHBonds`, `PyHBondStats`

- **Protein-ligand interaction analysis** (`descriptors` feature)
  - `binding_site(ligand_name, distance)` - Contact residues within distance of ligand
  - `ligand_interactions(ligand_name)` - Full interaction profile with H-bonds, salt bridges, hydrophobic contacts
  - `all_ligand_interactions()` - Analyze all ligands in structure
  - `BindingSite` with contact residues sorted by distance, ligand info
  - `LigandInteractionProfile` with categorized interactions and `total_interactions()`, `has_interactions()` methods
  - `ContactResidue` with min_distance and contact count
  - `ProteinLigandHBond` with donor/acceptor atoms and distance
  - `SaltBridge` with charged residue and ligand atom info
  - `HydrophobicContact` with non-polar atom contacts
  - Detection thresholds: H-bonds ≤3.5Å, salt bridges ≤4.0Å, hydrophobic ≤4.0Å
  - Full Python bindings: `PyBindingSite`, `PyLigandInteractionProfile`, `PyContactResidue`, `PyProteinLigandHBond`, `PySaltBridge`, `PyHydrophobicContact`

## [0.6.0] - 2025-01-21

### Added
- **B-factor (temperature factor) analysis** (`descriptors` feature)
  - `b_factor_mean()` - Mean B-factor across all atoms
  - `b_factor_mean_ca()` - Mean B-factor for CA atoms only
  - `b_factor_min()`, `b_factor_max()` - Min/max B-factor values
  - `b_factor_std()` - Standard deviation of B-factors
  - `b_factor_profile()` - Per-residue B-factor statistics (ResidueBFactor struct)
  - `flexible_residues(threshold)` - Identify mobile/disordered regions (B > threshold)
  - `rigid_residues(threshold)` - Identify well-ordered regions (B < threshold)
  - `normalize_b_factors()` - Z-score normalization for cross-structure comparison
  - `b_factor_percentile(atom_serial)` - Get percentile rank of atom's B-factor
  - B-factor fields (mean, mean_ca, min, max, std) added to `StructureDescriptors`
  - Full Python bindings: `ResidueBFactor` class with all statistics
- **Async RCSB Downloads** (`rcsb-async` feature)
  - `download_multiple_async()` for concurrent bulk downloads with rate limiting
  - `download_structure_async()`, `download_pdb_string_async()`, `download_to_file_async()` for single async downloads
  - `AsyncDownloadOptions` for controlling concurrency (max_concurrent), rate limiting (rate_limit_ms), timeout, and retries
  - Preset options: `conservative()` (2 concurrent, 500ms delay) and `fast()` (20 concurrent, 25ms delay)
  - Builder pattern: `.with_max_concurrent()`, `.with_rate_limit_ms()`, `.with_timeout_secs()`, `.with_retries()`
  - Python bindings: `download_multiple()` function with `AsyncDownloadOptions` and `DownloadResult` types
  - Uses tokio for async runtime and futures for concurrent execution
  - Automatic retry with exponential backoff on transient failures
- **DSSP 4-like secondary structure assignment** (`dssp` feature)
  - Implements Kabsch & Sander algorithm with DSSP 4 updates (Hekkelman et al., 2025)
  - 9-state classification: H (α-helix), G (3₁₀-helix), I (π-helix), P (κ-helix/PPII), E (extended strand), B (β-bridge), T (turn), S (bend), C (coil)
  - `assign_secondary_structure()` returns full `SecondaryStructureAssignment` with per-residue details
  - `secondary_structure_string()` returns compact string representation (e.g., "HHHHEEEECCCC")
  - `secondary_structure_composition()` returns (helix_fraction, sheet_fraction, coil_fraction)
  - H-bond detection using Kabsch-Sander electrostatic energy model
  - PPII/κ-helix detection using backbone dihedral angles (φ = -75° ± 29°, ψ = +145° ± 29°)
  - Full Python bindings with `SecondaryStructure`, `ResidueSSAssignment`, `SecondaryStructureAssignment` classes
  - Iterator and indexing support for residue assignments in Python
- **PyMOL/VMD-style selection language** (`filter` feature)
  - `select()` method for flexible atom selection using familiar syntax
  - Supports: `chain A`, `name CA`, `resname ALA`, `resid 50`, `resid 1:100`, `element C`
  - Keywords: `backbone`, `protein`, `nucleic`, `water`, `hetero`, `hydrogen`, `all`
  - Numeric comparisons: `bfactor < 30.0`, `occupancy >= 0.5`
  - Boolean operators: `and`, `or`, `not`, parentheses for grouping
  - Examples: `chain A and name CA`, `(chain A or chain B) and bfactor < 30.0`
  - `validate_selection()` for syntax checking without execution
- **mmCIF writing support**
  - `write_mmcif_file()` for writing structures to mmCIF format files
  - `write_mmcif_string()` for writing to strings
  - `write_mmcif()` for writing to any `Write` implementor
  - `write_gzip_mmcif_file()` for compressed output (`gzip` feature)
  - Proper ATOM/HETATM record distinction based on residue type
  - Support for SEQRES (sequence) and SSBOND (disulfide bonds) categories
  - Full Python bindings: `write_mmcif_file()`, `write_mmcif_string()`, `write_gzip_mmcif_file()`

### Changed
- `analysis` feature now includes `dssp` feature
- **Python 3.13 support** - Pre-built wheels now available for Python 3.9-3.13 (closes #9)

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

## [0.1.0] - 2024-12-26

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
