# PDBRust Roadmap

Future development ideas for PDBRust. Priority will be determined based on user feedback.

## Completed

### Python Bindings (PyO3) ✅
- Created `pdbrust-python` package, published to PyPI (`pip install pdbrust`)
- Full API: parsing, filtering, descriptors, quality, RCSB search/download
- Numpy integration: `get_coords_array()`, `get_ca_coords_array()`
- Multi-platform wheels: Linux, macOS, Windows (Python 3.9-3.12)
- GitHub Actions CI/CD for automated releases

### Contact Maps / Distance Matrices ✅
- `structure.distance_matrix_ca()` → 2D f64 matrix (CA atoms)
- `structure.contact_map_ca(threshold)` → 2D boolean matrix (default: 8.0 Å)
- `structure.distance_matrix()` → 2D f64 matrix (all atoms)
- `structure.contact_map(threshold)` → 2D boolean matrix (default: 4.5 Å)
- Python bindings return numpy arrays (N×N)
- Essential for ML applications (GNNs, protein transformers)

### RMSD / Structure Superposition ✅
- Kabsch algorithm for optimal alignment using nalgebra SVD
- `structure.rmsd_to(other)` → f64 (CA atoms)
- `structure.align_to(target)` → (aligned_structure, AlignmentResult)
- `structure.per_residue_rmsd_to(target)` → Vec<PerResidueRmsd>
- Flexible atom selection: CA only (default), backbone, all atoms, custom
- Python bindings with AtomSelection, AlignmentResult, PerResidueRmsd types
- Under `geometry` feature flag (requires nalgebra)

### mmCIF Writing ✅
- `write_mmcif_file()` function for file output
- `write_mmcif_string()` for string output
- `write_mmcif()` for generic writer output
- `write_gzip_mmcif_file()` for compressed output (gzip feature)
- Full Python bindings included
- Supports ATOM/HETATM, SEQRES, SSBOND data

### Selection Language (Query DSL) ✅
- PyMOL/VMD-style selection syntax: `structure.select("chain A and name CA")`
- Basic selectors: `chain`, `name`, `resname`, `resid`, `element`
- Range selections: `resid 1:100`
- Keywords: `backbone`, `protein`, `nucleic`, `water`, `hetero`, `hydrogen`
- Boolean operators: `and`, `or`, `not`, parentheses for grouping
- Numeric comparisons: `bfactor < 30.0`, `occupancy >= 0.5`
- Zero external dependencies (hand-written recursive descent parser)
- Full Python bindings included

### DSSP 4-like Secondary Structure Assignment ✅
- Implements Kabsch & Sander algorithm with DSSP 4 updates (Hekkelman et al., 2025)
- H-bond detection using electrostatic energy model (E < -0.5 kcal/mol threshold)
- 9-state classification: H (α-helix), G (3₁₀-helix), I (π-helix), P (κ-helix/PPII), E (extended), B (β-bridge), T (turn), S (bend), C (coil)
- PPII/κ-helix detection using backbone dihedral angles (φ = -75° ± 29°, ψ = +145° ± 29°)
- `structure.assign_secondary_structure()` → `SecondaryStructureAssignment`
- `structure.secondary_structure_string()` → String (e.g., "HHHHEEEECCCC")
- `structure.secondary_structure_composition()` → (helix_fraction, sheet_fraction, coil_fraction)
- Full Python bindings with iterator and indexing support
- Under `dssp` feature flag (included in `analysis`)

### Async RCSB Downloads ✅
- Async variants of download functions for efficient bulk downloading
- `download_multiple_async()` with concurrency control via `AsyncDownloadOptions`
- Configurable: max_concurrent (default: 5), rate_limit_ms (default: 100ms), timeout, retries
- Preset options: `conservative()` for rate-limited scenarios, `fast()` for high-throughput
- Automatic retry with exponential backoff on transient failures
- Python bindings: `download_multiple()` with `AsyncDownloadOptions` and `DownloadResult`
- Under `rcsb-async` feature flag (included in `full`)

### B-factor Analysis ✅
- `structure.b_factor_mean()` → mean B-factor across all atoms
- `structure.b_factor_mean_ca()` → mean B-factor for CA atoms only
- `structure.b_factor_min()`, `b_factor_max()`, `b_factor_std()` → basic statistics
- `structure.b_factor_profile()` → per-residue B-factor statistics (mean, min, max)
- `structure.flexible_residues(threshold)` → identify mobile/disordered regions
- `structure.rigid_residues(threshold)` → identify well-ordered regions
- `structure.normalize_b_factors()` → Z-score normalization for cross-structure comparison
- `structure.b_factor_percentile(atom_serial)` → get percentile rank of atom's B-factor
- `ResidueBFactor` struct with chain_id, residue_seq, residue_name, and B-factor stats
- Full Python bindings with `ResidueBFactor` class
- B-factor fields added to `StructureDescriptors`
- Under existing `descriptors` feature flag

### AlphaFold/pLDDT Support ✅
- Detect AlphaFold/ESMFold predicted structures from B-factor range heuristics
- Interpret B-factor column as pLDDT confidence scores (0-100)
- `structure.is_predicted()` → detect AI-predicted structures
- `structure.plddt_mean()` → mean confidence score
- `structure.per_residue_plddt()` → per-residue pLDDT with confidence categories
- `structure.low_confidence_regions(threshold)` → identify disordered regions (pLDDT < threshold)
- `structure.high_confidence_regions(threshold)` → identify well-predicted regions
- `structure.plddt_distribution()` → fraction in each confidence category (VeryHigh, Confident, Low, VeryLow)
- `ConfidenceCategory` enum with `is_reliable()` and `needs_caution()` methods
- `ResiduePlddt` struct with plddt, plddt_min, plddt_max, confidence_category
- Full Python bindings with `ConfidenceCategory`, `ResiduePlddt` classes
- Under existing `descriptors` feature flag

### Phi/Psi Dihedral Angles & Ramachandran Analysis ✅
- Expose DSSP's internal dihedral calculations to users
- `structure.phi_psi_angles()` → Vec<ResidueDihedrals> for all backbone dihedrals
- `structure.ramachandran_outliers()` → residues in unfavored regions
- `structure.ramachandran_statistics()` → RamachandranStats with favored/allowed/outlier counts and fractions
- Cis peptide bond detection via `ResidueDihedrals.is_cis_peptide()` and `is_trans_peptide()`
- `RamachandranRegion` enum: Core, Allowed, Generous, Outlier, Glycine, Proline, PrePro, Unknown
- Proper IUPAC sign convention for phi/psi angles
- Full Python bindings with `ResidueDihedrals`, `RamachandranRegion`, `RamachandranStats` classes
- Requires both `descriptors` and `dssp` feature flags

### Hydrogen Bond Network API ✅
- Expose DSSP's H-bond detection with user-friendly API
- `structure.mainchain_hbonds()` → Vec<MainchainHBond> for all backbone H-bonds
- `structure.hbonds_for_residue(chain, resid)` → ResidueHBonds with donated/accepted lists
- `structure.hbond_statistics()` → HBondStats with counts by type and mean energy
- `HBondType` enum: IntraHelical, BetaSheet, Turn, LongRange, InterChain
- `MainchainHBond` struct with donor/acceptor info, energy, distance, sequence separation
- Methods: `is_strong()` (E < -1.5), `is_helical()`, `is_beta_sheet()`
- Full Python bindings with `MainchainHBond`, `ResidueHBonds`, `HBondStats`, `HBondType` classes
- Requires both `descriptors` and `dssp` feature flags

### Protein-Ligand Interaction Analysis ✅
- `structure.binding_site(ligand_name, distance)` → BindingSite with contact residues
- `structure.ligand_interactions(ligand_name)` → LigandInteractionProfile
- `structure.all_ligand_interactions()` → analyze all ligands in structure
- Detection of H-bonds (≤3.5 Å), salt bridges (≤4.0 Å), hydrophobic contacts (≤4.0 Å)
- `ContactResidue` with min_distance and num_contacts
- `ProteinLigandHBond` with donor/acceptor identification
- `SaltBridge` with charge polarity information
- `HydrophobicContact` for non-polar interactions
- Full Python bindings with all interaction types
- Under existing `descriptors` feature flag

### LDDT (Local Distance Difference Test) ✅
- `structure.lddt_to(reference)` → f64 (0.0 to 1.0, higher is better)
- `structure.lddt_to_with_options(reference, selection, options)` → LddtResult with detailed statistics
- `structure.per_residue_lddt_to(reference)` → Vec<PerResidueLddt> for quality analysis
- Superposition-free metric (invariant to rotation/translation)
- Configurable inclusion radius (default: 15.0 Å) and thresholds (default: 0.5, 1.0, 2.0, 4.0 Å)
- Same metric used by AlphaFold (pLDDT) and CASP structure prediction evaluations
- `LddtResult` with global score, per-threshold scores, and pair counts
- `PerResidueLddt` for identifying poorly modeled regions
- Full Python bindings: `LddtOptions`, `LddtResult`, `PerResidueLddt` classes
- Under `geometry` feature flag (requires nalgebra)

## Future Ideas

### Surface Area Calculation
- Solvent accessible surface area (SASA)
- Per-residue and per-atom breakdown
- Requires rolling ball algorithm

### Symmetry Operations
- Parse and apply BIOMT records
- Generate biological assemblies
- `structure.biological_assembly()`

### Trajectory Support
- Parse multi-frame trajectories (e.g., from MD simulations)
- Memory-efficient streaming for large trajectories
- Basic trajectory analysis (RMSD over time, etc.)

## Community Requested

*This section will be populated based on GitHub issues and user feedback.*

---

## Contributing

If you'd like to work on any of these features, please:
1. Open an issue to discuss the approach
2. Reference this roadmap in your PR

Feedback and feature requests welcome at: https://github.com/HFooladi/pdbrust/issues
