# PDBRust Roadmap

Development plan for PDBRust. Last updated: 2026-09-16.

## Current Focus: Hardening & Validation (v0.7.1 → v0.8.0)

Through v0.7.0 PDBRust gained a broad feature set. The next cycle makes those features **correct, robust, and
validated against reference tools** before adding new ones. A code review in September 2026 found issues that can
make results silently wrong on real-world files. Fixing them comes first.

### Known Issues Being Fixed

| Area | Issue | Target |
|------|-------|--------|
| PDB parser | SSBOND symmetry operators misread (`1555` → `15`); panics on some short HEADER/TITLE/REMARK lines; non-numeric REMARKs (e.g. GROMACS output) abort parsing | ✅ Fixed in v0.7.1 |
| Writers | mmCIF writer does not quote values (a blank chain ID corrupts columns); PDB writer column-alignment issues | ✅ Fixed in v0.7.1 |
| mmCIF parser | Single-value items after the first `loop_` are dropped (title, resolution, cell missing); no multi-line `;` text fields | ✅ Fixed in v0.7.1 |
| mmCIF parser | Single-quoted values in loop rows are split at spaces; `pdbx_PDB_model_num` ignored | v0.7.2 |
| mmCIF IDs | Chain/residue IDs come from `label_*` fields, so they differ from the same entry in PDB format | v0.8.0 (author IDs by default, label IDs kept) |
| Multi-model files | Atoms are stored twice, and analyses run on all NMR models combined | v0.7.2 / v0.8.0 |
| Selections | `AtomSelection::CaOnly` also matches calcium ions; `Backbone` matches water oxygens | v0.7.2 |
| DSSP | Diverges from mkdssp: virtual-H placement, β-bridge patterns, helix-start rule, PPII dihedral sign | v0.8.0 (sign fix in v0.7.2) |
| DockQ | LRMSD matches atoms by index; interfaces that fail are skipped instead of scored | v0.8.0 |
| lDDT | Index-based atom matching; same-residue pairs included; no symmetric-atom handling | v0.8.0 |
| Python wheels | Linux wheels are built without the RCSB feature ([#8](https://github.com/HFooladi/pdbrust/issues/8)) | ✅ Fixed in v0.7.1 (switch to rustls) |

### Milestones

**Maintenance (no release) ✅:** modernized CI (current GitHub Actions, feature-combination checks, a Python test
job), committed `Cargo.lock`, switched `reqwest` to rustls (fixes #8), and merged PR #16. Removing unused
dependencies moves to v0.7.2.

**v0.7.1 — Safety patch (non-breaking) ✅**
- Panic-free PDB parsing (safe fixed-column access), with tolerant handling of blank or odd fields; hybrid-36
  numbers read and written
- Writer fixes: mmCIF value quoting and PDB column layout; files are checked by reading them back with gemmi
- mmCIF title, resolution and other single-value items read correctly
- No-panic property tests and PDB/mmCIF round-trip tests; first Python test suite; Python 3.14 wheels (3.9 dropped)
- Ships the unreleased molecular inventory

**v0.7.2 — CIF tokenizer + validation harness (non-breaking)**
- Multi-member gzip, `.gz` auto-detection, `from_file` format auto-detection (moved from v0.7.1)
- Validate PDB IDs used as download file names; remove unused dependencies
- A spec-compliant CIF tokenizer (text fields, quoting rules, multiple data blocks) behind the existing API
- Multi-model data kept per model, plus forward-compatible `atoms()` / `models()` accessors
- A `validation/` harness with a parser differential test against gemmi
- Fuzzing (cargo-fuzz), plus small science fixes (dihedral sign, hydrogen detection, Ramachandran wrap)

**v0.8.0 — Correctness release (the only planned breaking release before 1.0)**
- Data model: `models` as the single atom store (`atoms()` = first model), author IDs by default with label IDs
  kept, new atom fields (formal charge, segment ID), and structure metadata (resolution, methods, unit cell,
  entities, links)
- `ParseOptions` (lenient by default with warnings; strict mode) and error types with line numbers
- Residue/chain views, a shared `ResidueId`, and an internal spatial index
- DSSP parity with mkdssp, lDDT parity with OpenStructure, DockQ parity with DockQ v2
- API consistency pass with deprecation shims and a `MIGRATION.md`; Python bindings released the same day, with
  type stubs

**Validation targets for v0.8.0** (reference outputs committed as golden tests):

| Component | Reference | Target |
|-----------|-----------|--------|
| Parser (PDB + mmCIF) | gemmi | 100% agreement on a stratified sample (documented differences only) |
| DSSP | mkdssp 4.4 | ≥ 99.5% 8-state, ≥ 99.8% 3-state per-residue agreement |
| lDDT | OpenStructure | abs(Δ global score) ≤ 0.005 |
| DockQ | DockQ v2 | abs(ΔDockQ) ≤ 0.01 for ≥ 99% of interfaces |
| Ligand pose checks | PoseBusters (implemented checks) | ≥ 99% pass/fail agreement |

Performance claims will be re-measured against gemmi, Biotite, Biopython, and pdbtbx using a documented protocol.

**v0.9.x — Structure completeness (additive):** see [Future Work](#future-work-after-v080).

**v1.0 — API freeze:** semver checks in CI, a documentation site, an MSRV policy, and removal of the v0.8
deprecation shims.

### API Stability Policy

v0.8.0 is the single planned breaking release. v0.7.2 already adds the new accessors, so code can migrate early.
After v0.8.0, changes are additive only until 1.0.

## Completed

> Some completed features have known issues (DSSP, DockQ, lDDT, mmCIF parsing). They are listed under
> [Known Issues Being Fixed](#known-issues-being-fixed) and will be fixed in v0.7.x–v0.8.0.

### Python Bindings (PyO3) ✅
- Created `pdbrust-python` package, published to PyPI (`pip install pdbrust`)
- Full API: parsing, filtering, descriptors, quality, RCSB search/download
- Numpy integration: `get_coords_array()`, `get_ca_coords_array()`
- Multi-platform wheels: Linux, macOS, Windows (Python 3.9-3.13)
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

### Ligand Pose Quality (PoseBusters-style Geometry Checks) ✅
- Validate protein-ligand complex geometry using PoseBusters-inspired criteria
- VDW radii-based clash detection (0.75 × sum of vdW radii threshold)
- Grid-based volume overlap calculation (7.5% threshold, 0.8 vdW scaling)
- Cofactor clash detection with metal coordination support
- CONECT record handling for covalent ligands
- `structure.ligand_pose_quality(ligand_name)` → Option<LigandPoseReport>
- `structure.all_ligand_pose_quality()` → Vec<LigandPoseReport>
- `structure.get_ligand_names()` → Vec<String>
- `LigandPoseReport` with clashes, overlap %, and pass/fail verdicts
- `AtomClash` with severity scoring and detailed atom information
- Van der Waals radii (Bondi/Alvarez) and covalent radii (Cordero) tables
- Full Python bindings: `PyLigandPoseReport`, `PyAtomClash` classes
- Under `ligand-quality` feature flag (included in `analysis`)

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

### Molecular Inventory ✅ (v0.7.1)
- One-call breakdown of structure contents — chains, ligands, water, ions
- `structure.molecular_inventory()` → `MolecularInventory`
- Per-chain summary: `ChainInventory` with type (Protein, NucleicAcid, Mixed, Water, Other), residue/atom counts
- Ligand detection: `LigandInfo` with name, chain, residue sequence, atom count (ions excluded)
- Global counts: protein/nucleic/water/het atoms, chain counts, water molecules
- Pretty-print via `Display` trait
- No feature flags required — works on any `PdbStructure`
- Full Python bindings: `MolecularInventory`, `ChainInventory`, `ChainType`, `LigandInfo`

### DockQ v2 Interface Quality Assessment ✅
- Standard metric for CAPRI/CASP-multimer protein-protein interface evaluation
- `structure.dockq_to(native)` → DockQResult with per-interface scores and overall DockQ
- `structure.dockq_to_with_options(native, options)` → DockQResult with custom thresholds/mapping
- `find_chain_mapping(model, native)` → auto-detect chain correspondence via sequence alignment
- `calculate_interface_dockq(model, native, chains, options)` → score specific interfaces
- Needleman-Wunsch sequence alignment for chain and residue matching
- Interface contact detection: fnat (fraction native contacts), fnonnat, F1 score
- iRMSD (interface RMSD) with Kabsch superposition of interface backbone atoms
- LRMSD (ligand RMSD) after receptor alignment
- DockQ formula: (fnat + 1/(1+(iRMSD/1.5)²) + 1/(1+(LRMSD/8.5)²)) / 3
- Quality classification: Incorrect (<0.23), Acceptable (0.23–0.49), Medium (0.49–0.80), High (≥0.80)
- Multi-interface support with contact-weighted averaging
- Automatic chain permutation search (optimal for ≤8 chains, greedy for larger)
- `DockQResult`, `InterfaceResult`, `DockQOptions`, `DockQQuality`, `ChainMappingStrategy` types
- Under `dockq` feature flag (requires `geometry`/nalgebra, included in `analysis`)

## Future Work (after v0.8.0)

### v0.9.x — Structure Completeness (additive; the fields these need are reserved in v0.8.0)

Listed in dependency order.

#### Spatial Index and Neighbor Queries
- Cell-list neighbor search shared by contacts, clashes, interactions, lDDT, DSSP, and DockQ (replaces O(n²) loops)
- `structure.neighbors_within(point, radius)`, `structure.contacts_between(sel1, sel2, cutoff)`
- Selection language: `within X of ...`, `around`, `byres`

#### Unit Cell, Symmetry, and Biological Assemblies
- Unit cell and space group from CRYST1 (PDB) and `_cell`/`_symmetry` (mmCIF)
- Parse and apply REMARK 350 BIOMT records (PDB) and `_pdbx_struct_assembly_gen` + `_pdbx_struct_oper_list`
  (mmCIF), including operator expressions such as `(1-60)(61-88)`
- `structure.biological_assembly(id)` → full biological unit with deterministic chain naming
- `structure.symmetry_mates()` → crystallographic neighbors
- Critical for homo-oligomers (most proteins function as multimers); validated against gemmi

#### Solvent Accessible Surface Area (SASA)
- Shrake-Rupley algorithm with per-atom and per-residue breakdown
- Relative SASA (Tien et al. 2013 maximum values) and buried/exposed classification
- `structure.sasa()`, `structure.per_residue_sasa()`, `structure.buried_residues(threshold)`
- Validated against FreeSASA

#### Missing Residues, Missing Atoms, and Chain Breaks
- REMARK 465/470 (PDB) and `_pdbx_unobs_or_zero_occ_residues` (mmCIF)
- SEQRES/entity sequence ↔ modeled residue alignment
- Chain-break detection from C–N distances

#### Clashscore / Steric Clashes
- Heavy-atom van der Waals overlap (≥ 0.4 Å, MolProbity-like), documented as an approximation without
  hydrogens/Probe
- `structure.clashscore()` → clashes per 1000 atoms; `structure.steric_clashes()` → atom pairs with overlap
- Shares clash code with `ligand-quality`

#### Author Annotations and Links
- HELIX/SHEET and `_struct_conf`/`_struct_sheet_range` (author-assigned secondary structure)
- LINK and `_struct_conn` covalent/metal links (also gives ligands bond topology)
- MODRES, ANISOU

#### BinaryCIF
- `.bcif` reader behind a `bcif` feature, reusing the CIF data model; BinaryCIF downloads from RCSB/PDBe/AFDB

### Later / Parking Lot

Deferred during the hardening cycle; to be re-prioritized from user feedback.

- **Sequence alignment API**: expose the Needleman-Wunsch used by DockQ as `structure.align_sequence_to(other)` and
  `structure.rmsd_aligned(other)`
- **Ligand chemistry**: Chemical Component Dictionary (bond orders/topology) → the remaining PoseBusters checks
  (bond lengths, angles, planarity, chirality)
- **Electrostatics / partial charges**: residue-template charges (AMBER, CHARMM), potential on a grid
- **mmCIF dictionary validation** against PDBx/mmCIF
- **Trajectory support** (DCD, XTC, TRR), with streaming and per-frame RMSD/lDDT
- **Symmetry-expanded RMSD/lDDT**
- **Domain detection**
- **ML/data pipelines**: featurization, Arrow/Parquet export, interop with Biotite/gemmi
- **Additional sources**: AlphaFold DB and PDB-REDO clients

## Community Requested

- [#8](https://github.com/HFooladi/pdbrust/issues/8) — Enable the RCSB feature in Linux Python wheels (scheduled:
  switch `reqwest` to rustls)
- [PR #16](https://github.com/HFooladi/pdbrust/pull/16) — Fix a panic on REMARK lines with no content (to be merged
  during the maintenance phase)

---

## Contributing

If you'd like to work on any of these features, please:
1. Open an issue to discuss the approach
2. Reference this roadmap in your PR

Feedback and feature requests welcome at: https://github.com/HFooladi/pdbrust/issues
