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

## Medium Impact

### DSSP-like Secondary Structure Assignment
- Replace current heuristic with proper H-bond based assignment
- Compute secondary structure from coordinates alone
- Return per-residue SS codes (H, E, C, etc.)

### Async RCSB Downloads
- Add async variants of download functions
- Better for bulk downloading hundreds of structures
- `download_multiple_async()` with concurrency control

### B-factor Analysis
- `structure.b_factor_profile()` → per-residue B-factors
- Identify flexible regions
- Normalize B-factors across structures

## Lower Priority

### Binding Site Detection
- Identify residues near ligands/HETATM
- `structure.binding_site(ligand_name, distance_cutoff)`
- Useful for drug discovery applications

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
