# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PDBRust is a Rust library for parsing and analyzing Protein Data Bank (PDB) and mmCIF files. It provides comprehensive support for molecular structure data with robust error handling, structural analysis, quality assessment, and RCSB PDB integration.

### Core Architecture

- **Dual Format Support**: Supports both PDB (legacy text format) and mmCIF (modern dictionary format) with automatic format detection
- **Unified Data Model**: Both formats convert to a common `PdbStructure` representation
- **Feature-Gated Modules**: Optional functionality behind feature flags for minimal compile times
- **High Performance**: 40-260x faster than equivalent Python implementations

### Module Structure

```
src/
├── core/           # Core data structures
│   ├── pdb.rs          # PdbStructure definition
│   ├── mmcif.rs        # mmCIF parsing
│   └── mmcif_converter.rs
├── parser/         # Format-specific parsers with auto-detection
│   ├── pdb/            # PDB format parser
│   ├── mmcif/          # mmCIF format parser
│   └── gzip.rs         # [feature: gzip] Gzip-compressed file support
├── records/        # Record types (Atom, Residue, Chain, Model, SSBond, etc.)
├── writer/         # PDB file output
├── error.rs        # Error handling with thiserror
├── guide.rs        # Comprehensive user guide (accessible via cargo doc)
├── filter/         # [feature: filter] Filtering and cleaning operations
├── descriptors/    # [feature: descriptors] Structural descriptors
├── quality/        # [feature: quality] Quality assessment
├── summary/        # [feature: summary] Unified structure summaries
├── rcsb/           # [feature: rcsb] RCSB PDB search and download
└── geometry/       # [feature: geometry] RMSD and structure superposition

pdbrust-python/     # Python bindings (PyO3)
├── Cargo.toml          # Rust dependencies for Python bindings
├── pyproject.toml      # Python package configuration (maturin)
├── src/
│   ├── lib.rs          # PyO3 module entry point
│   ├── structure.rs    # PyPdbStructure wrapper
│   ├── atom.rs         # PyAtom wrapper
│   ├── records.rs      # PySSBond, PySeqRes, etc.
│   ├── parsing.rs      # Parsing functions
│   ├── error.rs        # Python exception mapping
│   ├── descriptors.rs  # StructureDescriptors bindings
│   ├── quality.rs      # QualityReport bindings
│   ├── summary.rs      # StructureSummary bindings
│   ├── rcsb.rs         # RCSB search/download bindings
│   ├── geometry.rs     # RMSD/alignment bindings
│   └── numpy_support.rs # Numpy array integration
└── python/pdbrust/
    ├── __init__.py     # Python exports
    └── py.typed        # PEP 561 type hint marker

docs/
└── GETTING_STARTED.md  # Quick start guide for new users

examples/
├── pdb_files/          # Sample PDB files for testing
├── analysis_workflow.rs
├── filtering_demo.rs
├── rcsb_workflow.rs
├── batch_processing.rs
├── full_pdb_benchmark.rs  # Full PDB archive benchmark
└── ...                 # Additional examples

.github/workflows/
├── rust.yml            # Rust CI/CD
└── python-publish.yml  # Python wheel building and PyPI publishing

benchmark_results/      # Results from full PDB archive benchmark
├── benchmark_report.txt
├── failures.tsv
└── timing_histogram.txt
```

### Feature Flags

| Feature | Description | Dependencies |
|---------|-------------|--------------|
| `filter` | Filtering, extraction, and cleaning operations | - |
| `descriptors` | Structural descriptors (Rg, composition, geometry) | - |
| `quality` | Quality assessment and reports | - |
| `summary` | Unified summaries combining quality + descriptors | `descriptors`, `quality` |
| `rcsb` | RCSB PDB search API and file download | `reqwest`, `serde`, `serde_json` |
| `gzip` | Parse gzip-compressed files (.ent.gz, .pdb.gz) | `flate2` |
| `parallel` | Parallel processing with Rayon | `rayon` |
| `geometry` | Geometric analysis with nalgebra | `nalgebra` |
| `analysis` | All analysis features | `filter`, `descriptors`, `quality`, `summary` |
| `full` | Everything | `parallel`, `geometry`, `analysis`, `rcsb`, `gzip` |

## Development Commands

### Building and Testing
```bash
# Build the project
cargo build

# Run all tests
cargo test

# Run tests with all features
cargo test --all-features

# Run tests for a specific feature
cargo test --features filter
cargo test --features rcsb

# Run network tests (RCSB download/search)
cargo test --features rcsb -- --ignored

# Run benchmarks
cargo bench
```

### Code Quality
```bash
# Format code
cargo fmt

# Run clippy lints (CI uses -D warnings flag)
cargo clippy --all-targets --all-features -- -D warnings

# Check documentation
cargo doc --no-deps --all-features

# Security audit
cargo audit
```

### GitHub Actions / CI
**Important:** Use `gh_cli` instead of `gh` for GitHub CLI commands in this environment.

```bash
# List recent workflow runs
gh_cli run list --limit 5

# Watch a specific run (with exit status)
gh_cli run watch <run_id> --exit-status

# View logs for failed jobs
gh_cli run view <run_id> --log-failed

# Rerun failed jobs
gh_cli run rerun <run_id> --failed
```

### Benchmarking
```bash
# Run Rust benchmark
cargo run --release --features "filter,descriptors" --example rust_benchmark

# Run Python comparison benchmark
python3 benchmarks/python_benchmark.py
```

## Testing Strategy

- **Unit tests**: In individual modules (`#[cfg(test)]` blocks)
- **Integration tests**: In `tests/` directory, one file per feature module
- **Property-based tests**: Using `proptest` for robust validation
- **Benchmark tests**: Using `criterion` for performance tracking
- **Network tests**: Marked with `#[ignore]` for RCSB API tests

Test files:
- `tests/filter_tests.rs` - 18 tests
- `tests/selection_tests.rs` - 30 tests (selection language)
- `tests/descriptors_tests.rs` - 17 tests
- `tests/quality_tests.rs` - 18 tests
- `tests/summary_tests.rs` - 18 tests
- `tests/rcsb_tests.rs` - 28 tests (11 network tests ignored by default)

## Key Data Structures

### Core Types
- `PdbStructure`: Main structure containing all parsed data
- `Atom`: Individual atom records with coordinates and metadata
- `Residue`: Amino acid/nucleotide residue information
- `Chain`: Protein/nucleic acid chain organization
- `Model`: Multi-model structures (NMR ensembles)
- `SSBond`: Disulfide bond connectivity

### Feature-Specific Types
- `SelectionError` (filter): Selection language parsing errors
- `QualityReport` (quality): Structure quality assessment
- `StructureDescriptors` (descriptors): Computed structural metrics
- `StructureSummary` (summary): Combined quality + descriptors
- `SearchQuery` (rcsb): RCSB search query builder
- `FileFormat` (rcsb): PDB/CIF format selection
- `AlignmentResult` (geometry): RMSD and transformation from alignment
- `PerResidueRmsd` (geometry): Per-residue RMSD for flexibility analysis
- `AtomSelection` (geometry): Atom selection for RMSD/alignment

## Common Patterns

### Parsing
```rust
// Auto-detect format
let structure = parse_structure_file("file.pdb")?;

// Explicit format
let structure = parse_pdb_file("file.pdb")?;
let structure = parse_mmcif_file("file.cif")?;

// From string
let structure = parse_pdb_string(content)?;

// Gzip-compressed files (feature: gzip)
let structure = parse_gzip_pdb_file("pdb1ubq.ent.gz")?;
let structure = parse_gzip_structure_file("structure.pdb.gz")?;  // Auto-detect
```

### Filtering (feature: filter)
```rust
// Fluent API with method chaining
let cleaned = structure
    .remove_ligands()
    .keep_only_chain("A")
    .keep_only_ca();

// In-place modifications
structure.normalize_chain_ids();
structure.reindex_residues();
structure.center_structure();
```

### Selection Language (feature: filter)
```rust
// PyMOL/VMD-style selection language
let selected = structure.select("chain A and name CA")?;
let backbone = structure.select("backbone and not hydrogen")?;
let residues = structure.select("resid 1:100 and protein")?;
let complex = structure.select("(chain A or chain B) and bfactor < 30.0")?;

// Validate without executing
PdbStructure::validate_selection("chain A and name CA")?;
```

**Selection syntax:**
- Basic: `chain A`, `name CA`, `resname ALA`, `resid 50`, `resid 1:100`, `element C`
- Keywords: `backbone`, `protein`, `nucleic`, `water`, `hetero`, `hydrogen`, `all`
- Numeric: `bfactor < 30.0`, `occupancy >= 0.5`
- Boolean: `and`, `or`, `not`, `()`

### Descriptors (feature: descriptors)
```rust
let rg = structure.radius_of_gyration();
let max_dist = structure.max_ca_distance();
let composition = structure.aa_composition();
let descriptors = structure.structure_descriptors(); // All at once
```

### Quality (feature: quality)
```rust
let report = structure.quality_report();
if report.is_analysis_ready() {
    // Single model, no altlocs, full atoms
}
```

### RCSB (feature: rcsb)
```rust
// Search
let query = SearchQuery::new()
    .with_text("kinase")
    .with_resolution_max(2.0);
let results = rcsb_search(&query, 10)?;

// Download
let structure = download_structure("1UBQ", FileFormat::Pdb)?;
```

### Geometry (feature: geometry)
```rust
use pdbrust::geometry::AtomSelection;

// RMSD calculation (without alignment)
let rmsd = structure1.rmsd_to(&structure2)?;  // CA atoms by default

// RMSD with different atom selection
let rmsd = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;

// Structure alignment (Kabsch algorithm)
let (aligned, result) = mobile.align_to(&target)?;
println!("RMSD: {:.4} Angstroms ({} atoms)", result.rmsd, result.num_atoms);

// Per-residue RMSD for flexibility analysis
let per_res = mobile.per_residue_rmsd_to(&target)?;
for r in &per_res {
    if r.rmsd > 2.0 {
        println!("Flexible: {}{} {}", r.residue_id.0, r.residue_id.1, r.rmsd);
    }
}
```

## Examples

The `examples/` directory contains runnable examples demonstrating common workflows:

| Example | Features Required | Description |
|---------|-------------------|-------------|
| `analysis_workflow.rs` | filter, descriptors, quality, summary | Complete pipeline: load → clean → analyze → export |
| `filtering_demo.rs` | filter | Fluent filtering API, method chaining, in-place modifications |
| `selection_demo.rs` | filter | PyMOL/VMD-style selection language: chain A and name CA |
| `geometry_demo.rs` | geometry | RMSD calculation, Kabsch alignment, per-residue RMSD |
| `rcsb_workflow.rs` | rcsb, descriptors | RCSB search queries, download, analyze (requires network) |
| `batch_processing.rs` | descriptors, summary | Process multiple files, compute summaries, export CSV |
| `full_pdb_benchmark.rs` | gzip, parallel, descriptors, quality, summary | Full PDB archive benchmark (230K structures) |
| `read_pdb.rs` | (none) | Basic PDB file reading and structure inspection |
| `write_pdb.rs` | (none) | Creating and writing PDB files |
| `basic_usage.rs` | (none) | Creating structures programmatically |
| `atom_interactive.rs` | (none) | Atom operations, distances, angles |
| `rust_benchmark.rs` | filter, descriptors | Performance benchmarking |

### Python Examples

The `pdbrust-python/examples/` directory contains Python examples:

| Example | Description |
|---------|-------------|
| `basic_usage.py` | Parsing, accessing atoms/residues, basic filtering |
| `writing_files.py` | Write PDB/mmCIF files, round-trip demonstration |
| `geometry_rmsd.py` | RMSD calculation, structure alignment, per-residue RMSD |
| `numpy_integration.py` | Coordinate arrays, distance matrices, contact maps |
| `rcsb_search.py` | RCSB search queries and structure downloads |

Run Python examples:
```bash
cd pdbrust-python/examples
python basic_usage.py
python geometry_rmsd.py
python numpy_integration.py
```

### Running Rust Examples
```bash
# Complete analysis workflow
cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"

# Filtering operations
cargo run --example filtering_demo --features "filter"

# Geometry: RMSD and alignment
cargo run --example geometry_demo --features "geometry"

# RCSB search and download
cargo run --example rcsb_workflow --features "rcsb,descriptors"

# Batch processing
cargo run --example batch_processing --features "descriptors,summary"

# Full PDB archive benchmark
cargo run --release --example full_pdb_benchmark \
    --features "gzip,parallel,descriptors,quality,summary" \
    -- /path/to/pdb/archive --output-dir ./benchmark_results

# Basic file reading
cargo run --example read_pdb -- examples/pdb_files/1UBQ.pdb
```

## File Structure Notes

- Test PDB files: `examples/pdb_files/`
- Examples: `examples/`
- Documentation: `docs/GETTING_STARTED.md`
- Benchmarks: `benchmarks/`
- MSRV: 1.70.0
- CI: Ubuntu, Windows, macOS + stable, MSRV

## Performance Notes

Rust vs Python benchmarks show:
- Parsing: 2-3x faster
- In-memory operations: 40-270x faster
- O(n²) operations (max_ca_distance): 260x faster

The speedup comes from:
1. Parse once, reuse structure (Python re-parses each call)
2. Zero-cost abstractions
3. No GIL
4. CPU cache locality

### Full PDB Archive Validation

PDBRust has been validated against the entire Protein Data Bank (230,655 structures):

| Metric | Value |
|--------|-------|
| Total Structures | 230,655 |
| Success Rate | 100% |
| Failed Parses | 0 |
| Total Atoms Parsed | 2,057,302,767 |
| Processing Rate | ~92 files/sec (128 threads) |
| Largest Structure | 2ku2 (1,290,100 atoms) |
| Smallest Structure | 5zmz (31 atoms) |

Results are stored in `benchmark_results/` directory.

## Python Bindings

The Python package is published to PyPI as `pdbrust`. It uses PyO3 for Rust-Python bindings and Maturin for building.

### Python Development Commands

```bash
# Install maturin
pip install maturin

# Build and install in development mode
cd pdbrust-python
maturin develop --release

# Build wheel for distribution
maturin build --release

# Publish to PyPI (requires MATURIN_PYPI_TOKEN or ~/.pypirc)
maturin publish --no-sdist
```

### Python Package Features

All features are enabled by default in the Python package:
- Parsing: PDB, mmCIF, gzip-compressed files
- Filtering: remove_ligands, keep_only_chain, keep_only_ca, etc.
- Descriptors: radius_of_gyration, max_ca_distance, aa_composition
- Quality: quality_report, has_altlocs, has_multiple_models
- RCSB: download_structure, rcsb_search with SearchQuery
- Numpy: get_coords_array, get_ca_coords_array (returns numpy.ndarray)

### Python API Pattern

```python
import pdbrust
import numpy as np

# Parse
structure = pdbrust.parse_pdb_file("protein.pdb")

# Filter (method chaining)
cleaned = structure.remove_ligands().keep_only_chain("A")

# Descriptors
rg = structure.radius_of_gyration()
desc = structure.structure_descriptors()

# Numpy arrays (efficient coordinate access)
coords = structure.get_coords_array()        # Shape: (N, 3)
ca_coords = structure.get_ca_coords_array()  # Shape: (CA, 3)

# RCSB
from pdbrust import download_structure, FileFormat, SearchQuery, rcsb_search
structure = download_structure("1UBQ", FileFormat.pdb())
results = rcsb_search(SearchQuery().with_text("kinase"), 10)
```

### CI/CD for Python

The `.github/workflows/python-publish.yml` workflow:
- Builds wheels for Linux (x86_64, aarch64), macOS (x86_64, aarch64), Windows (x64)
- Supports Python 3.9, 3.10, 3.11, 3.12
- Automatically publishes to PyPI on version tags (v*)
- Uses PyPI trusted publishing (no token needed in CI)

### Releasing a New Version

**IMPORTANT:** Always update CHANGELOG.md before releasing!

1. **Update CHANGELOG.md** with new version entry:
   - Move items from `[Unreleased]` to new version section
   - Follow [Keep a Changelog](https://keepachangelog.com/) format
   - Include: Added, Changed, Fixed, etc. sections as needed

2. **Update version numbers** in all files:
   - `Cargo.toml` (main library)
   - `pdbrust-python/Cargo.toml`
   - `pdbrust-python/pyproject.toml`
   - `pdbrust-python/python/pdbrust/__init__.py` (`__version__`)
   - `README.md` (version in dependency examples)

3. **Commit all changes**:
   ```bash
   git add -A
   git commit -m "chore: release vX.Y.Z"
   ```

4. **Wait for CI to pass** on the commit

5. **Create and push tag**:
   ```bash
   git tag vX.Y.Z
   git push origin main vX.Y.Z
   ```

6. **Publish to crates.io** (after GitHub Actions wheel builds pass):
   ```bash
   cargo publish
   ```

7. **Verify releases**:
   - Check https://crates.io/crates/pdbrust
   - Check https://pypi.org/project/pdbrust/
   - GitHub Actions will auto-publish to PyPI on tag push

### Zenodo Integration

PDBRust uses Zenodo for DOI-based academic citations. Once enabled, each GitHub release automatically gets archived with a unique DOI.

**Initial Setup (one-time):**
1. Go to https://zenodo.org and log in with GitHub
2. Navigate to Settings > GitHub
3. Enable the repository (HFooladi/pdbrust)
4. Create a GitHub release - Zenodo will automatically archive it

**After First Zenodo Release:**
1. Get your DOI from https://zenodo.org/account/settings/github/
2. Update `README.md`: uncomment and update the Zenodo badge with your DOI
3. Update `README.md`: uncomment and update the BibTeX citation with your DOI
4. Update `CITATION.cff`: add the `doi:` field with your DOI
5. Commit these changes

**Files for Zenodo:**
- `CITATION.cff` - Machine-readable citation metadata (GitHub recognizes this)
- `.zenodo.json` - Zenodo-specific metadata (keywords, related identifiers, etc.)

**Version-Specific vs Concept DOI:**
- Each release gets a unique version DOI (e.g., `10.5281/zenodo.1234567`)
- There's also a "concept DOI" that always resolves to the latest version
- Use the concept DOI in documentation for always-current citations
