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
└── rcsb/           # [feature: rcsb] RCSB PDB search and download

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
- `QualityReport` (quality): Structure quality assessment
- `StructureDescriptors` (descriptors): Computed structural metrics
- `StructureSummary` (summary): Combined quality + descriptors
- `SearchQuery` (rcsb): RCSB search query builder
- `FileFormat` (rcsb): PDB/CIF format selection

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

## Examples

The `examples/` directory contains runnable examples demonstrating common workflows:

| Example | Features Required | Description |
|---------|-------------------|-------------|
| `analysis_workflow.rs` | filter, descriptors, quality, summary | Complete pipeline: load → clean → analyze → export |
| `filtering_demo.rs` | filter | Fluent filtering API, method chaining, in-place modifications |
| `rcsb_workflow.rs` | rcsb, descriptors | RCSB search queries, download, analyze (requires network) |
| `batch_processing.rs` | descriptors, summary | Process multiple files, compute summaries, export CSV |
| `full_pdb_benchmark.rs` | gzip, parallel, descriptors, quality, summary | Full PDB archive benchmark (230K structures) |
| `read_pdb.rs` | (none) | Basic PDB file reading and structure inspection |
| `write_pdb.rs` | (none) | Creating and writing PDB files |
| `basic_usage.rs` | (none) | Creating structures programmatically |
| `atom_interactive.rs` | (none) | Atom operations, distances, angles |
| `rust_benchmark.rs` | filter, descriptors | Performance benchmarking |

### Running Examples
```bash
# Complete analysis workflow
cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"

# Filtering operations
cargo run --example filtering_demo --features "filter"

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
