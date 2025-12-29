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
├── records/        # Record types (Atom, Residue, Chain, Model, SSBond, etc.)
├── writer/         # PDB file output
├── error.rs        # Error handling with thiserror
├── filter/         # [feature: filter] Filtering and cleaning operations
├── descriptors/    # [feature: descriptors] Structural descriptors
├── quality/        # [feature: quality] Quality assessment
├── summary/        # [feature: summary] Unified structure summaries
└── rcsb/           # [feature: rcsb] RCSB PDB search and download
```

### Feature Flags

| Feature | Description | Dependencies |
|---------|-------------|--------------|
| `filter` | Filtering, extraction, and cleaning operations | - |
| `descriptors` | Structural descriptors (Rg, composition, geometry) | - |
| `quality` | Quality assessment and reports | - |
| `summary` | Unified summaries combining quality + descriptors | `descriptors`, `quality` |
| `rcsb` | RCSB PDB search API and file download | `reqwest`, `serde`, `serde_json` |
| `parallel` | Parallel processing with Rayon | `rayon` |
| `geometry` | Geometric analysis with nalgebra | `nalgebra` |
| `analysis` | All analysis features | `filter`, `descriptors`, `quality`, `summary` |
| `full` | Everything | `parallel`, `geometry`, `analysis`, `rcsb` |

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

## File Structure Notes

- Test PDB files: `examples/pdb_files/`
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
