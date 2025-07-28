# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PDBRust is a Rust library for parsing and analyzing Protein Data Bank (PDB) and mmCIF files. It provides comprehensive support for molecular structure data with robust error handling and utility functions for structural analysis.

### Core Architecture

- **Dual Format Support**: The library supports both PDB (legacy text format) and mmCIF (modern dictionary format) with automatic format detection
- **Unified Data Model**: Both formats are converted to a common `PdbStructure` representation for consistent API usage
- **Modular Design**: Clean separation between parsing, data models, and utility functions

### Key Modules

- `src/core/`: Core data structures and format-specific implementations
  - `pdb.rs`: PDB format handling and `PdbStructure` definition
  - `mmcif.rs`: mmCIF format parsing and data structures
  - `mmcif_converter.rs`: Converts mmCIF data to unified PdbStructure format
- `src/parser/`: Format-specific parsers with auto-detection
- `src/records/`: Individual record type definitions (Atom, Residue, Chain, etc.)
- `src/writer/`: Output functionality for writing structure data
- `src/error.rs`: Comprehensive error handling with `thiserror`

### Format Detection Logic

The library auto-detects format based on file content, not just extension. Use `parse_structure_file()` for automatic detection or format-specific functions like `parse_pdb_file()` and `parse_mmcif_file()`.

## Development Commands

### Building and Testing
```bash
# Build the project
cargo build

# Run all tests
cargo test

# Run tests with all features enabled
cargo test --all-features

# Run benchmarks
cargo bench
```

### Code Quality
```bash
# Format code
cargo fmt

# Run clippy lints (CI uses -D warnings flag)
cargo clippy --all-targets --all-features

# Check documentation
cargo doc --no-deps --all-features

# Security audit
cargo audit
```

### Features
The project uses optional features defined in Cargo.toml:
- `parallel`: Enables Rayon-based parallel processing
- `geometry`: Adds nalgebra-based geometric analysis functions

Enable features during development:
```bash
cargo test --features "parallel,geometry"
cargo build --all-features
```

## Testing Strategy

- Unit tests in individual modules
- Integration tests in `tests/` directory
- Property-based testing using `proptest` for robust validation
- Benchmark tests using `criterion` for performance tracking
- CI runs tests on multiple platforms (Ubuntu, Windows, macOS) and Rust versions (stable, MSRV 1.70.0)

## File Structure Notes

- Test data is located in `examples/pdb_files/` with sample PDB and mmCIF files
- Examples in `examples/` demonstrate common usage patterns
- The library maintains MSRV (Minimum Supported Rust Version) of 1.70.0

## Error Handling

The codebase uses `thiserror` for structured error handling. All parsing functions return `Result<T, PdbError>` with detailed error context for debugging malformed files.

## Key Data Structures

- `PdbStructure`: Main structure containing all parsed data
- `Atom`: Individual atom records with coordinates and metadata
- `Residue`: Amino acid/nucleotide residue information
- `Chain`: Protein/nucleic acid chain organization
- `Model`: Support for multi-model structures (NMR ensembles)
- `SSBond`: Disulfide bond connectivity information