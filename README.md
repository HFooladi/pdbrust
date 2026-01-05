[![Crates.io](https://img.shields.io/crates/v/pdbrust.svg)](https://crates.io/crates/pdbrust)
[![PyPI](https://img.shields.io/pypi/v/pdbrust.svg)](https://pypi.org/project/pdbrust/)
[![Documentation](https://docs.rs/pdbrust/badge.svg)](https://docs.rs/pdbrust)
[![Rust CI/CD](https://github.com/hfooladi/pdbrust/actions/workflows/rust.yml/badge.svg)](https://github.com/hfooladi/pdbrust/actions/workflows/rust.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# PDBRust

A fast Rust library for parsing and analyzing PDB and mmCIF protein structure files. Also available as a Python package with **40-260x speedups** over pure Python implementations.

## Installation

### Python

```bash
pip install pdbrust
```

### Rust

```toml
[dependencies]
pdbrust = "0.3"
```

With optional features:

```toml
[dependencies]
pdbrust = { version = "0.3", features = ["filter", "descriptors", "rcsb", "gzip"] }
```

## Quick Start

### Python

```python
import pdbrust

# Parse a PDB file
structure = pdbrust.parse_pdb_file("protein.pdb")
print(f"Atoms: {structure.num_atoms}")
print(f"Chains: {structure.get_chain_ids()}")

# Filter and analyze
cleaned = structure.remove_ligands().keep_only_chain("A")
rg = cleaned.radius_of_gyration()
print(f"Radius of gyration: {rg:.2f} Å")

# Get coordinates as numpy arrays (fast!)
import numpy as np
coords = structure.get_coords_array()      # Shape: (N, 3)
ca_coords = structure.get_ca_coords_array() # Shape: (CA, 3)

# Download from RCSB PDB
from pdbrust import download_structure, FileFormat
structure = download_structure("1UBQ", FileFormat.pdb())
```

### Rust

```rust
use pdbrust::{parse_pdb_file, PdbStructure};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let structure = parse_pdb_file("protein.pdb")?;

    println!("Atoms: {}", structure.atoms.len());
    println!("Chains: {:?}", structure.get_chain_ids());

    Ok(())
}
```

## Features

| Feature | Description |
|---------|-------------|
| `filter` | Filter atoms, extract chains, remove ligands, clean structures |
| `descriptors` | Radius of gyration, amino acid composition, geometric metrics |
| `quality` | Structure quality assessment (altlocs, missing residues, etc.) |
| `summary` | Combined quality + descriptors in one call |
| `rcsb` | Search and download structures from RCSB PDB |
| `gzip` | Parse gzip-compressed files (.ent.gz, .pdb.gz, .cif.gz) |
| `parallel` | Parallel processing with Rayon |

## Examples

### Filter and Clean Structures

```rust
use pdbrust::parse_pdb_file;

let structure = parse_pdb_file("protein.pdb")?;

// Extract CA coordinates
let ca_coords = structure.get_ca_coords(None);

// Chain operations with fluent API
let chain_a = structure
    .remove_ligands()
    .keep_only_chain("A")
    .keep_only_ca();
```

### Compute Structural Descriptors

```rust
let structure = parse_pdb_file("protein.pdb")?;

let rg = structure.radius_of_gyration();
let max_dist = structure.max_ca_distance();
let composition = structure.aa_composition();

// Or get everything at once
let descriptors = structure.structure_descriptors();
```

### Parse Gzip-Compressed Files

```rust
use pdbrust::parse_gzip_pdb_file;

// Parse gzip-compressed PDB files from the PDB archive
let structure = parse_gzip_pdb_file("pdb1ubq.ent.gz")?;
println!("Atoms: {}", structure.atoms.len());
```

### Download from RCSB PDB

```rust
use pdbrust::rcsb::{download_structure, rcsb_search, SearchQuery, FileFormat};

// Download a structure
let structure = download_structure("1UBQ", FileFormat::Pdb)?;

// Search RCSB
let query = SearchQuery::new()
    .with_text("kinase")
    .with_organism("Homo sapiens")
    .with_resolution_max(2.0);

let results = rcsb_search(&query, 10)?;
```

## Common Workflows

See the [examples/](examples/) directory for complete working code:

| Workflow | Example | Features Used |
|----------|---------|---------------|
| Load, clean, analyze, export | [analysis_workflow.rs](examples/analysis_workflow.rs) | filter, descriptors, quality, summary |
| Filter and clean structures | [filtering_demo.rs](examples/filtering_demo.rs) | filter |
| Search and download from RCSB | [rcsb_workflow.rs](examples/rcsb_workflow.rs) | rcsb, descriptors |
| Process multiple files | [batch_processing.rs](examples/batch_processing.rs) | descriptors, summary |

Run examples with:

```bash
cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"
cargo run --example filtering_demo --features "filter"
cargo run --example rcsb_workflow --features "rcsb,descriptors"
cargo run --example batch_processing --features "descriptors,summary"
```

For a complete getting started guide, see [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md).

## Performance

Benchmarks against equivalent Python code show **40-260x speedups** for in-memory operations:

| Operation | Speedup |
|-----------|---------|
| Parsing | 2-3x |
| get_ca_coords | 240x |
| max_ca_distance | 260x |
| radius_of_gyration | 100x |

### Full PDB Archive Validation

PDBRust has been validated against the **entire Protein Data Bank**:

| Metric | Value |
|--------|-------|
| Total Structures Tested | 230,655 |
| Success Rate | **100%** |
| Failed Parses | 0 |
| Total Atoms Parsed | 2,057,302,767 |
| Processing Rate | ~92 files/sec |
| Largest Structure | 2ku2 (1,290,100 atoms) |

Run the full benchmark yourself:

```bash
cargo run --release --example full_pdb_benchmark \
    --features "gzip,parallel,descriptors,quality,summary" \
    -- /path/to/pdb/archive --output-dir ./results
```

## Python Package

Pre-built wheels available for Linux, macOS, and Windows (Python 3.9-3.12):

```bash
pip install pdbrust
```

See [pdbrust-python/README.md](pdbrust-python/README.md) for full Python API documentation.

## Documentation

- [API Documentation (Rust)](https://docs.rs/pdbrust)
- [PyPI Package](https://pypi.org/project/pdbrust/)
- [Examples](examples/)

## Citation

If you use PDBRust in your research, please cite:

```bibtex
@software{pdbrust,
  author = {Fooladi, Hosein},
  title = {PDBRust: A High-Performance Rust Library for PDB/mmCIF Parsing and Analysis},
  year = {2025},
  url = {https://github.com/HFooladi/pdbrust},
  version = {0.3.0}
}
```

Or in text format:

> Fooladi, H. (2025). PDBRust: A High-Performance Rust Library for PDB/mmCIF Parsing and Analysis. https://github.com/HFooladi/pdbrust

## License

MIT
