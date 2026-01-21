[![Crates.io](https://img.shields.io/crates/v/pdbrust.svg)](https://crates.io/crates/pdbrust)
[![PyPI](https://img.shields.io/pypi/v/pdbrust.svg)](https://pypi.org/project/pdbrust/)
[![Documentation](https://docs.rs/pdbrust/badge.svg)](https://docs.rs/pdbrust)
[![Rust CI/CD](https://github.com/hfooladi/pdbrust/actions/workflows/rust.yml/badge.svg)](https://github.com/hfooladi/pdbrust/actions/workflows/rust.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18232203.svg)](https://doi.org/10.5281/zenodo.18232203)

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
pdbrust = "0.5"
```

With optional features:

```toml
[dependencies]
pdbrust = { version = "0.5", features = ["filter", "descriptors", "rcsb", "gzip"] }
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
| `geometry` | RMSD calculation, structure alignment (Kabsch), per-residue RMSD |
| `dssp` | DSSP 4-like secondary structure assignment (H, G, I, P, E, B, T, S, C) |
| `rcsb` | Search and download structures from RCSB PDB |
| `rcsb-async` | Async/concurrent bulk downloads with rate limiting |
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

### Geometry: RMSD and Alignment

```rust
use pdbrust::{parse_pdb_file, geometry::AtomSelection};

let structure1 = parse_pdb_file("model1.pdb")?;
let structure2 = parse_pdb_file("model2.pdb")?;

// Calculate RMSD (without alignment)
let rmsd = structure1.rmsd_to(&structure2)?;
println!("RMSD: {:.3} Å", rmsd);

// Align structures using Kabsch algorithm
let (aligned, result) = structure1.align_to(&structure2)?;
println!("Alignment RMSD: {:.3} Å ({} atoms)", result.rmsd, result.num_atoms);

// Per-residue RMSD for flexibility analysis
let per_res = structure1.per_residue_rmsd_to(&structure2)?;
for r in per_res.iter().filter(|r| r.rmsd > 2.0) {
    println!("Flexible: {}{} {:.2} Å", r.chain_id, r.residue_seq, r.rmsd);
}

// Different atom selections
let rmsd_bb = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;
let rmsd_all = structure1.rmsd_to_with_selection(&structure2, AtomSelection::AllAtoms)?;
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

### Bulk Downloads with Async

```rust
use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};

#[tokio::main]
async fn main() {
    let pdb_ids = vec!["1UBQ", "8HM2", "4INS", "1HHB", "2MBP"];

    // Download with default options (5 concurrent, 100ms rate limit)
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;

    // Or with custom options
    let options = AsyncDownloadOptions::default()
        .with_max_concurrent(10)
        .with_rate_limit_ms(50);
    let results = download_multiple_async(&pdb_ids, FileFormat::Cif, Some(options)).await;

    for (pdb_id, result) in results {
        match result {
            Ok(structure) => println!("{}: {} atoms", pdb_id, structure.atoms.len()),
            Err(e) => eprintln!("{}: {}", pdb_id, e),
        }
    }
}
```

**Python:**

```python
from pdbrust import download_multiple, AsyncDownloadOptions, FileFormat

# Download multiple structures concurrently
results = download_multiple(["1UBQ", "8HM2", "4INS"], FileFormat.pdb())

# With custom options
options = AsyncDownloadOptions(max_concurrent=10, rate_limit_ms=50)
results = download_multiple(pdb_ids, FileFormat.cif(), options)

for r in results:
    if r.success:
        print(f"{r.pdb_id}: {len(r.get_structure().atoms)} atoms")
    else:
        print(f"{r.pdb_id}: {r.error}")
```

### Secondary Structure Assignment (DSSP)

```rust
use pdbrust::parse_pdb_file;

let structure = parse_pdb_file("protein.pdb")?;

// Compute DSSP-like secondary structure
let ss = structure.assign_secondary_structure();
println!("Helix: {:.1}%", ss.helix_fraction * 100.0);
println!("Sheet: {:.1}%", ss.sheet_fraction * 100.0);
println!("Coil:  {:.1}%", ss.coil_fraction * 100.0);

// Get as compact string (e.g., "HHHHEEEECCCC")
let ss_string = structure.secondary_structure_string();

// Get composition tuple
let (helix, sheet, coil) = structure.secondary_structure_composition();
```

**Python:**

```python
import pdbrust

structure = pdbrust.parse_pdb_file("protein.pdb")

# Get secondary structure assignment
ss = structure.assign_secondary_structure()
print(f"Helix: {ss.helix_fraction*100:.1f}%")
print(f"Sheet: {ss.sheet_fraction*100:.1f}%")

# Compact string representation
ss_string = structure.secondary_structure_string()
print(f"SS: {ss_string}")  # e.g., "CCCCHHHHHHHCCEEEEEECCC"

# Iterate over residue assignments
for res in ss:
    print(f"{res.chain_id}{res.residue_seq}: {res.ss.code()}")
```

## Common Workflows

See the [examples/](examples/) directory for complete working code:

| Workflow | Example | Features Used |
|----------|---------|---------------|
| Load, clean, analyze, export | [analysis_workflow.rs](examples/analysis_workflow.rs) | filter, descriptors, quality, summary |
| Filter and clean structures | [filtering_demo.rs](examples/filtering_demo.rs) | filter |
| RMSD and structure alignment | [geometry_demo.rs](examples/geometry_demo.rs) | geometry |
| Search and download from RCSB | [rcsb_workflow.rs](examples/rcsb_workflow.rs) | rcsb, descriptors |
| Async bulk downloads | [async_download_demo.rs](examples/async_download_demo.rs) | rcsb-async, descriptors |
| Process multiple files | [batch_processing.rs](examples/batch_processing.rs) | descriptors, summary |

**Python examples** are available in [pdbrust-python/examples/](pdbrust-python/examples/):
- `basic_usage.py` - Parsing and structure access
- `writing_files.py` - Write PDB/mmCIF files
- `geometry_rmsd.py` - RMSD and alignment
- `numpy_integration.py` - Numpy arrays, distance matrices, contact maps
- `rcsb_search.py` - RCSB search and download

Run Rust examples with:

```bash
cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"
cargo run --example filtering_demo --features "filter"
cargo run --example geometry_demo --features "geometry"
cargo run --example rcsb_workflow --features "rcsb,descriptors"
cargo run --example async_download_demo --features "rcsb-async,descriptors"
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

### Platform Notes

The Python package includes full functionality on **macOS** and **Windows**.
On **Linux**, the RCSB download/search features are not available in the pre-built wheels
due to cross-compilation constraints. All other features (parsing, filtering, analysis,
geometry, numpy arrays, etc.) work on all platforms.

| Platform | Parsing | Filtering | Descriptors | Geometry | DSSP | RCSB |
|----------|---------|-----------|-------------|----------|------|------|
| macOS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Windows | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Linux | ✓ | ✓ | ✓ | ✓ | ✓ | - |

See [pdbrust-python/README.md](pdbrust-python/README.md) for full Python API documentation.

## Documentation

- [API Documentation (Rust)](https://docs.rs/pdbrust)
- [PyPI Package](https://pypi.org/project/pdbrust/)
- [Examples](examples/)

## Citation

If you use PDBRust in your research, please cite it using the metadata in our [CITATION.cff](CITATION.cff) file:

```bibtex
@software{pdbrust,
  author = {Fooladi, Hosein},
  title = {PDBRust: A High-Performance Rust Library for PDB/mmCIF Parsing and Analysis},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18232203},
  url = {https://doi.org/10.5281/zenodo.18232203},
  version = {0.5.0}
}
```

Or in text format:

> Fooladi, H. (2025). PDBRust: A High-Performance Rust Library for PDB/mmCIF Parsing and Analysis. Zenodo. https://doi.org/10.5281/zenodo.18232203

## License

MIT
