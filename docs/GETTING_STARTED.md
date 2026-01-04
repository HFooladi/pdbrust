# Getting Started with PDBRust

This guide will help you get up and running with PDBRust quickly.

## Installation

Add PDBRust to your `Cargo.toml`:

```toml
[dependencies]
pdbrust = "0.3"
```

### Choosing Features

PDBRust uses feature flags to keep the core library lightweight. Enable only what you need:

```toml
# Minimal: just parsing
pdbrust = "0.3"

# Common setup: parsing + filtering + analysis
pdbrust = { version = "0.3", features = ["filter", "descriptors", "quality"] }

# Full analysis suite
pdbrust = { version = "0.3", features = ["analysis"] }

# Everything including RCSB download
pdbrust = { version = "0.3", features = ["full"] }
```

### Which Features Do I Need?

| If you want to... | Enable these features |
|-------------------|----------------------|
| Just parse PDB/mmCIF files | (none - included by default) |
| Filter atoms, extract chains, clean structures | `filter` |
| Compute Rg, composition, distances | `descriptors` |
| Assess structure quality | `quality` |
| Get all metrics in one call | `summary` |
| Download from RCSB PDB | `rcsb` |
| Process files in parallel | `parallel` |
| All analysis features | `analysis` |
| Everything | `full` |

## Quick Start (5 minutes)

### 1. Parse a Structure

```rust
use pdbrust::{parse_pdb_file, parse_structure_file};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse PDB format explicitly
    let structure = parse_pdb_file("protein.pdb")?;

    // Or auto-detect format (works with .pdb, .cif, .ent)
    let structure = parse_structure_file("protein.cif")?;

    // Basic info
    println!("Atoms: {}", structure.atoms.len());
    println!("Chains: {:?}", structure.get_chain_ids());
    println!("Residues: {}", structure.get_residues().len());

    Ok(())
}
```

### 2. Filter and Clean (requires `filter` feature)

```rust
use pdbrust::parse_pdb_file;

let structure = parse_pdb_file("protein.pdb")?;

// Method chaining for clean code
let cleaned = structure
    .remove_ligands()      // Remove waters, ions, ligands
    .remove_hydrogens()    // Remove H atoms
    .keep_only_chain("A"); // Keep only chain A

// Extract CA coordinates for distance calculations
let ca_coords = structure.get_ca_coords(None);

// Get CA atoms as full Atom structs
let ca_atoms = structure.get_ca_atoms(None);
```

### 3. Compute Descriptors (requires `descriptors` feature)

```rust
let structure = parse_pdb_file("protein.pdb")?;

// Individual metrics
let rg = structure.radius_of_gyration();
let max_dist = structure.max_ca_distance();
let composition = structure.aa_composition();

// Or get everything at once
let descriptors = structure.structure_descriptors();
println!("Rg: {:.2} A", descriptors.radius_of_gyration);
println!("Hydrophobic: {:.1}%", descriptors.hydrophobic_ratio * 100.0);
```

### 4. Quality Assessment (requires `quality` feature)

```rust
let structure = parse_pdb_file("protein.pdb")?;

let report = structure.quality_report();
println!("Models: {}", report.num_models);
println!("Has altlocs: {}", report.has_altlocs);
println!("Has HETATM: {}", report.has_hetatm);

// Check if suitable for typical analysis
if report.is_analysis_ready() {
    println!("Ready for analysis!");
}
```

### 5. Download from RCSB (requires `rcsb` feature)

```rust
use pdbrust::rcsb::{download_structure, rcsb_search, SearchQuery, FileFormat};

// Download by PDB ID
let structure = download_structure("1UBQ", FileFormat::Pdb)?;

// Search RCSB
let query = SearchQuery::new()
    .with_text("kinase")
    .with_organism("Homo sapiens")
    .with_resolution_max(2.0);

let results = rcsb_search(&query, 10)?;
println!("Found {} structures", results.total_count);
```

## Common Workflows

### Workflow 1: Clean Structure for MD

```rust
let structure = parse_pdb_file("raw_structure.pdb")?;

let cleaned = structure
    .remove_ligands()      // Remove non-protein
    .remove_hydrogens()    // Will be added by MD software
    .keep_only_chain("A"); // Single chain

// Center at origin
let mut centered = cleaned.clone();
centered.center_structure();
centered.renumber_atoms();

write_pdb_file(&centered, "cleaned_for_md.pdb")?;
```

### Workflow 2: Dataset Characterization

```rust
use pdbrust::summary::batch_summarize;

let structures = vec![
    parse_pdb_file("protein1.pdb")?,
    parse_pdb_file("protein2.pdb")?,
    parse_pdb_file("protein3.pdb")?,
];

let summaries = batch_summarize(&structures);

for (i, summary) in summaries.iter().enumerate() {
    println!("Structure {}: {} residues, Rg={:.1}A",
        i+1, summary.num_residues, summary.radius_of_gyration);
}
```

### Workflow 3: RCSB Search Pipeline

```rust
use pdbrust::rcsb::{rcsb_search, download_structure, SearchQuery, FileFormat};

// Find human kinases with high resolution
let query = SearchQuery::new()
    .with_text("kinase")
    .with_organism("Homo sapiens")
    .with_resolution_max(2.0);

let results = rcsb_search(&query, 5)?;

for pdb_id in &results.pdb_ids {
    let structure = download_structure(pdb_id, FileFormat::Pdb)?;
    let rg = structure.radius_of_gyration();
    println!("{}: Rg = {:.1} A", pdb_id, rg);
}
```

## Running Examples

The repository includes several example files:

```bash
# Complete analysis workflow
cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"

# Filtering operations
cargo run --example filtering_demo --features "filter"

# RCSB search and download
cargo run --example rcsb_workflow --features "rcsb,descriptors"

# Batch processing
cargo run --example batch_processing --features "descriptors,summary"

# Basic file reading
cargo run --example read_pdb -- examples/pdb_files/1UBQ.pdb
```

## Error Handling

All parsing functions return `Result` with detailed error types:

```rust
use pdbrust::{parse_pdb_file, PdbError};

match parse_pdb_file("structure.pdb") {
    Ok(structure) => {
        println!("Loaded {} atoms", structure.atoms.len());
    }
    Err(PdbError::Io(e)) => {
        eprintln!("File error: {}", e);
    }
    Err(PdbError::InvalidRecord(msg)) => {
        eprintln!("Parse error: {}", msg);
    }
    Err(e) => {
        eprintln!("Other error: {}", e);
    }
}
```

## Performance Tips

1. **Parse once, reuse**: Parse the structure once and perform multiple analyses
2. **Filter early**: Apply filters before expensive computations
3. **Use iterators**: Prefer `.iter()` over collecting to vectors when possible
4. **Batch operations**: Use `batch_summarize` for multiple structures

```rust
// Good: filter first, then compute
let chain_a = structure.keep_only_chain("A");
let rg = chain_a.radius_of_gyration();

// Good: use iterators
let ca_count = structure.atoms.iter()
    .filter(|a| a.name.trim() == "CA")
    .count();
```

## Next Steps

- See the [API Documentation](https://docs.rs/pdbrust) for full details
- Check out the [examples/](../examples/) directory for more code
- Read the [guide module](https://docs.rs/pdbrust/latest/pdbrust/guide/) for comprehensive documentation

## Troubleshooting

### "feature X not found"

Make sure you've enabled the required feature in `Cargo.toml`:

```toml
pdbrust = { version = "0.3", features = ["filter", "descriptors"] }
```

### Network errors with RCSB

The `rcsb` feature requires internet access. Check your connection and firewall settings.

### Parse errors

- Ensure the file is a valid PDB or mmCIF format
- Check for file encoding issues (should be UTF-8 or ASCII)
- Try using `parse_structure_file()` for auto-detection
