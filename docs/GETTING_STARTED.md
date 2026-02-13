# Getting Started with PDBRust

This guide will help you get up and running with PDBRust quickly.

## Installation

Add PDBRust to your `Cargo.toml`:

```toml
[dependencies]
pdbrust = "0.7"
```

### Choosing Features

PDBRust uses feature flags to keep the core library lightweight. Enable only what you need:

```toml
# Minimal: just parsing
pdbrust = "0.7"

# Common setup: parsing + filtering + analysis
pdbrust = { version = "0.7", features = ["filter", "descriptors", "quality"] }

# Full analysis suite
pdbrust = { version = "0.7", features = ["analysis"] }

# Everything including RCSB download
pdbrust = { version = "0.7", features = ["full"] }
```

### Which Features Do I Need?

| If you want to... | Enable these features |
|-------------------|----------------------|
| Just parse PDB/mmCIF files | (none - included by default) |
| Filter atoms, extract chains, clean structures | `filter` |
| Use selection language (chain A and name CA) | `filter` |
| Compute Rg, composition, B-factor analysis | `descriptors` |
| Validate ligand pose geometry (clashes, overlap) | `ligand-quality` |
| Assess protein-protein docking quality (DockQ) | `dockq` |
| Assess structure quality | `quality` |
| Get all metrics in one call | `summary` |
| Calculate RMSD and align structures | `geometry` |
| Compute secondary structure (DSSP) | `dssp` |
| Download from RCSB PDB | `rcsb` |
| Async bulk downloads with rate limiting | `rcsb-async` |
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

### 2b. Selection Language (requires `filter` feature)

```rust
let structure = parse_pdb_file("protein.pdb")?;

// PyMOL/VMD-style selections
let chain_a = structure.select("chain A")?;
let ca_atoms = structure.select("name CA")?;
let backbone = structure.select("backbone")?;

// Combine with boolean operators
let chain_a_ca = structure.select("chain A and name CA")?;
let heavy_atoms = structure.select("protein and not hydrogen")?;
let complex = structure.select("(chain A or chain B) and bfactor < 30.0")?;

// Residue ranges and numeric comparisons
let active_site = structure.select("resid 50:60")?;
let flexible = structure.select("bfactor > 40.0")?;
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

### 3b. B-factor Analysis (requires `descriptors` feature)

```rust
let structure = parse_pdb_file("protein.pdb")?;

// B-factor statistics
let mean_b = structure.b_factor_mean();
let mean_ca = structure.b_factor_mean_ca();
let std_b = structure.b_factor_std();
println!("Mean B-factor: {:.2} Å²", mean_b);

// Per-residue B-factor profile
let profile = structure.b_factor_profile();
for res in &profile {
    println!("{}{}: mean={:.2}", res.chain_id, res.residue_seq, res.mean);
}

// Identify flexible/rigid regions
let flexible = structure.flexible_residues(50.0);  // B > 50 Å²
let rigid = structure.rigid_residues(15.0);        // B < 15 Å²

// Normalize for cross-structure comparison
let normalized = structure.normalize_b_factors();
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

### 5b. Async Bulk Downloads (requires `rcsb-async` feature)

```rust
use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};

#[tokio::main]
async fn main() {
    let pdb_ids = vec!["1UBQ", "8HM2", "4INS", "1HHB"];

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

### 6. Geometry: RMSD and Alignment (requires `geometry` feature)

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
    println!("Flexible: {} {:.2} Å", r.residue_name, r.rmsd);
}

// Different atom selections
let rmsd_bb = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;
```

### 7. Secondary Structure (requires `dssp` feature)

```rust
let structure = parse_pdb_file("protein.pdb")?;

// Compute DSSP-like secondary structure
let ss = structure.assign_secondary_structure();
println!("Helix: {:.1}%", ss.helix_fraction * 100.0);
println!("Sheet: {:.1}%", ss.sheet_fraction * 100.0);
println!("Coil:  {:.1}%", ss.coil_fraction * 100.0);

// Get compact string (e.g., "HHHHEEEECCCC")
let ss_string = structure.secondary_structure_string();

// Get composition tuple
let (helix, sheet, coil) = structure.secondary_structure_composition();

// Iterate over per-residue assignments
for res in &ss.residue_assignments {
    println!("{}{}: {} ({})",
        res.chain_id, res.residue_seq, res.residue_name, res.ss.code());
}
```

### 8. Ligand Pose Quality (requires `ligand-quality` feature)

```rust
let structure = parse_pdb_file("protein_ligand.pdb")?;

// List all ligands in the structure
let ligands = structure.get_ligand_names();
println!("Found ligands: {:?}", ligands);

// Validate a specific ligand
if let Some(report) = structure.ligand_pose_quality("LIG") {
    println!("Ligand: {} ({}{})", report.ligand_name,
             report.ligand_chain_id, report.ligand_residue_seq);
    println!("Atoms: {}", report.ligand_atom_count);
    println!("Min distance to protein: {:.2} Å", report.min_protein_ligand_distance);
    println!("Protein clashes: {}", report.num_clashes);
    println!("Volume overlap: {:.1}%", report.protein_volume_overlap_pct);

    if report.is_geometry_valid {
        println!("✓ Pose passes geometry checks");
    } else {
        println!("✗ Pose fails geometry checks");
        // Show clash details
        for clash in report.clashes.iter().take(3) {
            println!("  Clash: {} {} - {} {}: {:.2}Å (expected >{:.2}Å)",
                clash.protein_residue_name, clash.protein_atom_name,
                clash.ligand_atom_name, clash.ligand_element,
                clash.distance, clash.expected_min_distance);
        }
    }
}

// Validate all ligands at once
let reports = structure.all_ligand_pose_quality();
for report in &reports {
    let status = if report.is_geometry_valid { "PASS" } else { "FAIL" };
    println!("{}: {}", report.ligand_name, status);
}
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

# Selection language
cargo run --example selection_demo --features "filter"

# B-factor analysis
cargo run --example b_factor_demo --features "descriptors"

# Secondary structure (DSSP)
cargo run --example secondary_structure_demo --features "dssp"

# Ligand pose quality (PoseBusters-style checks)
cargo run --example ligand_quality_demo --features "ligand-quality"

# RMSD and structure alignment
cargo run --example geometry_demo --features "geometry"

# RCSB search and download
cargo run --example rcsb_workflow --features "rcsb,descriptors"

# Async bulk downloads
cargo run --example async_download_demo --features "rcsb-async,descriptors"

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
pdbrust = { version = "0.7", features = ["filter", "descriptors"] }
```

### Network errors with RCSB

The `rcsb` feature requires internet access. Check your connection and firewall settings.

### Parse errors

- Ensure the file is a valid PDB or mmCIF format
- Check for file encoding issues (should be UTF-8 or ASCII)
- Try using `parse_structure_file()` for auto-detection
