# PDBRust

A comprehensive Rust library for parsing and analyzing PDB (Protein Data Bank) files. This library provides a robust and efficient way to work with protein structure data in PDB format.

## Features

- **Complete Record Support**
  - ATOM/HETATM records with full coordinate and metadata support
  - MODEL/ENDMDL for multi-model structures
  - SEQRES records for sequence information
  - CONECT records for connectivity data
  - SSBOND records for disulfide bonds
  - REMARK records with categorization
  - HEADER and TITLE metadata

- **Robust Error Handling**
  - Custom error types for different parsing scenarios
  - Detailed error messages for debugging
  - Safe handling of malformed PDB files

- **Utility Functions**
  - Chain identification and analysis
  - Residue sequence extraction
  - Connectivity analysis
  - Model-based structure organization
  - Disulfide bond analysis

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
pdbrust = "0.1.0"
```

## Usage

### Basic Structure Loading

```rust
use pdbrust::PdbStructure;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load a PDB file
    let structure = PdbStructure::from_file("protein.pdb")?;
    
    // Access basic information
    if let Some(header) = &structure.header {
        println!("Structure header: {}", header);
    }
    
    println!("Number of atoms: {}", structure.atoms.len());
    println!("Number of models: {}", structure.models.len());
    
    Ok(())
}
```

### Chain and Residue Analysis

```rust
use pdbrust::PdbStructure;

fn analyze_chains(structure: &PdbStructure) {
    // Get all chain IDs
    let chain_ids = structure.get_chain_ids();
    
    for chain_id in chain_ids {
        // Get residues in this chain
        let residues = structure.get_residues_for_chain(&chain_id);
        println!("Chain {} has {} residues", chain_id, residues.len());
        
        // Get sequence if available
        let sequence = structure.get_sequence(&chain_id);
        if !sequence.is_empty() {
            println!("Sequence: {}", sequence.join("-"));
        }
    }
}
```

### Connectivity Analysis

```rust
use pdbrust::PdbStructure;

fn analyze_connectivity(structure: &PdbStructure) {
    // Analyze disulfide bonds
    for bond in &structure.ssbonds {
        println!("Disulfide bond between:");
        println!("  Residue 1: {} {} {}", bond.residue1_name, bond.chain1_id, bond.residue1_seq);
        println!("  Residue 2: {} {} {}", bond.residue2_name, bond.chain2_id, bond.residue2_seq);
        if let Some(dist) = bond.distance {
            println!("  Distance: {:.2} Å", dist);
        }
    }
    
    // Analyze atom connectivity
    for atom in &structure.atoms {
        let connected = structure.get_connected_atoms(atom.serial);
        println!("Atom {} ({}) is connected to {} other atoms", 
            atom.serial, atom.name, connected.len());
    }
}
```

### Working with Models

```rust
use pdbrust::PdbStructure;

fn analyze_models(structure: &PdbStructure) {
    for model in &structure.models {
        println!("Model {} contains {} atoms", model.serial, model.atoms.len());
        
        // Access model-specific remarks
        for remark in &model.remarks {
            println!("Model-specific remark {}: {}", remark.number, remark.content);
        }
    }
}
```

## Error Handling

The library provides detailed error handling through the `PdbError` enum:

```rust
use pdbrust::{PdbStructure, PdbError};

fn load_structure(path: &str) -> Result<PdbStructure, PdbError> {
    let structure = PdbStructure::from_file(path)?;
    Ok(structure)
}
```

## Running Examples

The repository includes example programs demonstrating various features:

```bash
# Run the basic PDB reader example
cargo run --example read_pdb path/to/your/protein.pdb
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 