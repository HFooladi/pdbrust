# PDBRust

A Rust library for parsing PDB (Protein Data Bank) files. This library provides a simple and efficient way to read and work with protein structure data in PDB format.

## Features

- Parse ATOM and HETATM records from PDB files
- Extract header and title information
- Support for all standard PDB atom fields
- Error handling for malformed PDB files
- Simple and intuitive API

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
pdbrust = "0.1.0"
```

### Example

```rust
use pdbrust::PdbStructure;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load a PDB file
    let structure = PdbStructure::from_file("protein.pdb")?;
    
    // Access structure information
    println!("Number of atoms: {}", structure.atoms.len());
    
    // Iterate through atoms
    for atom in &structure.atoms {
        println!("Atom {} of residue {} {}", 
            atom.name,
            atom.residue_name,
            atom.residue_seq
        );
    }
    
    Ok(())
}
```

## Running the Example

The repository includes an example program that demonstrates how to use the library. To run it:

```bash
cargo run --example read_pdb path/to/your/protein.pdb
```

## License

This project is licensed under the MIT License. 