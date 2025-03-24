use pdbrust::PdbStructure;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <pdb_file>", args[0]);
        std::process::exit(1);
    }

    let structure = PdbStructure::from_file(&args[1])?;
    
    // Print basic information about the structure
    if let Some(header) = &structure.header {
        println!("Header: {}", header);
    }
    if let Some(title) = &structure.title {
        println!("Title: {}", title);
    }
    
    println!("Number of atoms: {}", structure.atoms.len());
    
    // Print the first few atoms
    for (i, atom) in structure.atoms.iter().take(5).enumerate() {
        println!("\nAtom {}:", i + 1);
        println!("  Name: {} (Element: {})", atom.name, atom.element);
        println!("  Residue: {} {}", atom.residue_name, atom.residue_seq);
        println!("  Chain: {}", atom.chain_id);
        println!("  Position: ({:.3}, {:.3}, {:.3})", atom.x, atom.y, atom.z);
    }

    Ok(())
} 