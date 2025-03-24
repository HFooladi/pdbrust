use pdbrust::PdbStructure;
use std::env;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Get the PDB file path from command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <pdb_file>", args[0]);
        std::process::exit(1);
    }

    let pdb_file = &args[1];
    println!("Reading PDB file: {}", pdb_file);

    // Read the PDB file
    let structure = PdbStructure::from_file(pdb_file)?;

    // Print basic information about the structure
    if let Some(header) = &structure.header {
        println!("Header: {}", header);
    }
    if let Some(title) = &structure.title {
        println!("Title: {}", title);
    }

    // Print statistics
    println!("\nStructure Statistics:");
    println!("Number of atoms: {}", structure.atoms.len());
    println!("Number of models: {}", structure.models.len());
    
    // Print chain information
    let chains = structure.get_chain_ids();
    println!("\nChains in structure: {}", chains.len());
    for chain_id in &chains {
        let residues = structure.get_residues_for_chain(chain_id);
        let sequence = structure.get_sequence(chain_id);
        println!("Chain {}: {} residues, {} residues in SEQRES", 
                chain_id, residues.len(), sequence.len());
    }

    // Print connectivity information
    println!("\nConnectivity Information:");
    println!("Number of CONECT records: {}", structure.connects.len());
    println!("Number of disulfide bonds: {}", structure.ssbonds.len());

    // Print some remarks if present
    if !structure.remarks.is_empty() {
        println!("\nFirst few remarks:");
        for remark in structure.remarks.iter().take(3) {
            println!("REMARK {}: {}", remark.number, remark.content);
        }
    }

    Ok(())
}