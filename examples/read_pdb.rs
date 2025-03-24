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
    
    // Print model information
    println!("\nModels: {}", structure.models.len());
    for model in &structure.models {
        println!("  Model {}: {} atoms", model.serial, model.atoms.len());
    }
    
    // Print chain information
    let chain_ids = structure.get_chain_ids();
    println!("\nChains: {}", chain_ids.len());
    for chain_id in &chain_ids {
        let residues = structure.get_residues_for_chain(chain_id);
        println!("  Chain {}: {} residues", chain_id, residues.len());
        
        // Print sequence if available
        let sequence = structure.get_sequence(chain_id);
        if !sequence.is_empty() {
            println!("    Sequence: {} residues", sequence.len());
            // Print first 10 residues if available
            if sequence.len() > 10 {
                println!("    First 10 residues: {}", sequence[0..10].join("-"));
            } else {
                println!("    Sequence: {}", sequence.join("-"));
            }
        }
    }
    
    // Print disulfide bonds
    println!("\nDisulfide bonds: {}", structure.ssbonds.len());
    for (i, bond) in structure.ssbonds.iter().enumerate() {
        println!("  Bond {}: {} {} {} -- {} {} {}", 
            i + 1,
            bond.residue1_name, bond.chain1_id, bond.residue1_seq,
            bond.residue2_name, bond.chain2_id, bond.residue2_seq
        );
        if let Some(dist) = bond.distance {
            println!("    Distance: {:.2} Å", dist);
        }
    }
    
    // Print some remarks (if any)
    if !structure.remarks.is_empty() {
        println!("\nRemarks:");
        
        // Get unique remark numbers
        let remark_numbers: Vec<i32> = {
            let mut numbers: Vec<i32> = structure.remarks.iter()
                .map(|r| r.number)
                .collect();
            numbers.sort();
            numbers.dedup();
            numbers
        };
        
        for num in remark_numbers.iter().take(5) {
            let remarks = structure.get_remarks_by_number(*num);
            println!("  REMARK {}: {} entries", num, remarks.len());
            
            // Print first remark of this type
            if !remarks.is_empty() {
                println!("    First entry: {}", remarks[0].content);
            }
        }
        
        if remark_numbers.len() > 5 {
            println!("  ... and {} more remark types", remark_numbers.len() - 5);
        }
    }
    
    // Print connectivity information for the first few atoms
    println!("\nConnectivity examples:");
    for i in 0..5 {
        if i < structure.atoms.len() {
            let atom = &structure.atoms[i];
            let connected = structure.get_connected_atoms(atom.serial);
            println!("  Atom {} ({}): connected to {} atoms", 
                atom.serial, atom.name, connected.len());
            
            for (j, conn_atom) in connected.iter().take(3).enumerate() {
                println!("    {}: Atom {} ({})", 
                    j + 1, conn_atom.serial, conn_atom.name);
            }
            
            if connected.len() > 3 {
                println!("    ... and {} more atoms", connected.len() - 3);
            }
        }
    }

    Ok(())
}