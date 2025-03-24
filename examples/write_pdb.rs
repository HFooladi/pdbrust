use pdbrust::{PdbStructure, Atom, SeqRes, Conect, SSBond, Remark};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Create a new PDB structure
    let mut structure = PdbStructure::new();

    // Add header and title information
    structure.header = Some("EXAMPLE STRUCTURE".to_string());
    structure.title = Some("DEMONSTRATION OF PDB FILE WRITING".to_string());

    // Add some remarks
    structure.remarks.push(Remark {
        number: 1,
        content: "THIS IS AN EXAMPLE STRUCTURE".to_string(),
    });

    // Add some atoms (creating a small peptide)
    let atoms = vec![
        Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "N".to_string(),
            ins_code: None,
        },
        Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.5,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
            ins_code: None,
        },
        Atom {
            serial: 3,
            name: "C".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 2.5,
            y: 1.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
            ins_code: None,
        },
    ];

    // Add atoms to structure
    structure.atoms.extend(atoms);

    // Add sequence information
    structure.seqres.push(SeqRes {
        serial: 1,
        chain_id: "A".to_string(),
        num_residues: 1,
        residues: vec!["ALA".to_string()],
    });

    // Add connectivity information
    structure.connects.push(Conect {
        atom_serial: 1,
        bonded_atoms: vec![2],
    });
    structure.connects.push(Conect {
        atom_serial: 2,
        bonded_atoms: vec![1, 3],
    });

    // Add a disulfide bond (even though this example doesn't have cysteines)
    structure.ssbonds.push(SSBond {
        serial: 1,
        residue1_name: "CYS".to_string(),
        chain1_id: "A".to_string(),
        residue1_seq: 1,
        residue2_name: "CYS".to_string(),
        chain2_id: "A".to_string(),
        residue2_seq: 2,
        distance: Some(2.05),
    });

    // Write the structure to a file
    structure.to_file("example_output.pdb")?;
    println!("Successfully wrote PDB file to example_output.pdb");

    // Demonstrate that we can read it back
    let read_structure = PdbStructure::from_file("example_output.pdb")?;
    println!("Successfully read back the PDB file");
    println!("Number of atoms: {}", read_structure.atoms.len());
    println!("Number of SEQRES records: {}", read_structure.seqres.len());
    println!("Number of CONECT records: {}", read_structure.connects.len());

    Ok(())
} 