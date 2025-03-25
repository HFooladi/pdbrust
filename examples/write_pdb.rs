use pdbrust::{PdbStructure, Atom, SeqRes, Conect, SSBond, Remark, write_pdb_file};
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
            x: 2.0,
            y: 1.5,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
            ins_code: None,
        },
        Atom {
            serial: 4,
            name: "O".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 3.0,
            y: 1.5,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "O".to_string(),
            ins_code: None,
        },
    ];
    structure.atoms.extend(atoms);

    // Add SEQRES records
    structure.seqres.push(SeqRes {
        serial: 1,
        chain_id: "A".to_string(),
        num_residues: 1,
        residues: vec!["ALA".to_string()],
    });

    // Add CONECT records
    structure.connects.push(Conect {
        atom1: 1,
        atom2: 2,
        atom3: 0,
        atom4: 0,
    });
    structure.connects.push(Conect {
        atom1: 2,
        atom2: 3,
        atom3: 0,
        atom4: 0,
    });
    structure.connects.push(Conect {
        atom1: 3,
        atom2: 4,
        atom3: 0,
        atom4: 0,
    });

    // Add an SSBOND record (example)
    structure.ssbonds.push(SSBond {
        ser_num: 1,
        res1: "CYS".to_string(),
        chain1: "A".to_string(),
        resseq1: 1,
        icode1: " ".to_string(),
        res2: "CYS".to_string(),
        chain2: "A".to_string(),
        resseq2: 2,
        icode2: " ".to_string(),
        sym1: "1_555".to_string(),
        sym2: "1_555".to_string(),
        length: 2.03,
    });

    // Write the structure to a file using the new write_pdb_file function
    write_pdb_file(&structure, "example.pdb")?;
    println!("Successfully wrote example.pdb");

    // Demonstrate some structure operations
    println!("\nStructure Information:");
    println!("Number of atoms: {}", structure.atoms.len());
    println!("Number of chains: {}", structure.get_chain_ids().len());
    
    // Example of translating the structure
    println!("\nTranslating structure by (1.0, 1.0, 1.0)...");
    structure.translate(1.0, 1.0, 1.0);
    
    // Write the translated structure
    write_pdb_file(&structure, "example_translated.pdb")?;
    println!("Successfully wrote example_translated.pdb");

    Ok(())
} 