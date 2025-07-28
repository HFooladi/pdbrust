use pdbrust::core::PdbStructure;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a new empty structure
    let mut structure = PdbStructure::new();

    // Set some basic information
    structure.header = Some("Example PDB Structure".to_string());
    structure.title = Some("Basic Usage Example".to_string());

    // Add some atoms
    structure.atoms.push(pdbrust::records::Atom {
        serial: 1,
        name: "CA".to_string(),
        alt_loc: None,
        residue_name: "ALA".to_string(),
        chain_id: "A".to_string(),
        residue_seq: 1,
        ins_code: None,
        x: 0.0,
        y: 0.0,
        z: 0.0,
        occupancy: 1.0,
        temp_factor: 20.0,
        element: "C".to_string(),
    });

    // Add a sequence record
    structure.seqres.push(pdbrust::records::SeqRes {
        serial: 1,
        chain_id: "A".to_string(),
        num_residues: 1,
        residues: vec!["ALA".to_string()],
    });

    // Add a remark
    structure.remarks.push(pdbrust::records::Remark {
        number: 2,
        content: "RESOLUTION. 2.0 ANGSTROMS.".to_string(),
    });

    // Get chain IDs
    let chain_ids = structure.get_chain_ids();
    println!("Found chains: {:?}", chain_ids);

    // Get residues for chain A
    let residues = structure.get_residues_for_chain("A");
    println!("Residues in chain A: {:?}", residues);

    // Get sequence for chain A
    let sequence = structure.get_sequence("A");
    println!("Sequence of chain A: {:?}", sequence);

    // Get remarks with number 2
    let resolution_remarks = structure.get_remarks_by_number(2);
    for remark in resolution_remarks {
        println!("Resolution remark: {}", remark.content);
    }

    // Translate the structure
    structure.translate(1.0, 2.0, 3.0);
    println!(
        "Translated first atom coordinates: ({}, {}, {})",
        structure.atoms[0].x, structure.atoms[0].y, structure.atoms[0].z
    );

    // Save the structure to a file
    structure.to_file("example_output.pdb")?;
    println!("Structure saved to example_output.pdb");

    Ok(())
}
