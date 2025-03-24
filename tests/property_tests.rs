use proptest::prelude::*;
use pdbrust::{PdbStructure, PdbError};
use std::fs::File;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_pdb(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file
}

// Generate valid atom records with varying data
fn generate_atom_record(
    serial: i32,
    name: &str,
    residue_name: &str,
    chain_id: &str,
    residue_seq: i32,
    x: f64,
    y: f64,
    z: f64,
) -> String {
    format!(
        "ATOM  {:5} {:4} {:3} {:1}{:4}     {:8.3}{:8.3}{:8.3}  1.00  0.00           {:>2}  ",
        serial, name, residue_name, chain_id, residue_seq, x, y, z,
        name.chars().next().unwrap()
    )
}

proptest! {
    // Test that valid atom records are always parsed correctly
    #[test]
    fn test_valid_atom_records(
        serial in 1..99999,
        x in -999.999f64..999.999,
        y in -999.999f64..999.999,
        z in -999.999f64..999.999,
    ) {
        let atom_record = generate_atom_record(
            serial,
            "N",
            "ALA",
            "A",
            1,
            x,
            y,
            z,
        );
        
        let file = create_test_pdb(&atom_record);
        let result = PdbStructure::from_file(file.path());
        prop_assert!(result.is_ok());
        
        let structure = result.unwrap();
        prop_assert_eq!(structure.atoms.len(), 1);
        
        let atom = &structure.atoms[0];
        prop_assert_eq!(atom.serial, serial);
        prop_assert_eq!(&atom.name, "N");
        prop_assert_eq!(&atom.residue_name, "ALA");
        prop_assert_eq!(&atom.chain_id, "A");
        prop_assert_eq!(atom.residue_seq, 1);
        prop_assert!((atom.x - x).abs() < 1e-3);
        prop_assert!((atom.y - y).abs() < 1e-3);
        prop_assert!((atom.z - z).abs() < 1e-3);
    }

    // Test that chain IDs are always handled correctly
    #[test]
    fn test_chain_id_handling(
        chain_id in "[A-Z]",
        num_atoms in 1..10usize,
    ) {
        let mut content = String::new();
        for i in 0..num_atoms {
            content.push_str(&generate_atom_record(
                i as i32 + 1,
                "N",
                "ALA",
                &chain_id,
                1,
                0.0,
                0.0,
                0.0,
            ));
            content.push('\n');
        }
        
        let file = create_test_pdb(&content);
        let result = PdbStructure::from_file(file.path());
        prop_assert!(result.is_ok());
        
        let structure = result.unwrap();
        let chain_ids = structure.get_chain_ids();
        prop_assert_eq!(chain_ids.len(), 1);
        prop_assert!(chain_ids.contains(&chain_id.to_string()));
    }

    // Test that residue sequences are handled correctly
    #[test]
    fn test_residue_sequence(
        residue_seq in 1..9999,
        residue_name in prop::sample::select(vec!["ALA", "GLY", "VAL", "LEU", "ILE"]),
    ) {
        let atom_record = generate_atom_record(
            1,
            "N",
            residue_name,
            "A",
            residue_seq,
            0.0,
            0.0,
            0.0,
        );
        
        let file = create_test_pdb(&atom_record);
        let result = PdbStructure::from_file(file.path());
        prop_assert!(result.is_ok());
        
        let structure = result.unwrap();
        let residues = structure.get_residues_for_chain("A");
        prop_assert_eq!(residues.len(), 1);
        prop_assert_eq!(&residues[0], &(residue_seq, residue_name.to_string()));
    }

    // Test that multiple models are handled correctly
    #[test]
    fn test_multiple_models(
        num_models in 1..10usize,
        atoms_per_model in 1..5usize,
    ) {
        let mut content = String::new();
        for model in 1..=num_models {
            content.push_str(&format!("MODEL {:8}\n", model));
            for atom in 1..=atoms_per_model {
                content.push_str(&generate_atom_record(
                    atom as i32,
                    "N",
                    "ALA",
                    "A",
                    1,
                    0.0,
                    0.0,
                    0.0,
                ));
                content.push('\n');
            }
            content.push_str("ENDMDL\n");
        }
        
        let file = create_test_pdb(&content);
        let result = PdbStructure::from_file(file.path());
        prop_assert!(result.is_ok());
        
        let structure = result.unwrap();
        prop_assert_eq!(structure.models.len(), num_models);
        for model in &structure.models {
            prop_assert_eq!(model.atoms.len(), atoms_per_model);
        }
    }

    // Test that SEQRES records are handled correctly
    #[test]
    fn test_seqres_records(
        chain_id in "[A-Z]",
        num_residues in 1..20usize,
    ) {
        let residues = vec!["ALA", "GLY", "VAL", "LEU", "ILE", "PRO", "SER", "THR"];
        let selected_residues: Vec<&str> = (0..num_residues)
            .map(|i| residues[i % residues.len()])
            .collect();
        
        let content = format!(
            "SEQRES   1 {} {:3} {}\n",
            chain_id,
            num_residues,
            selected_residues.join(" ")
        );
        
        let file = create_test_pdb(&content);
        let result = PdbStructure::from_file(file.path());
        prop_assert!(result.is_ok());
        
        let structure = result.unwrap();
        let sequence = structure.get_sequence(&chain_id);
        prop_assert_eq!(sequence.len(), num_residues);
        for (i, residue) in sequence.iter().enumerate() {
            prop_assert_eq!(residue, selected_residues[i]);
        }
    }
} 