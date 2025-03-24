use pdbrust::{PdbStructure, PdbError};
use std::fs::File;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_pdb(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file
}

#[test]
fn test_parse_empty_file() {
    let file = create_test_pdb("");
    let result = PdbStructure::from_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert!(structure.atoms.is_empty());
    assert!(structure.models.is_empty());
}

#[test]
fn test_parse_single_atom() {
    let content = "ATOM      1  N   ALA A   1      27.409  24.517   5.258  1.00 39.29           N  ";
    let file = create_test_pdb(content);
    let result = PdbStructure::from_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    
    let atom = &structure.atoms[0];
    assert_eq!(atom.serial, 1);
    assert_eq!(&atom.name, "N");
    assert_eq!(&atom.residue_name, "ALA");
    assert_eq!(&atom.chain_id, "A");
    assert_eq!(atom.residue_seq, 1);
    assert!((atom.x - 27.409).abs() < 1e-6);
    assert!((atom.y - 24.517).abs() < 1e-6);
    assert!((atom.z - 5.258).abs() < 1e-6);
    assert!((atom.occupancy - 1.00).abs() < 1e-6);
    assert!((atom.temp_factor - 39.29).abs() < 1e-6);
    assert_eq!(&atom.element, "N");
}

#[test]
fn test_parse_seqres() {
    let content = "SEQRES   1 A  20  MET ALA CYS PRO THR GLY LYS ALA SER VAL ALA ARG GLY LEU                  \n\
                  SEQRES   2 A  20  ALA LYS THR PRO GLY ASN                                                    ";
    let file = create_test_pdb(content);
    let result = PdbStructure::from_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.seqres.len(), 1);
    
    let seqres = &structure.seqres[0];
    assert_eq!(&seqres.chain_id, "A");
    assert_eq!(seqres.num_residues, 20);
    assert_eq!(seqres.residues.len(), 20);
    assert_eq!(&seqres.residues[0], "MET");
    assert_eq!(&seqres.residues[19], "ASN");
}

#[test]
fn test_parse_conect() {
    let content = "ATOM      1  N   CYS A   1      27.409  24.517   5.258  1.00 39.29           N  \n\
                  ATOM      2  CA  CYS A   1      26.897  23.161   5.048  1.00 39.29           C  \n\
                  CONECT    1    2\n\
                  CONECT    2    1";
    let file = create_test_pdb(content);
    let result = PdbStructure::from_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.connects.len(), 2);
    
    let connect1 = structure.connects.iter().find(|c| c.atom_serial == 1).unwrap();
    assert_eq!(connect1.bonded_atoms, vec![2]);
    
    let connect2 = structure.connects.iter().find(|c| c.atom_serial == 2).unwrap();
    assert_eq!(connect2.bonded_atoms, vec![1]);
}

#[test]
fn test_parse_ssbond() {
    let content = "SSBOND   1 CYS A   22    CYS A   44                          1555   1555  2.03  ";
    let file = create_test_pdb(content);
    let result = PdbStructure::from_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.ssbonds.len(), 1);
    
    let ssbond = &structure.ssbonds[0];
    assert_eq!(ssbond.serial, 1);
    assert_eq!(&ssbond.residue1_name, "CYS");
    assert_eq!(&ssbond.chain1_id, "A");
    assert_eq!(ssbond.residue1_seq, 22);
    assert_eq!(&ssbond.residue2_name, "CYS");
    assert_eq!(&ssbond.chain2_id, "A");
    assert_eq!(ssbond.residue2_seq, 44);
    assert!((ssbond.distance.unwrap() - 2.03).abs() < 1e-6);
}

#[test]
fn test_parse_multiple_models() {
    let content = "MODEL        1\n\
                  ATOM      1  N   ALA A   1      27.409  24.517   5.258  1.00 39.29           N  \n\
                  ENDMDL\n\
                  MODEL        2\n\
                  ATOM      1  N   ALA A   1      27.415  24.521   5.261  1.00 39.30           N  \n\
                  ENDMDL";
    let file = create_test_pdb(content);
    let result = PdbStructure::from_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.models.len(), 2);
    assert_eq!(structure.models[0].atoms.len(), 1);
    assert_eq!(structure.models[1].atoms.len(), 1);
}

#[test]
fn test_invalid_atom_record() {
    let content = "ATOM   INVALID DATA";
    let file = create_test_pdb(content);
    let result = PdbStructure::from_file(file.path());
    assert!(matches!(result, Err(pdbrust::PdbError::InvalidRecord(_))));
}

#[test]
fn test_get_chain_ids() {
    let content = "ATOM      1  N   ALA A   1      27.409  24.517   5.258  1.00 39.29           N  \n\
                  ATOM      2  N   ALA B   1      27.409  24.517   5.258  1.00 39.29           N  \n\
                  SEQRES   1 C  20  MET ALA CYS PRO THR GLY LYS ALA SER VAL                         ";
    let file = create_test_pdb(content);
    let structure = PdbStructure::from_file(file.path()).unwrap();
    let chain_ids = structure.get_chain_ids();
    assert_eq!(chain_ids.len(), 3);
    assert!(chain_ids.contains(&"A".to_string()));
    assert!(chain_ids.contains(&"B".to_string()));
    assert!(chain_ids.contains(&"C".to_string()));
}

#[test]
fn test_get_residues_for_chain() {
    let content = "ATOM      1  N   ALA A   1      27.409  24.517   5.258  1.00 39.29           N  \n\
                  ATOM      2  N   GLY A   2      27.409  24.517   5.258  1.00 39.29           N  \n\
                  ATOM      3  N   PRO B   1      27.409  24.517   5.258  1.00 39.29           N  ";
    let file = create_test_pdb(content);
    let structure = PdbStructure::from_file(file.path()).unwrap();
    
    let residues_a = structure.get_residues_for_chain("A");
    assert_eq!(residues_a.len(), 2);
    assert_eq!(&residues_a[0], &(1, "ALA".to_string()));
    assert_eq!(&residues_a[1], &(2, "GLY".to_string()));
    
    let residues_b = structure.get_residues_for_chain("B");
    assert_eq!(residues_b.len(), 1);
    assert_eq!(&residues_b[0], &(1, "PRO".to_string()));
}

#[test]
fn test_get_sequence() {
    let content = "SEQRES   1 A  20  MET ALA CYS PRO THR GLY LYS ALA SER VAL                         \n\
                  SEQRES   2 A  20  ALA ARG GLY LEU ALA LYS THR PRO GLY ASN                         ";
    let file = create_test_pdb(content);
    let structure = PdbStructure::from_file(file.path()).unwrap();
    
    let sequence = structure.get_sequence("A");
    assert_eq!(sequence.len(), 20);
    assert_eq!(&sequence[0], "MET");
    assert_eq!(&sequence[19], "ASN");
}

#[test]
fn test_get_connected_atoms() {
    let content = "ATOM      1  N   CYS A   1      27.409  24.517   5.258  1.00 39.29           N  \n\
                  ATOM      2  CA  CYS A   1      26.897  23.161   5.048  1.00 39.29           C  \n\
                  ATOM      3  C   CYS A   1      27.847  22.059   5.504  1.00 39.29           C  \n\
                  CONECT    2    1    3\n\
                  CONECT    1    2\n\
                  CONECT    3    2";
    let file = create_test_pdb(content);
    let structure = PdbStructure::from_file(file.path()).unwrap();
    
    let connected = structure.get_connected_atoms(2);
    assert_eq!(connected.len(), 2);
    assert!(connected.iter().any(|a| a.serial == 1));
    assert!(connected.iter().any(|a| a.serial == 3));
} 