use pdbrust::{parse_pdb_file, Atom, Conect, PdbStructure, SSBond, SeqRes};
use std::fs::File;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_pdb(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(content.as_bytes()).unwrap();
    file
}

#[test]
fn test_parse_empty_file() {
    let file = create_test_pdb("");
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert!(structure.atoms.is_empty());
    assert!(structure.connects.is_empty());
    assert!(structure.ssbonds.is_empty());
    assert!(structure.seqres.is_empty());
}

#[test]
fn test_parse_single_atom() {
    let content =
        "ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    let atom = &structure.atoms[0];
    assert_eq!(atom.serial, 1);
    assert_eq!(atom.name, "N");
    assert_eq!(atom.residue_name, "ALA");
    assert_eq!(atom.chain_id, "A");
    assert_eq!(atom.residue_seq, 1);
    assert_eq!(atom.alt_loc, None);
    assert_eq!(atom.ins_code, None);
    assert_eq!(atom.occupancy, 1.00);
    assert_eq!(atom.temp_factor, 13.79);
    assert_eq!(atom.element, "N");
}

#[test]
fn test_parse_seqres() {
    let content = "SEQRES   1 A   21  MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.seqres.len(), 1);
    let seqres = &structure.seqres[0];
    assert_eq!(seqres.serial, 1);
    assert_eq!(seqres.chain_id, "A");
    assert_eq!(seqres.num_residues, 21);
    assert_eq!(seqres.residues.len(), 12);
    assert_eq!(seqres.residues[0], "MET");
    assert_eq!(seqres.residues[11], "THR");
}

#[test]
fn test_parse_conect() {
    let content = "CONECT    1    2    3    4\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.connects.len(), 1);
    let conect = &structure.connects[0];
    assert_eq!(conect.atom1, 1);
    assert_eq!(conect.atom2, 2);
    assert_eq!(conect.atom3, Some(3));
    assert_eq!(conect.atom4, Some(4));
}

#[test]
fn test_parse_ssbond() {
    let content = "SSBOND   1 CYS A    6    CYS A   11    1555   1555  2.04\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.ssbonds.len(), 1);
    let ssbond = &structure.ssbonds[0];
    assert_eq!(ssbond.serial, 1);
    assert_eq!(ssbond.residue1_name, "CYS");
    assert_eq!(ssbond.chain1_id, "A");
    assert_eq!(ssbond.residue1_seq, 6);
    assert_eq!(ssbond.icode1, None);
    assert_eq!(ssbond.residue2_name, "CYS");
    assert_eq!(ssbond.chain2_id, "A");
    assert_eq!(ssbond.residue2_seq, 11);
    assert_eq!(ssbond.icode2, None);
    assert_eq!(ssbond.sym1, 1555);
    assert_eq!(ssbond.sym2, 1555);
    assert_eq!(ssbond.length, 2.04);
}

#[test]
fn test_parse_complete_structure() {
    let content = "\
HEADER    PROTEIN                                01-JAN-01   1ABC
TITLE     EXAMPLE PROTEIN STRUCTURE
REMARK   1 THIS IS A TEST STRUCTURE
SEQRES   1 A   21  MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR
ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
CONECT    1    2
SSBOND   1 CYS A    6    CYS A   11    1555   1555  2.04
END
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(
        structure.header,
        Some("PROTEIN                                01-JAN-01   1ABC".to_string())
    );
    assert_eq!(
        structure.title,
        Some("EXAMPLE PROTEIN STRUCTURE".to_string())
    );
    assert_eq!(structure.remarks.len(), 1);
    assert_eq!(structure.remarks[0].content, "THIS IS A TEST STRUCTURE");
    assert_eq!(structure.seqres.len(), 1);
    assert_eq!(structure.atoms.len(), 2);
    assert_eq!(structure.connects.len(), 1);
    assert_eq!(structure.ssbonds.len(), 1);
}

// #[test]
// fn test_parse_multi_model() {
//     let content = "\
// MODEL     1
// ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
// ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
// ENDMDL
// MODEL     2
// ATOM      3  N   ALA A   1      27.147  14.199   3.215  1.00 13.79           N
// ATOM      4  CA  ALA A   1      26.147  13.199   2.215  1.00 12.79           C
// ENDMDL
// END
// ";
//     let file = create_test_pdb(content);
//     let result = parse_pdb_file(file.path());
//     assert!(result.is_ok());
//     let structure = result.unwrap();
//     assert_eq!(structure.models.len(), 2);
// 
//     let model1 = &structure.models[0];
//     assert_eq!(model1.serial, 1);
//     assert_eq!(model1.atoms.len(), 2);
//     assert_eq!(model1.atoms[0].serial, 1);
//     assert_eq!(model1.atoms[1].serial, 2);
// 
//     let model2 = &structure.models[1];
//     assert_eq!(model2.serial, 2);
//     assert_eq!(model2.atoms.len(), 2);
//     assert_eq!(model2.atoms[0].serial, 3);
//     assert_eq!(model2.atoms[1].serial, 4);
// 
//     // All atoms should be in the structure.atoms as well
//     assert_eq!(structure.atoms.len(), 4);
//     
//     // Check the coordinates to ensure the atoms are correctly assigned
//     assert!((structure.atoms[0].x - 27.047).abs() < 1e-6);
//     assert!((structure.atoms[2].x - 27.147).abs() < 1e-6);
// }
