use pdbrust::{parse_pdb_file, PdbStructure, Atom, SeqRes, Conect, SSBond};
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
    assert!(structure.header.is_empty());
    assert!(structure.title.is_empty());
    assert!(structure.remarks.is_empty());
    assert!(structure.seqres.is_empty());
    assert!(structure.atoms.is_empty());
    assert!(structure.conect.is_empty());
    assert!(structure.ssbond.is_empty());
}

#[test]
fn test_parse_single_atom() {
    let content = "ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    let atom = &structure.atoms[0];
    assert_eq!(atom.serial, 1);
    assert_eq!(atom.name, "N");
    assert_eq!(atom.res_name, "ALA");
    assert_eq!(atom.chain_id, 'A');
    assert_eq!(atom.res_seq, 1);
    assert_eq!(atom.x, 27.047);
    assert_eq!(atom.y, 14.099);
    assert_eq!(atom.z, 3.115);
    assert_eq!(atom.occupancy, 1.00);
    assert_eq!(atom.temp_factor, 13.79);
    assert_eq!(atom.element, "N");
    assert_eq!(atom.alt_loc, ' ');
    assert_eq!(atom.i_code, ' ');
    assert_eq!(atom.segment, "");
    assert_eq!(atom.charge, "");
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
    assert_eq!(seqres.ser_num, 1);
    assert_eq!(seqres.chain_id, 'A');
    assert_eq!(seqres.num_res, 21);
    assert_eq!(seqres.res_names.len(), 12);
    assert_eq!(seqres.res_names[0], "MET");
    assert_eq!(seqres.res_names[11], "THR");
}

#[test]
fn test_parse_conect() {
    let content = "CONECT    1    2    3    4\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.conect.len(), 1);
    let conect = &structure.conect[0];
    assert_eq!(conect.serial1, 1);
    assert_eq!(conect.serial2, 2);
    assert_eq!(conect.serial3, Some(3));
    assert_eq!(conect.serial4, Some(4));
}

#[test]
fn test_parse_ssbond() {
    let content = "SSBOND   1 CYS A    6    CYS A   11    1555   1555  2.04\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.ssbond.len(), 1);
    let ssbond = &structure.ssbond[0];
    assert_eq!(ssbond.ser_num, 1);
    assert_eq!(ssbond.res1, "CYS");
    assert_eq!(ssbond.chain_id1, 'A');
    assert_eq!(ssbond.seq1, 6);
    assert_eq!(ssbond.icode1, ' ');
    assert_eq!(ssbond.res2, "CYS");
    assert_eq!(ssbond.chain_id2, 'A');
    assert_eq!(ssbond.seq2, 11);
    assert_eq!(ssbond.icode2, ' ');
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
    assert_eq!(structure.header, "HEADER    PROTEIN                                01-JAN-01   1ABC");
    assert_eq!(structure.title, "TITLE     EXAMPLE PROTEIN STRUCTURE");
    assert_eq!(structure.remarks.len(), 1);
    assert_eq!(structure.remarks[0], "REMARK   1 THIS IS A TEST STRUCTURE");
    assert_eq!(structure.seqres.len(), 1);
    assert_eq!(structure.atoms.len(), 2);
    assert_eq!(structure.conect.len(), 1);
    assert_eq!(structure.ssbond.len(), 1);
} 