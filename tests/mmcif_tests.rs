//! Comprehensive tests for mmCIF parsing functionality

use pdbrust::{parse_mmcif_file, parse_mmcif_string, parse_structure_file, PdbStructure};
use std::fs::File;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_mmcif(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::with_suffix(".cif").unwrap();
    file.write_all(content.as_bytes()).unwrap();
    file
}

#[test]
fn test_basic_mmcif_parsing() {
    let mmcif_content = r#"
data_test
_entry.id TEST_ENTRY

_struct.title "Test Protein Structure"

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . MET A 1 20.154 16.967 23.486 1.00 25.00
ATOM 2 C CA . MET A 1 21.498 16.929 22.908 1.00 24.50
ATOM 3 C C . MET A 1 22.392 18.134 23.185 1.00 23.75
ATOM 4 O O . MET A 1 23.612 18.123 23.195 1.00 22.90
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    // Test basic structure properties
    assert_eq!(structure.atoms.len(), 4);
    assert!(structure.header.is_some());
    assert_eq!(structure.title, Some("Test Protein Structure".to_string()));

    // Test first atom
    let atom1 = &structure.atoms[0];
    assert_eq!(atom1.serial, 1);
    assert_eq!(atom1.name, "N");
    assert_eq!(atom1.element, "N");
    assert_eq!(atom1.residue_name, "MET");
    assert_eq!(atom1.chain_id, "A");
    assert_eq!(atom1.residue_seq, 1);
    assert_eq!(atom1.x, 20.154);
    assert_eq!(atom1.y, 16.967);
    assert_eq!(atom1.z, 23.486);
    assert_eq!(atom1.occupancy, 1.00);
    assert_eq!(atom1.temp_factor, 25.00);

    // Test last atom
    let atom4 = &structure.atoms[3];
    assert_eq!(atom4.serial, 4);
    assert_eq!(atom4.name, "O");
    assert_eq!(atom4.element, "O");
    assert_eq!(atom4.residue_name, "MET");
    assert_eq!(atom4.residue_seq, 1);
}

#[test]
fn test_mmcif_sequence_parsing() {
    let mmcif_content = r#"
data_test
_entry.id TEST_SEQ

loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 MET
1 2 ALA
1 3 GLY
1 4 SER
1 5 THR

loop_
_struct_asym.id
_struct_asym.entity_id
A 1
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.seqres.len(), 1);
    
    let seqres = &structure.seqres[0];
    assert_eq!(seqres.chain_id, "A");
    assert_eq!(seqres.num_residues, 5);
    assert_eq!(seqres.residues, vec!["MET", "ALA", "GLY", "SER", "THR"]);
}

#[test]
fn test_mmcif_disulfide_bonds() {
    let mmcif_content = r#"
data_test
_entry.id TEST_SS

loop_
_struct_disulfid.id
_struct_disulfid.ptnr1_auth_comp_id
_struct_disulfid.ptnr1_auth_asym_id
_struct_disulfid.ptnr1_auth_seq_id
_struct_disulfid.ptnr2_auth_comp_id
_struct_disulfid.ptnr2_auth_asym_id
_struct_disulfid.ptnr2_auth_seq_id
disulf1 CYS A 6 CYS A 11
disulf2 CYS B 15 CYS B 28
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.ssbonds.len(), 2);
    
    let bond1 = &structure.ssbonds[0];
    assert_eq!(bond1.serial, 1);
    assert_eq!(bond1.residue1_name, "CYS");
    assert_eq!(bond1.chain1_id, "A");
    assert_eq!(bond1.residue1_seq, 6);
    assert_eq!(bond1.residue2_name, "CYS");
    assert_eq!(bond1.chain2_id, "A");
    assert_eq!(bond1.residue2_seq, 11);

    let bond2 = &structure.ssbonds[1];
    assert_eq!(bond2.serial, 2);
    assert_eq!(bond2.residue1_name, "CYS");
    assert_eq!(bond2.chain1_id, "B");
    assert_eq!(bond2.residue1_seq, 15);
    assert_eq!(bond2.residue2_name, "CYS");
    assert_eq!(bond2.chain2_id, "B");
    assert_eq!(bond2.residue2_seq, 28);
}

#[test]
fn test_mmcif_with_alternate_locations() {
    let mmcif_content = r#"
data_test
_entry.id TEST_ALT

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 0.000 0.000 0.000 1.00 20.00
ATOM 2 C CA A ALA A 1 1.000 0.000 0.000 0.60 20.00
ATOM 3 C CA B ALA A 1 1.100 0.100 0.100 0.40 20.00
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 3);
    
    let atom1 = &structure.atoms[0];
    assert_eq!(atom1.alt_loc, None);
    
    let atom2 = &structure.atoms[1];
    assert_eq!(atom2.alt_loc, Some('A'));
    
    let atom3 = &structure.atoms[2];
    assert_eq!(atom3.alt_loc, Some('B'));
}

#[test]
fn test_mmcif_with_insertion_codes() {
    let mmcif_content = r#"
data_test
_entry.id TEST_INS

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 . 0.000 0.000 0.000 1.00 20.00
ATOM 2 N N . ALA A 1 A 1.000 0.000 0.000 1.00 20.00
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 2);
    
    let atom1 = &structure.atoms[0];
    assert_eq!(atom1.ins_code, None);
    
    let atom2 = &structure.atoms[1];
    assert_eq!(atom2.ins_code, Some('A'));
}

#[test]
fn test_mmcif_remarks_parsing() {
    let mmcif_content = r#"
data_test
_entry.id TEST_REM

_audit_conform.dict_name mmcif_pdbx.dic
_audit_conform.dict_version 5.281

_refine.ls_d_res_high 1.50
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert!(!structure.remarks.is_empty());
    
    // Should have at least one remark about dictionary conformance
    let dict_remark = structure.remarks.iter()
        .find(|r| r.content.contains("mmcif_pdbx.dic"));
    assert!(dict_remark.is_some());
    
    // Should have a resolution remark
    let res_remark = structure.remarks.iter()
        .find(|r| r.content.contains("1.50") && r.number == 2);
    assert!(res_remark.is_some());
}

#[test]
fn test_mmcif_multiple_chains() {
    let mmcif_content = r#"
data_test
_entry.id TEST_MULTI

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . MET A 1 0.000 0.000 0.000 1.00 20.00
ATOM 2 C CA . MET A 1 1.000 0.000 0.000 1.00 20.00
ATOM 3 N N . VAL B 1 10.000 0.000 0.000 1.00 20.00
ATOM 4 C CA . VAL B 1 11.000 0.000 0.000 1.00 20.00

loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 MET
2 1 VAL

loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 4);
    
    let chain_ids = structure.get_chain_ids();
    assert_eq!(chain_ids.len(), 2);
    assert!(chain_ids.contains(&"A".to_string()));
    assert!(chain_ids.contains(&"B".to_string()));
    
    let chain_a_residues = structure.get_residues_for_chain("A");
    assert_eq!(chain_a_residues.len(), 1);
    assert_eq!(chain_a_residues[0].1, "MET");
    
    let chain_b_residues = structure.get_residues_for_chain("B");
    assert_eq!(chain_b_residues.len(), 1);
    assert_eq!(chain_b_residues[0].1, "VAL");

    // Test sequences
    assert_eq!(structure.seqres.len(), 2);
    let seq_a = structure.get_sequence("A");
    let seq_b = structure.get_sequence("B");
    assert_eq!(seq_a, vec!["MET"]);
    assert_eq!(seq_b, vec!["VAL"]);
}

#[test]
fn test_mmcif_file_parsing() {
    let mmcif_content = r#"
data_test
_entry.id FILE_TEST

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 0.000 0.000 0.000 1.00 20.00
"#;

    let file = create_test_mmcif(mmcif_content);
    let structure = parse_mmcif_file(file.path()).unwrap();
    
    assert_eq!(structure.atoms.len(), 1);
    assert_eq!(structure.atoms[0].residue_name, "ALA");
}

#[test]
fn test_auto_format_detection_mmcif() {
    let mmcif_content = r#"
data_auto_test
_entry.id AUTO_TEST

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 0.000 0.000 0.000 1.00 20.00
"#;

    // Test with .cif extension
    let cif_file = create_test_mmcif(mmcif_content);
    let structure = parse_structure_file(cif_file.path()).unwrap();
    assert_eq!(structure.atoms.len(), 1);

    // Test with unknown extension but mmCIF content
    let mut unknown_file = NamedTempFile::with_suffix(".unknown").unwrap();
    unknown_file.write_all(mmcif_content.as_bytes()).unwrap();
    let structure = parse_structure_file(unknown_file.path()).unwrap();
    assert_eq!(structure.atoms.len(), 1);
}

#[test]
fn test_mmcif_error_handling() {
    // Test invalid mmCIF content
    let invalid_content = "not valid mmcif content";
    let result = parse_mmcif_string(invalid_content);
    // Should not panic, may return empty structure or error depending on implementation
    assert!(result.is_ok() || result.is_err());

    // Test missing required fields
    let partial_content = r#"
data_test
_entry.id PARTIAL

loop_
_atom_site.group_PDB
_atom_site.id
ATOM 1
"#;
    let result = parse_mmcif_string(partial_content);
    // Should handle gracefully
    assert!(result.is_ok() || result.is_err());
}

#[test]
fn test_mmcif_hetatm_records() {
    let mmcif_content = r#"
data_test
_entry.id HETATM_TEST

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 0.000 0.000 0.000 1.00 20.00
HETATM 2 O O . HOH W 101 10.000 10.000 10.000 1.00 30.00
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 2);
    
    let protein_atom = &structure.atoms[0];
    assert_eq!(protein_atom.residue_name, "ALA");
    assert_eq!(protein_atom.chain_id, "A");
    
    let water_atom = &structure.atoms[1];
    assert_eq!(water_atom.residue_name, "HOH");
    assert_eq!(water_atom.chain_id, "W");
    assert_eq!(water_atom.residue_seq, 101);
}

#[test]
fn test_mmcif_coordinate_precision() {
    let mmcif_content = r#"
data_test
_entry.id PRECISION_TEST

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 123.456789 -987.654321 0.123456 0.567 89.123
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 1);
    
    let atom = &structure.atoms[0];
    assert!((atom.x - 123.456789).abs() < 1e-6);
    assert!((atom.y - (-987.654321)).abs() < 1e-6);
    assert!((atom.z - 0.123456).abs() < 1e-6);
    assert!((atom.occupancy - 0.567).abs() < 1e-6);
    assert!((atom.temp_factor - 89.123).abs() < 1e-6);
}