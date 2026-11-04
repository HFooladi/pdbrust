//! Comprehensive tests for mmCIF parsing and writing functionality

use pdbrust::{
    Atom, PdbStructure, parse_mmcif_file, parse_mmcif_string, parse_pdb_file, parse_structure_file,
    write_mmcif_file, write_mmcif_string,
};
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

/// The disulfide pairs of crambin (1CRN) with the distances in its RCSB mmCIF.
const CRAMBIN_DISULFIDES: [(i32, i32, f64); 3] = [(3, 40, 2.004), (4, 32, 2.035), (16, 26, 2.051)];

fn assert_crambin_disulfides(structure: &PdbStructure) {
    let pairs: Vec<_> = structure
        .ssbonds
        .iter()
        .map(|b| (b.residue1_seq, b.residue2_seq, b.length))
        .collect();
    assert_eq!(pairs.len(), 3, "expected 3 disulfides, got {pairs:?}");
    for ((seq1, seq2, length), (exp1, exp2, exp_len)) in pairs.iter().zip(CRAMBIN_DISULFIDES) {
        assert_eq!((*seq1, *seq2), (exp1, exp2));
        assert!(
            (length - exp_len).abs() < 0.01,
            "length {length} != {exp_len}"
        );
    }
    for (i, bond) in structure.ssbonds.iter().enumerate() {
        assert_eq!(bond.serial, i as i32 + 1);
        assert_eq!(
            (bond.chain1_id.as_str(), bond.chain2_id.as_str()),
            ("A", "A")
        );
        assert_eq!(
            (bond.residue1_name.as_str(), bond.residue2_name.as_str()),
            ("CYS", "CYS")
        );
        assert_eq!((bond.sym1, bond.sym2), (1555, 1555));
    }
}

/// RCSB mmCIF files record disulfides in `_struct_conn` (conn_type_id `disulf`),
/// not in `_struct_disulfid`.
#[test]
fn test_mmcif_disulfides_from_struct_conn_rcsb_file() {
    let structure = parse_mmcif_file("examples/pdb_files/1CRN.cif").unwrap();
    assert_crambin_disulfides(&structure);
}

#[test]
fn test_mmcif_struct_conn_keeps_only_disulfides() {
    let mmcif_content = r#"
data_test
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.ptnr2_symmetry
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.pdbx_dist_value
covale1 covale A NAG 1 ? 1_555 B ASN 5 ? 1_555 C NAG 501 B ASN 45 1.44
disulf1 disulf B CYS 2 ? 1_555 B CYS 9 A 1_555 B CYS 22 B CYS 30 2.03
metalc1 metalc B CYS 3 ? 1_555 D ZN 1 ? 1_555 B CYS 23 B ZN 601 2.30
disulf2 DISULF B CYS 4 ? 1_555 E CYS 4 ? 2_655 B CYS 24 E CYS 24 2.07
"#;
    let structure = parse_mmcif_string(mmcif_content).unwrap();
    assert_eq!(structure.ssbonds.len(), 2);

    // Label numbering wins, matching how `_atom_site` rows are converted.
    let bond1 = &structure.ssbonds[0];
    assert_eq!(bond1.serial, 1);
    assert_eq!((bond1.chain1_id.as_str(), bond1.residue1_seq), ("B", 2));
    assert_eq!((bond1.chain2_id.as_str(), bond1.residue2_seq), ("B", 9));
    assert_eq!((bond1.icode1, bond1.icode2), (None, Some('A')));
    assert!((bond1.length - 2.03).abs() < 1e-9);

    let bond2 = &structure.ssbonds[1];
    assert_eq!(bond2.serial, 2);
    assert_eq!((bond2.chain2_id.as_str(), bond2.residue2_seq), ("E", 4));
    assert_eq!((bond2.sym1, bond2.sym2), (1555, 2655));
    assert!((bond2.length - 2.07).abs() < 1e-9);
}

#[test]
fn test_mmcif_disulfides_survive_write_read_round_trip() {
    let structure = parse_mmcif_file("examples/pdb_files/1CRN.cif").unwrap();
    let written = write_mmcif_string(&structure).unwrap();
    let reread = parse_mmcif_string(&written).unwrap();
    assert_crambin_disulfides(&reread);
}

#[test]
fn test_mmcif_disulfides_listed_in_both_categories_are_not_duplicated() {
    let mmcif_content = r#"
data_test
loop_
_struct_disulfid.id
_struct_disulfid.ptnr1_auth_comp_id
_struct_disulfid.ptnr1_auth_asym_id
_struct_disulfid.ptnr1_auth_seq_id
_struct_disulfid.ptnr2_auth_comp_id
_struct_disulfid.ptnr2_auth_asym_id
_struct_disulfid.ptnr2_auth_seq_id
disulf1 CYS A 6 CYS A 11
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.pdbx_dist_value
disulf1 disulf A CYS 11 A CYS 6 2.04
disulf2 disulf A CYS 20 A CYS 30 2.05
"#;
    let structure = parse_mmcif_string(mmcif_content).unwrap();
    let pairs: Vec<_> = structure
        .ssbonds
        .iter()
        .map(|b| (b.serial, b.residue1_seq, b.residue2_seq))
        .collect();
    assert_eq!(pairs, vec![(1, 6, 11), (2, 20, 30)]);
}

#[test]
fn test_mmcif_struct_disulfid_reads_distance_not_symmetry() {
    let mmcif_content = r#"
data_test
loop_
_struct_disulfid.id
_struct_disulfid.ptnr1_auth_comp_id
_struct_disulfid.ptnr1_auth_asym_id
_struct_disulfid.ptnr1_auth_seq_id
_struct_disulfid.ptnr1_symmetry
_struct_disulfid.ptnr2_auth_comp_id
_struct_disulfid.ptnr2_auth_asym_id
_struct_disulfid.ptnr2_auth_seq_id
_struct_disulfid.pdbx_dist_value
disulf1 CYS A 6 1_555 CYS A 11 2.10
"#;
    let structure = parse_mmcif_string(mmcif_content).unwrap();
    assert_eq!(structure.ssbonds.len(), 1);
    assert!((structure.ssbonds[0].length - 2.10).abs() < 1e-9);
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
    let dict_remark = structure
        .remarks
        .iter()
        .find(|r| r.content.contains("mmcif_pdbx.dic"));
    assert!(dict_remark.is_some());

    // Should have a resolution remark
    let res_remark = structure
        .remarks
        .iter()
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
    // Test with auth_seq_id fallback for non-polymer atoms (label_seq_id = ".")
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
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 1 0.000 0.000 0.000 1.00 20.00
HETATM 2 O O . HOH W . 101 10.000 10.000 10.000 1.00 30.00
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 2);

    let protein_atom = &structure.atoms[0];
    assert_eq!(protein_atom.residue_name, "ALA");
    assert_eq!(protein_atom.chain_id, "A");
    assert_eq!(protein_atom.residue_seq, 1);

    // Water uses auth_seq_id since label_seq_id is "."
    let water_atom = &structure.atoms[1];
    assert_eq!(water_atom.residue_name, "HOH");
    assert_eq!(water_atom.chain_id, "W");
    assert_eq!(water_atom.residue_seq, 101);
}

#[test]
fn test_mmcif_hetatm_dot_label_seq_id() {
    // Test that HETATM records with "." label_seq_id use auth_seq_id correctly
    let mmcif_content = r#"
data_test
_entry.id DOT_SEQ_TEST

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
HETATM 1 ZN ZN . ZN B . 401 -21.378 29.266 -1.191 1.00 24.04
HETATM 2 S S . SO4 D . 403 -45.892 41.335 -2.826 1.00 37.98
"#;

    let structure = parse_mmcif_string(mmcif_content).unwrap();

    assert_eq!(structure.atoms.len(), 2);

    let zn_atom = &structure.atoms[0];
    assert_eq!(zn_atom.residue_name, "ZN");
    assert_eq!(zn_atom.chain_id, "B");
    assert_eq!(zn_atom.residue_seq, 401);

    let so4_atom = &structure.atoms[1];
    assert_eq!(so4_atom.residue_name, "SO4");
    assert_eq!(so4_atom.chain_id, "D");
    assert_eq!(so4_atom.residue_seq, 403);
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

// ============================================================================
// mmCIF Writing Tests
// ============================================================================

#[test]
fn test_write_mmcif_basic() {
    // Create a simple structure programmatically
    let mut structure = PdbStructure::new();
    structure.header = Some("TEST STRUCTURE".to_string());
    structure.title = Some("Test Protein Structure".to_string());

    structure.atoms.push(Atom::new(
        1,
        "N".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        0.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "N".to_string(),
        None,
    ));
    structure.atoms.push(Atom::new(
        2,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        1.458,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    ));

    // Write to string
    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Verify basic structure
    assert!(mmcif_string.contains("data_"));
    assert!(mmcif_string.contains("_entry.id"));
    assert!(mmcif_string.contains("_struct.title"));
    assert!(mmcif_string.contains("loop_"));
    assert!(mmcif_string.contains("_atom_site."));
    assert!(mmcif_string.contains("ATOM"));
    assert!(mmcif_string.contains("ALA"));
}

#[test]
fn test_write_mmcif_file() {
    let mut structure = PdbStructure::new();
    structure.atoms.push(Atom::new(
        1,
        "CA".to_string(),
        None,
        "GLY".to_string(),
        "A".to_string(),
        1,
        10.5,
        20.3,
        30.1,
        1.0,
        25.0,
        "C".to_string(),
        None,
    ));

    // Write to temp file
    let temp_file = NamedTempFile::with_suffix(".cif").unwrap();
    write_mmcif_file(&structure, temp_file.path()).unwrap();

    // Read back and verify
    let content = std::fs::read_to_string(temp_file.path()).unwrap();
    assert!(content.contains("data_"));
    assert!(content.contains("_atom_site."));
    assert!(content.contains("GLY"));
    assert!(content.contains("10.500"));
}

#[test]
fn test_mmcif_roundtrip() {
    // Parse an existing mmCIF
    let mmcif_content = r#"
data_test
_entry.id ROUNDTRIP

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
"#;

    let structure1 = parse_mmcif_string(mmcif_content).unwrap();

    // Write to string
    let written_mmcif = write_mmcif_string(&structure1).unwrap();

    // Parse again
    let structure2 = parse_mmcif_string(&written_mmcif).unwrap();

    // Compare
    assert_eq!(structure1.atoms.len(), structure2.atoms.len());

    for (atom1, atom2) in structure1.atoms.iter().zip(structure2.atoms.iter()) {
        assert_eq!(atom1.name, atom2.name);
        assert_eq!(atom1.residue_name, atom2.residue_name);
        assert_eq!(atom1.chain_id, atom2.chain_id);
        assert_eq!(atom1.residue_seq, atom2.residue_seq);
        assert!((atom1.x - atom2.x).abs() < 0.001);
        assert!((atom1.y - atom2.y).abs() < 0.001);
        assert!((atom1.z - atom2.z).abs() < 0.001);
    }
}

#[test]
fn test_write_mmcif_with_multiple_chains() {
    let mut structure = PdbStructure::new();

    // Chain A
    structure.atoms.push(Atom::new(
        1,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        0.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    ));

    // Chain B
    structure.atoms.push(Atom::new(
        2,
        "CA".to_string(),
        None,
        "GLY".to_string(),
        "B".to_string(),
        1,
        10.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    ));

    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Should contain both chains
    assert!(mmcif_string.contains(" A "));
    assert!(mmcif_string.contains(" B "));
    assert!(mmcif_string.contains("ALA"));
    assert!(mmcif_string.contains("GLY"));
}

#[test]
fn test_write_mmcif_with_altlocs() {
    let mut structure = PdbStructure::new();

    structure.atoms.push(Atom::new(
        1,
        "CA".to_string(),
        Some('A'),
        "SER".to_string(),
        "A".to_string(),
        1,
        0.0,
        0.0,
        0.0,
        0.6,
        20.0,
        "C".to_string(),
        None,
    ));
    structure.atoms.push(Atom::new(
        2,
        "CA".to_string(),
        Some('B'),
        "SER".to_string(),
        "A".to_string(),
        1,
        0.1,
        0.1,
        0.1,
        0.4,
        20.0,
        "C".to_string(),
        None,
    ));

    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Should contain alternate location identifiers
    assert!(mmcif_string.contains(" A SER")); // altloc A
    assert!(mmcif_string.contains(" B SER")); // altloc B
}

#[test]
fn test_write_mmcif_with_hetatm() {
    let mut structure = PdbStructure::new();

    // Standard amino acid (should be ATOM)
    structure.atoms.push(Atom::new(
        1,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        0.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        None,
    ));

    // Water (should be HETATM)
    structure.atoms.push(Atom::new_hetatm(
        2,
        "O".to_string(),
        None,
        "HOH".to_string(),
        "W".to_string(),
        101,
        10.0,
        10.0,
        10.0,
        1.0,
        30.0,
        "O".to_string(),
        None,
    ));

    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Should have both ATOM and HETATM
    assert!(mmcif_string.contains("ATOM 1"));
    assert!(mmcif_string.contains("HETATM 2"));

    // Should have auth_seq_id header
    assert!(mmcif_string.contains("_atom_site.auth_seq_id"));

    // ATOM should have numeric label_seq_id (1) and auth_seq_id (1)
    // HETATM should have "." for label_seq_id and numeric auth_seq_id (101)
    assert!(mmcif_string.contains("ATOM 1 C CA . ALA A 1 1"));
    assert!(mmcif_string.contains("HETATM 2 O O . HOH W . 101"));
}

#[test]
fn test_write_mmcif_empty_structure() {
    let structure = PdbStructure::new();
    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Should still produce valid mmCIF header
    assert!(mmcif_string.contains("data_"));
    assert!(mmcif_string.contains("_entry.id"));
    assert!(mmcif_string.contains("loop_"));
    assert!(mmcif_string.contains("_atom_site."));
}

#[test]
fn test_write_mmcif_with_insertion_code() {
    let mut structure = PdbStructure::new();

    structure.atoms.push(Atom::new(
        1,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        52,
        0.0,
        0.0,
        0.0,
        1.0,
        20.0,
        "C".to_string(),
        Some('A'), // Insertion code
    ));

    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Should contain insertion code
    assert!(mmcif_string.contains(" A ")); // The insertion code in the output
}

#[test]
fn test_pdb_to_mmcif_conversion() {
    // Test conversion from PDB file to mmCIF
    let pdb_path = "examples/pdb_files/1UBQ.pdb";
    if std::path::Path::new(pdb_path).exists() {
        let structure = parse_pdb_file(pdb_path).unwrap();

        // Write to mmCIF
        let mmcif_string = write_mmcif_string(&structure).unwrap();

        // Parse the mmCIF back
        let structure2 = parse_mmcif_string(&mmcif_string).unwrap();

        // Verify atom count is preserved
        assert_eq!(structure.atoms.len(), structure2.atoms.len());

        // Verify first atom
        assert_eq!(structure.atoms[0].name, structure2.atoms[0].name);
        assert_eq!(
            structure.atoms[0].residue_name,
            structure2.atoms[0].residue_name
        );
    }
}

#[test]
fn test_write_mmcif_coordinate_precision() {
    let mut structure = PdbStructure::new();

    structure.atoms.push(Atom::new(
        1,
        "CA".to_string(),
        None,
        "ALA".to_string(),
        "A".to_string(),
        1,
        123.456789,
        -987.654321,
        0.001,
        0.567,
        89.123,
        "C".to_string(),
        None,
    ));

    let mmcif_string = write_mmcif_string(&structure).unwrap();

    // Verify coordinate precision (3 decimal places)
    assert!(mmcif_string.contains("123.457")); // x rounded
    assert!(mmcif_string.contains("-987.654")); // y rounded
    assert!(mmcif_string.contains("0.001")); // z

    // Parse back and verify precision is maintained
    let structure2 = parse_mmcif_string(&mmcif_string).unwrap();
    let atom = &structure2.atoms[0];

    assert!((atom.x - 123.457).abs() < 0.001);
    assert!((atom.y - (-987.654)).abs() < 0.001);
    assert!((atom.z - 0.001).abs() < 0.001);
}

/// When label and author numbering differ, disulfide partners must name the
/// residues as they appear in `structure.atoms`.
#[test]
fn test_mmcif_disulfide_partners_match_parsed_atoms() {
    let mmcif_content = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 S SG . CYS A 1 1 ? 0.000 0.000 0.000 1.00 10.00 22 CYS H SG 1
ATOM 2 S SG . CYS A 1 7 ? 2.030 0.000 0.000 1.00 10.00 95 CYS H SG 1
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.pdbx_dist_value
disulf1 disulf A CYS 1 A CYS 7 H CYS 22 H CYS 95 2.03
"#;
    let structure = parse_mmcif_string(mmcif_content).unwrap();
    assert_eq!(structure.ssbonds.len(), 1);
    let bond = &structure.ssbonds[0];
    for (chain, seq) in [
        (&bond.chain1_id, bond.residue1_seq),
        (&bond.chain2_id, bond.residue2_seq),
    ] {
        assert!(
            structure
                .atoms
                .iter()
                .any(|a| &a.chain_id == chain && a.residue_seq == seq && a.name == "SG"),
            "no SG atom for {chain}{seq}"
        );
    }
}
