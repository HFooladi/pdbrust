use pdbrust::parse_pdb_file;
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

#[test]
fn test_parse_multi_model() {
    let content = "\
MODEL        1
ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
ENDMDL
MODEL        2
ATOM      3  N   ALA A   1      27.147  14.199   3.215  1.00 13.79           N
ATOM      4  CA  ALA A   1      26.147  13.199   2.215  1.00 12.79           C
ENDMDL
END
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.models.len(), 2);

    let model1 = &structure.models[0];
    assert_eq!(model1.serial, 1);
    assert_eq!(model1.atoms.len(), 2);
    assert_eq!(model1.atoms[0].serial, 1);
    assert_eq!(model1.atoms[1].serial, 2);

    let model2 = &structure.models[1];
    assert_eq!(model2.serial, 2);
    assert_eq!(model2.atoms.len(), 2);
    assert_eq!(model2.atoms[0].serial, 3);
    assert_eq!(model2.atoms[1].serial, 4);

    // All atoms should be in the structure.atoms as well
    assert_eq!(structure.atoms.len(), 4);

    // Check the coordinates to ensure the atoms are correctly assigned
    assert!((structure.atoms[0].x - 27.047).abs() < 1e-6);
    assert!((structure.atoms[2].x - 27.147).abs() < 1e-6);
}

// ============================================================================
// Edge Case Tests - HETATM and Ligands
// ============================================================================

#[test]
fn test_parse_hetatm() {
    let content =
        "HETATM    1  C1  LIG A   1      10.000  20.000  30.000  1.00 25.00           C\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    let atom = &structure.atoms[0];
    assert_eq!(atom.serial, 1);
    assert_eq!(atom.name, "C1");
    assert_eq!(atom.residue_name, "LIG");
    assert_eq!(atom.element, "C");
    assert!(atom.is_hetatm); // Verify is_hetatm is true for HETATM records
}

#[test]
fn test_parse_water_molecules() {
    let content = "\
HETATM  100  O   HOH A 201      15.000  25.000  35.000  1.00 30.00           O
HETATM  101  O   HOH A 202      16.000  26.000  36.000  1.00 31.00           O
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 2);
    assert!(structure.atoms.iter().all(|a| a.residue_name == "HOH"));
    assert!(structure.atoms.iter().all(|a| a.element == "O"));
    assert!(structure.atoms.iter().all(|a| a.is_hetatm)); // All water molecules are HETATM
}

#[test]
fn test_parse_mixed_atom_hetatm() {
    let content = "\
ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
HETATM    3  C1  LIG A 100      10.000  20.000  30.000  1.00 25.00           C
HETATM    4  O   HOH A 201      15.000  25.000  35.000  1.00 30.00           O
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 4);
    // Standard ATOM records
    assert_eq!(structure.atoms[0].residue_name, "ALA");
    assert!(!structure.atoms[0].is_hetatm);
    assert_eq!(structure.atoms[1].residue_name, "ALA");
    assert!(!structure.atoms[1].is_hetatm);
    // HETATM records (ligand and water)
    assert_eq!(structure.atoms[2].residue_name, "LIG");
    assert!(structure.atoms[2].is_hetatm);
    assert_eq!(structure.atoms[3].residue_name, "HOH");
    assert!(structure.atoms[3].is_hetatm);
}

// ============================================================================
// Edge Case Tests - Alternative Locations (Altloc)
// ============================================================================

#[test]
fn test_parse_altloc() {
    let content = "\
ATOM      1  N  AALA A   1      27.047  14.099   3.115  0.60 13.79           N
ATOM      2  N  BALA A   1      27.147  14.199   3.215  0.40 14.79           N
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 2);
    assert_eq!(structure.atoms[0].alt_loc, Some('A'));
    assert_eq!(structure.atoms[1].alt_loc, Some('B'));
    assert!((structure.atoms[0].occupancy - 0.60).abs() < 1e-6);
    assert!((structure.atoms[1].occupancy - 0.40).abs() < 1e-6);
}

// ============================================================================
// Edge Case Tests - Insertion Codes
// ============================================================================

#[test]
fn test_parse_insertion_code() {
    let content = "\
ATOM      1  N   ALA A  10      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  N   ALA A  10A     27.147  14.199   3.215  1.00 14.79           N
ATOM      3  N   ALA A  10B     27.247  14.299   3.315  1.00 15.79           N
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 3);
    assert_eq!(structure.atoms[0].ins_code, None);
    assert_eq!(structure.atoms[1].ins_code, Some('A'));
    assert_eq!(structure.atoms[2].ins_code, Some('B'));
    // All have same residue_seq
    assert!(structure.atoms.iter().all(|a| a.residue_seq == 10));
}

// ============================================================================
// Edge Case Tests - Negative Residue Numbers
// ============================================================================

#[test]
fn test_parse_negative_residue_number() {
    let content = "\
ATOM      1  N   ALA A  -5      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  ALA A  -5      26.047  13.099   2.115  1.00 12.79           C
ATOM      3  N   ALA A  -4      28.047  15.099   4.115  1.00 13.79           N
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 3);
    assert_eq!(structure.atoms[0].residue_seq, -5);
    assert_eq!(structure.atoms[1].residue_seq, -5);
    assert_eq!(structure.atoms[2].residue_seq, -4);
}

// ============================================================================
// Edge Case Tests - Multiple Chains
// ============================================================================

#[test]
fn test_parse_multiple_chains() {
    let content = "\
ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
ATOM      3  N   GLY B   1      30.047  17.099   6.115  1.00 15.79           N
ATOM      4  CA  GLY B   1      29.047  16.099   5.115  1.00 14.79           C
ATOM      5  N   VAL C   1      33.047  20.099   9.115  1.00 17.79           N
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 5);

    let chain_ids = structure.get_chain_ids();
    assert_eq!(chain_ids.len(), 3);
    assert!(chain_ids.contains(&"A".to_string()));
    assert!(chain_ids.contains(&"B".to_string()));
    assert!(chain_ids.contains(&"C".to_string()));
}

// ============================================================================
// Edge Case Tests - Coordinates
// ============================================================================

#[test]
fn test_parse_negative_coordinates() {
    let content =
        "ATOM      1  N   ALA A   1     -27.047 -14.099  -3.115  1.00 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    let atom = &structure.atoms[0];
    assert!((atom.x - (-27.047)).abs() < 1e-6);
    assert!((atom.y - (-14.099)).abs() < 1e-6);
    assert!((atom.z - (-3.115)).abs() < 1e-6);
}

#[test]
fn test_parse_large_coordinates() {
    let content =
        "ATOM      1  N   ALA A   1     999.999 999.999 999.999  1.00 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    let atom = &structure.atoms[0];
    assert!((atom.x - 999.999).abs() < 1e-6);
    assert!((atom.y - 999.999).abs() < 1e-6);
    assert!((atom.z - 999.999).abs() < 1e-6);
}

#[test]
fn test_parse_zero_coordinates() {
    let content =
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
    let atom = &structure.atoms[0];
    assert!((atom.x).abs() < 1e-6);
    assert!((atom.y).abs() < 1e-6);
    assert!((atom.z).abs() < 1e-6);
}

// ============================================================================
// Edge Case Tests - Various Atom Names
// ============================================================================

#[test]
fn test_parse_four_char_atom_name() {
    // Some atom names are 4 characters (e.g., calcium CA vs alpha carbon CA)
    let content = "\
ATOM      1  CA  ALA A   1      27.047  14.099   3.115  1.00 13.79           C
ATOM      2 CA   ALA A   2      28.047  15.099   4.115  1.00 14.79          CA
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 2);
    assert_eq!(structure.atoms[0].name, "CA");
    assert_eq!(structure.atoms[1].name, "CA");
}

#[test]
fn test_parse_hydrogen_atoms() {
    let content = "\
ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  H   ALA A   1      27.547  14.599   3.615  1.00 15.79           H
ATOM      3 1H   ALA A   1      27.147  14.199   3.215  1.00 14.79           H
ATOM      4 2H   ALA A   1      27.247  14.299   3.315  1.00 14.79           H
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 4);
    assert_eq!(structure.atoms[0].element, "N");
    assert_eq!(structure.atoms[1].element, "H");
    assert_eq!(structure.atoms[2].element, "H");
    assert_eq!(structure.atoms[3].element, "H");
}

// ============================================================================
// Edge Case Tests - TER Records
// ============================================================================

#[test]
fn test_parse_with_ter_record() {
    let content = "\
ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
TER       3      ALA A   1
ATOM      4  N   GLY B   1      30.047  17.099   6.115  1.00 15.79           N
ATOM      5  CA  GLY B   1      29.047  16.099   5.115  1.00 14.79           C
TER       6      GLY B   1
END
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    // TER records should not add atoms
    assert_eq!(structure.atoms.len(), 4);
    let chain_ids = structure.get_chain_ids();
    assert_eq!(chain_ids.len(), 2);
}

// ============================================================================
// Edge Case Tests - Various REMARK Types
// ============================================================================

#[test]
fn test_parse_multiple_remarks() {
    let content = "\
REMARK   1 REFERENCE 1
REMARK   2 RESOLUTION.    2.00 ANGSTROMS.
REMARK   3 REFINEMENT.
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B
REMARK 465 MISSING RESIDUES
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.remarks.len(), 5);
    assert_eq!(structure.remarks[0].number, 1);
    assert_eq!(structure.remarks[1].number, 2);
    assert_eq!(structure.remarks[2].number, 3);
    assert_eq!(structure.remarks[3].number, 350);
    assert_eq!(structure.remarks[4].number, 465);
}

// ============================================================================
// Edge Case Tests - CONECT with Partial Connections
// ============================================================================

#[test]
fn test_parse_conect_partial() {
    let content = "\
CONECT    1    2
CONECT    3    4    5
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.connects.len(), 2);

    let conect1 = &structure.connects[0];
    assert_eq!(conect1.atom1, 1);
    assert_eq!(conect1.atom2, 2);
    assert_eq!(conect1.atom3, None);
    assert_eq!(conect1.atom4, None);

    let conect2 = &structure.connects[1];
    assert_eq!(conect2.atom1, 3);
    assert_eq!(conect2.atom2, 4);
    assert_eq!(conect2.atom3, Some(5));
    assert_eq!(conect2.atom4, None);
}

// ============================================================================
// Edge Case Tests - Large Serial Numbers
// ============================================================================

#[test]
fn test_parse_large_serial_numbers() {
    let content = "\
ATOM  99998  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM  99999  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 2);
    assert_eq!(structure.atoms[0].serial, 99998);
    assert_eq!(structure.atoms[1].serial, 99999);
}

// ============================================================================
// Edge Case Tests - Special Residue Names
// ============================================================================

#[test]
fn test_parse_non_standard_residues() {
    let content = "\
ATOM      1  N   MSE A   1      27.047  14.099   3.115  1.00 13.79           N
ATOM      2  CA  MSE A   1      26.047  13.099   2.115  1.00 12.79           C
ATOM      3  SE  MSE A   1      25.047  12.099   1.115  1.00 11.79          SE
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 3);
    assert!(structure.atoms.iter().all(|a| a.residue_name == "MSE"));
    assert_eq!(structure.atoms[2].element, "SE");
}

#[test]
fn test_parse_nucleic_acid_residues() {
    let content = "\
ATOM      1  P     A A   1      27.047  14.099   3.115  1.00 13.79           P
ATOM      2  O5'   A A   1      26.047  13.099   2.115  1.00 12.79           O
ATOM      3  C5'   A A   1      25.047  12.099   1.115  1.00 11.79           C
ATOM      4  P     G A   2      30.047  17.099   6.115  1.00 15.79           P
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 4);
    assert_eq!(structure.atoms[0].residue_name, "A");
    assert_eq!(structure.atoms[3].residue_name, "G");
}

// ============================================================================
// Edge Case Tests - Whitespace Handling
// ============================================================================

#[test]
fn test_parse_lines_with_trailing_whitespace() {
    let content =
        "ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N   \n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 1);
}

#[test]
fn test_parse_file_with_blank_lines() {
    let content = "\
HEADER    TEST

ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 13.79           N

ATOM      2  CA  ALA A   1      26.047  13.099   2.115  1.00 12.79           C

END
";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert_eq!(structure.atoms.len(), 2);
}

// ============================================================================
// Edge Case Tests - B-factor Edge Cases
// ============================================================================

#[test]
fn test_parse_zero_bfactor() {
    let content =
        "ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00  0.00           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert!((structure.atoms[0].temp_factor).abs() < 1e-6);
}

#[test]
fn test_parse_high_bfactor() {
    let content =
        "ATOM      1  N   ALA A   1      27.047  14.099   3.115  1.00 99.99           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert!((structure.atoms[0].temp_factor - 99.99).abs() < 1e-6);
}

// ============================================================================
// Edge Case Tests - Occupancy Edge Cases
// ============================================================================

#[test]
fn test_parse_partial_occupancy() {
    let content =
        "ATOM      1  N   ALA A   1      27.047  14.099   3.115  0.50 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert!((structure.atoms[0].occupancy - 0.50).abs() < 1e-6);
}

#[test]
fn test_parse_zero_occupancy() {
    let content =
        "ATOM      1  N   ALA A   1      27.047  14.099   3.115  0.00 13.79           N\n";
    let file = create_test_pdb(content);
    let result = parse_pdb_file(file.path());
    assert!(result.is_ok());
    let structure = result.unwrap();
    assert!((structure.atoms[0].occupancy).abs() < 1e-6);
}
