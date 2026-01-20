//! Integration tests for the DSSP secondary structure assignment module.
//!
//! These tests verify secondary structure assignment functionality using real PDB files.

#![cfg(feature = "dssp")]

use pdbrust::dssp::{SecondaryStructure, SecondaryStructureAssignment};
use pdbrust::{PdbStructure, parse_mmcif_file, parse_pdb_file};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Basic Assignment Tests
// ============================================================================

#[test]
fn test_dssp_1ubq_pdb() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss = structure.assign_secondary_structure();

    // 1UBQ is ubiquitin, a small protein with both helix and sheet content
    assert!(!ss.is_empty(), "Assignment should not be empty");

    // 1UBQ has known secondary structure content
    // Alpha helix ~15-20%, Beta sheet ~40-50%
    println!("1UBQ SS assignment:");
    println!("  Total residues: {}", ss.len());
    println!("  Helix fraction: {:.1}%", ss.helix_fraction * 100.0);
    println!("  Sheet fraction: {:.1}%", ss.sheet_fraction * 100.0);
    println!("  Coil fraction:  {:.1}%", ss.coil_fraction * 100.0);
    println!("  SS string: {}", ss.as_codes());

    // Basic sanity checks
    assert!(ss.helix_fraction >= 0.0 && ss.helix_fraction <= 1.0);
    assert!(ss.sheet_fraction >= 0.0 && ss.sheet_fraction <= 1.0);
    assert!(ss.coil_fraction >= 0.0 && ss.coil_fraction <= 1.0);

    // Fractions should sum to 1
    let total = ss.helix_fraction + ss.sheet_fraction + ss.coil_fraction;
    assert!(
        (total - 1.0).abs() < 0.01,
        "Fractions should sum to 1, got {}",
        total
    );
}

#[test]
fn test_dssp_1crn_mmcif() {
    let path = get_test_file("1CRN.cif");
    let structure = parse_mmcif_file(&path).expect("Failed to parse 1CRN.cif");

    let ss = structure.assign_secondary_structure();

    // Crambin (1CRN) is a small plant protein with helix and sheet
    assert!(!ss.is_empty(), "Assignment should not be empty");

    let ss_string = ss.as_codes();
    println!("1CRN SS string: {}", ss_string);

    // All characters should be valid DSSP codes (H, G, I, P, E, B, T, S, C)
    for c in ss_string.chars() {
        assert!("HGIPEBTSC".contains(c), "Invalid SS code: {}", c);
    }
}

#[test]
fn test_ss_string_length_matches_residues() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss = structure.assign_secondary_structure();
    let ss_string = structure.secondary_structure_string();

    assert_eq!(
        ss_string.len(),
        ss.len(),
        "SS string length should match number of assignments"
    );
}

// ============================================================================
// Secondary Structure Type Tests
// ============================================================================

#[test]
fn test_secondary_structure_codes() {
    // Test that all codes are correct
    assert_eq!(SecondaryStructure::AlphaHelix.code(), 'H');
    assert_eq!(SecondaryStructure::Helix310.code(), 'G');
    assert_eq!(SecondaryStructure::PiHelix.code(), 'I');
    assert_eq!(SecondaryStructure::KappaHelix.code(), 'P');
    assert_eq!(SecondaryStructure::ExtendedStrand.code(), 'E');
    assert_eq!(SecondaryStructure::BetaBridge.code(), 'B');
    assert_eq!(SecondaryStructure::Turn.code(), 'T');
    assert_eq!(SecondaryStructure::Bend.code(), 'S');
    assert_eq!(SecondaryStructure::Coil.code(), 'C');
}

#[test]
fn test_secondary_structure_from_code() {
    assert_eq!(
        SecondaryStructure::from_code('H'),
        Some(SecondaryStructure::AlphaHelix)
    );
    assert_eq!(
        SecondaryStructure::from_code('E'),
        Some(SecondaryStructure::ExtendedStrand)
    );
    assert_eq!(
        SecondaryStructure::from_code('P'),
        Some(SecondaryStructure::KappaHelix)
    );
    assert_eq!(SecondaryStructure::from_code('X'), None);
    assert_eq!(SecondaryStructure::from_code('h'), None); // Case sensitive
}

#[test]
fn test_is_helix() {
    assert!(SecondaryStructure::AlphaHelix.is_helix());
    assert!(SecondaryStructure::Helix310.is_helix());
    assert!(SecondaryStructure::PiHelix.is_helix());
    assert!(SecondaryStructure::KappaHelix.is_helix());
    assert!(!SecondaryStructure::ExtendedStrand.is_helix());
    assert!(!SecondaryStructure::Coil.is_helix());
}

#[test]
fn test_is_sheet() {
    assert!(SecondaryStructure::ExtendedStrand.is_sheet());
    assert!(SecondaryStructure::BetaBridge.is_sheet());
    assert!(!SecondaryStructure::AlphaHelix.is_sheet());
    assert!(!SecondaryStructure::Coil.is_sheet());
}

#[test]
fn test_is_coil() {
    assert!(SecondaryStructure::Coil.is_coil());
    assert!(SecondaryStructure::Turn.is_coil());
    assert!(SecondaryStructure::Bend.is_coil());
    assert!(!SecondaryStructure::AlphaHelix.is_coil());
    assert!(!SecondaryStructure::ExtendedStrand.is_coil());
}

// ============================================================================
// Composition Tests
// ============================================================================

#[test]
fn test_secondary_structure_composition() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let (helix, sheet, coil) = structure.secondary_structure_composition();

    // All fractions should be between 0 and 1
    assert!((0.0..=1.0).contains(&helix));
    assert!((0.0..=1.0).contains(&sheet));
    assert!((0.0..=1.0).contains(&coil));

    // Fractions should sum to 1
    assert!((helix + sheet + coil - 1.0).abs() < 0.01);
}

// ============================================================================
// Empty Structure Tests
// ============================================================================

#[test]
fn test_empty_structure() {
    let structure = PdbStructure::new();
    let ss = structure.assign_secondary_structure();

    assert!(ss.is_empty());
    assert_eq!(ss.len(), 0);
    assert!(ss.warnings.iter().any(|w| w.contains("No atoms")));
}

#[test]
fn test_empty_structure_composition() {
    let structure = PdbStructure::new();
    let (helix, sheet, coil) = structure.secondary_structure_composition();

    // Empty structure should have zero fractions
    assert_eq!(helix, 0.0);
    assert_eq!(sheet, 0.0);
    assert_eq!(coil, 0.0);
}

// ============================================================================
// Edge Case Tests
// ============================================================================

#[test]
fn test_single_residue_structure() {
    use pdbrust::records::Atom;

    let mut structure = PdbStructure::new();
    structure.atoms = vec![
        Atom {
            serial: 1,
            name: "N".to_string(),
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
            element: "N".to_string(),
        },
        Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            ins_code: None,
            x: 1.5,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        },
        Atom {
            serial: 3,
            name: "C".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            ins_code: None,
            x: 2.5,
            y: 1.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        },
        Atom {
            serial: 4,
            name: "O".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            ins_code: None,
            x: 2.5,
            y: 2.2,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "O".to_string(),
        },
    ];

    let ss = structure.assign_secondary_structure();

    // Should have one residue assignment
    assert_eq!(ss.len(), 1);
    // Single residue should be coil (no H-bonds possible)
    assert_eq!(ss.residue_assignments[0].ss, SecondaryStructure::Coil);
}

// ============================================================================
// CA-Only Structure Tests
// ============================================================================

#[test]
fn test_ca_only_structure() {
    use pdbrust::records::Atom;

    let mut structure = PdbStructure::new();
    structure.atoms = vec![
        Atom {
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
        },
        Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "GLY".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 2,
            ins_code: None,
            x: 3.8,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        },
        Atom {
            serial: 3,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "VAL".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 3,
            ins_code: None,
            x: 7.6,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        },
    ];

    let ss = structure.assign_secondary_structure();

    // Should have 3 residue assignments
    assert_eq!(ss.len(), 3);

    // CA-only structure should warn and return all coil
    assert!(ss.warnings.iter().any(|w| w.contains("CA-only")));

    for assignment in &ss.residue_assignments {
        assert_eq!(assignment.ss, SecondaryStructure::Coil);
    }
}

// ============================================================================
// Assignment Consistency Tests
// ============================================================================

#[test]
fn test_assignment_consistency() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Run assignment twice - should get same results
    let ss1 = structure.assign_secondary_structure();
    let ss2 = structure.assign_secondary_structure();

    assert_eq!(ss1.len(), ss2.len());
    assert_eq!(ss1.to_string(), ss2.to_string());
    assert_eq!(ss1.helix_fraction, ss2.helix_fraction);
    assert_eq!(ss1.sheet_fraction, ss2.sheet_fraction);
}

// ============================================================================
// Multiple Chain Tests
// ============================================================================

#[test]
fn test_multiple_chains() {
    let path = get_test_file("8HM2.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 8HM2.pdb");

    let ss = structure.assign_secondary_structure();

    // Should have assignments
    assert!(!ss.is_empty());

    // Check that multiple chains are handled
    let chain_ids: std::collections::HashSet<_> =
        ss.residue_assignments.iter().map(|a| &a.chain_id).collect();

    println!("Chains in 8HM2: {:?}", chain_ids);
}

// ============================================================================
// Statistics Tests
// ============================================================================

#[test]
fn test_statistics_computation() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss = structure.assign_secondary_structure();

    // Verify counts match fractions
    let total = ss.len() as f64;
    if total > 0.0 {
        let expected_helix_frac = ss.helix_count as f64 / total;
        let expected_sheet_frac = ss.sheet_count as f64 / total;
        let expected_coil_frac = ss.coil_count as f64 / total;

        assert!((ss.helix_fraction - expected_helix_frac).abs() < 0.001);
        assert!((ss.sheet_fraction - expected_sheet_frac).abs() < 0.001);
        assert!((ss.coil_fraction - expected_coil_frac).abs() < 0.001);
    }

    // Counts should sum to total
    assert_eq!(ss.helix_count + ss.sheet_count + ss.coil_count, ss.len());
}

// ============================================================================
// Per-Residue Assignment Tests
// ============================================================================

#[test]
fn test_per_residue_assignments() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss = structure.assign_secondary_structure();

    // Check that each residue has valid data
    for assignment in &ss.residue_assignments {
        assert!(!assignment.chain_id.is_empty());
        assert!(!assignment.residue_name.is_empty());
        // SS should be one of the valid types
        assert!(
            "HGIPEBTSC".contains(assignment.ss.code()),
            "Invalid SS: {}",
            assignment.ss.code()
        );
    }
}

// ============================================================================
// PPII Detection Tests
// ============================================================================

#[test]
fn test_ppii_detection() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss = structure.assign_secondary_structure();
    let ss_string = ss.as_codes();

    // PPII regions may or may not be present
    // Just verify the algorithm handles it without errors
    let ppii_count = ss_string.chars().filter(|&c| c == 'P').count();

    println!("PPII residues in 1UBQ: {}", ppii_count);

    // P codes should only appear in regions that would otherwise be coil
    // This is just a sanity check that the algorithm runs
}

// ============================================================================
// Display/Format Tests
// ============================================================================

#[test]
fn test_display_trait() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss = structure.assign_secondary_structure();

    // Test Display trait
    let display_string = format!("{}", ss);
    let to_string_result = ss.as_codes();

    assert_eq!(display_string, to_string_result);
}

#[test]
fn test_secondary_structure_display() {
    // Test Display trait for SecondaryStructure
    assert_eq!(format!("{}", SecondaryStructure::AlphaHelix), "H");
    assert_eq!(format!("{}", SecondaryStructure::ExtendedStrand), "E");
    assert_eq!(format!("{}", SecondaryStructure::Coil), "C");
    assert_eq!(format!("{}", SecondaryStructure::KappaHelix), "P");
}

// ============================================================================
// Default Trait Tests
// ============================================================================

#[test]
fn test_default_secondary_structure() {
    let ss = SecondaryStructure::default();
    assert_eq!(ss, SecondaryStructure::Coil);
}

#[test]
fn test_default_assignment() {
    let assignment = SecondaryStructureAssignment::default();
    assert!(assignment.is_empty());
    assert_eq!(assignment.helix_fraction, 0.0);
    assert_eq!(assignment.sheet_fraction, 0.0);
    assert_eq!(assignment.coil_fraction, 0.0);
}
