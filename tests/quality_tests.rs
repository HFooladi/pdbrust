//! Integration tests for the quality module.
//!
//! These tests verify quality check functionality using real PDB files.

#![cfg(feature = "quality")]

use pdbrust::{parse_pdb_file, PdbStructure};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Basic Quality Check Tests
// ============================================================================

#[test]
fn test_has_ca_only_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // 1UBQ is a full-atom structure, not CA-only
    assert!(
        !structure.has_ca_only(),
        "1UBQ should not be CA-only"
    );
}

#[test]
fn test_has_multiple_models_single_model() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // 1UBQ is a single-model X-ray structure
    assert!(
        !structure.has_multiple_models(),
        "1UBQ should be single-model"
    );
    assert_eq!(structure.num_models(), 1);
}

#[test]
fn test_has_multiple_models_nmr() {
    let path = get_test_file("multi_model.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse multi_model.pdb");

    // multi_model.pdb has multiple models
    assert!(
        structure.has_multiple_models(),
        "multi_model.pdb should have multiple models"
    );
    assert!(structure.num_models() > 1);
}

#[test]
fn test_has_altlocs_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Check if structure has alternate locations
    let has_altlocs = structure.has_altlocs();
    let altloc_ids = structure.get_altloc_ids();

    // Consistency check
    assert_eq!(has_altlocs, !altloc_ids.is_empty());

    if has_altlocs {
        println!("Found altlocs: {:?}", altloc_ids);
    }
}

#[test]
fn test_has_hydrogens() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let has_h = structure.has_hydrogens();
    let h_count = structure.count_hydrogens();

    // Consistency check
    assert_eq!(has_h, h_count > 0);
}

#[test]
fn test_has_ssbonds() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Ubiquitin doesn't have disulfide bonds (no cysteines involved in SS)
    // Just check that the method works
    let _ = structure.has_ssbonds();
}

#[test]
fn test_has_conect() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Check connectivity records
    let _ = structure.has_conect();
}

// ============================================================================
// Quality Report Tests
// ============================================================================

#[test]
fn test_quality_report_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let report = structure.quality_report();

    // Basic sanity checks
    assert!(report.num_atoms > 0, "Should have atoms");
    assert!(report.num_chains > 0, "Should have chains");
    assert!(report.num_residues > 0, "Should have residues");
    assert_eq!(report.num_models, 1, "1UBQ should have 1 model");
    assert!(!report.has_ca_only, "1UBQ should not be CA-only");
}

#[test]
fn test_quality_report_consistency() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let report = structure.quality_report();

    // Verify consistency with individual method calls
    assert_eq!(report.has_ca_only, structure.has_ca_only());
    assert_eq!(report.has_multiple_models, structure.has_multiple_models());
    assert_eq!(report.has_altlocs, structure.has_altlocs());
    assert_eq!(report.num_chains, structure.get_num_chains());
    assert_eq!(report.num_atoms, structure.get_num_atoms());
    assert_eq!(report.num_residues, structure.get_num_residues());
    assert_eq!(report.has_hydrogens, structure.has_hydrogens());
    assert_eq!(report.has_ssbonds, structure.has_ssbonds());
    assert_eq!(report.has_conect, structure.has_conect());
}

#[test]
fn test_quality_report_is_analysis_ready() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let report = structure.quality_report();

    // 1UBQ should be analysis-ready (single model, no issues)
    if !report.has_altlocs {
        assert!(
            report.is_analysis_ready(),
            "1UBQ should be analysis-ready if no altlocs"
        );
    }
}

#[test]
fn test_quality_report_multi_model() {
    let path = get_test_file("multi_model.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse multi_model.pdb");

    let report = structure.quality_report();

    assert!(report.has_multiple_models);
    assert!(report.num_models > 1);

    // Multi-model structures are not analysis-ready by default
    assert!(
        !report.is_analysis_ready(),
        "Multi-model structure should not be analysis-ready"
    );
}

// ============================================================================
// Empty Structure Tests
// ============================================================================

#[test]
fn test_empty_structure_quality() {
    let structure = PdbStructure::new();

    assert!(structure.is_empty());
    assert!(!structure.has_ca_only());
    assert!(!structure.has_multiple_models());
    assert!(!structure.has_altlocs());
    assert!(!structure.has_hydrogens());
    assert!(!structure.has_ssbonds());
    assert!(!structure.has_conect());

    let report = structure.quality_report();
    assert_eq!(report.num_atoms, 0);
    assert_eq!(report.num_chains, 0);
    assert_eq!(report.num_residues, 0);
}

// ============================================================================
// Resolution Tests
// ============================================================================

#[test]
fn test_get_resolution() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Try to get resolution
    let resolution = structure.get_resolution();

    // 1UBQ should have resolution around 1.8 Å
    if let Some(res) = resolution {
        assert!(
            res > 0.5 && res < 5.0,
            "Resolution should be reasonable, got {}",
            res
        );
        println!("Resolution: {} Å", res);
    }
}

// ============================================================================
// Altloc Counting Tests
// ============================================================================

#[test]
fn test_count_altloc_atoms() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let altloc_count = structure.count_altloc_atoms();
    let has_altlocs = structure.has_altlocs();

    // Consistency check
    assert_eq!(has_altlocs, altloc_count > 0);
}

// ============================================================================
// Compare Multiple Structures
// ============================================================================

#[test]
fn test_compare_structure_quality() {
    let ubq = parse_pdb_file(get_test_file("1UBQ.pdb")).expect("Failed to parse 1UBQ.pdb");
    let hm2 = parse_pdb_file(get_test_file("8HM2.pdb")).expect("Failed to parse 8HM2.pdb");

    let ubq_report = ubq.quality_report();
    let hm2_report = hm2.quality_report();

    // Both should have valid reports
    assert!(ubq_report.num_atoms > 0);
    assert!(hm2_report.num_atoms > 0);

    // Different structures should have different atom counts
    // (unless they happen to be the same size)
    println!("1UBQ atoms: {}", ubq_report.num_atoms);
    println!("8HM2 atoms: {}", hm2_report.num_atoms);
}

// ============================================================================
// HETATM Detection Tests
// ============================================================================

#[test]
fn test_has_hetatm() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Just check that the method works
    let _ = structure.has_hetatm();
}

// ============================================================================
// Default Report Tests
// ============================================================================

#[test]
fn test_quality_report_default() {
    let report = pdbrust::quality::QualityReport::default();

    assert!(!report.has_ca_only);
    assert!(!report.has_multiple_models);
    assert!(!report.has_altlocs);
    assert_eq!(report.num_chains, 0);
    assert_eq!(report.num_models, 0);
    assert_eq!(report.num_atoms, 0);
}

#[test]
fn test_quality_report_is_clean() {
    let mut report = pdbrust::quality::QualityReport::default();
    report.num_atoms = 100;
    report.has_ca_only = false;
    report.has_altlocs = false;

    assert!(report.is_clean());

    // CA-only is not clean
    report.has_ca_only = true;
    assert!(!report.is_clean());

    // Altlocs are not clean
    report.has_ca_only = false;
    report.has_altlocs = true;
    assert!(!report.is_clean());
}
