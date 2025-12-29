//! Integration tests for the summary module.
//!
//! These tests verify the unified summary functionality using real PDB files.

#![cfg(feature = "summary")]

use pdbrust::{parse_pdb_file, PdbStructure};
use pdbrust::summary::{StructureSummary, batch_summarize, summaries_to_csv};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Basic Summary Tests
// ============================================================================

#[test]
fn test_summary_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();

    // Basic sanity checks
    assert!(summary.num_atoms > 0, "Should have atoms");
    assert!(summary.num_residues > 0, "Should have residues");
    assert!(summary.num_chains > 0, "Should have chains");
    assert_eq!(summary.num_models, 1, "1UBQ should have 1 model");
}

#[test]
fn test_summary_quality_fields() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();

    // 1UBQ should have reasonable quality
    assert!(!summary.has_ca_only, "1UBQ is full-atom");
    assert!(!summary.has_multiple_models, "1UBQ is X-ray");
}

#[test]
fn test_summary_descriptor_fields() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();

    // 1UBQ should have valid descriptors
    assert!(summary.radius_of_gyration > 0.0, "Rg should be positive");
    assert!(summary.max_ca_distance > 0.0, "Max distance should be positive");
    assert!(summary.compactness_index > 0.0, "Compactness should be positive");

    // Composition should be valid
    assert!(summary.hydrophobic_ratio >= 0.0 && summary.hydrophobic_ratio <= 1.0);
    assert!(summary.glycine_ratio >= 0.0 && summary.glycine_ratio <= 1.0);
}

#[test]
fn test_summary_consistency_with_parts() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();
    let quality = structure.quality_report();
    let descriptors = structure.structure_descriptors();

    // Quality fields should match
    assert_eq!(summary.has_ca_only, quality.has_ca_only);
    assert_eq!(summary.has_multiple_models, quality.has_multiple_models);
    assert_eq!(summary.has_altlocs, quality.has_altlocs);
    assert_eq!(summary.num_chains, quality.num_chains);
    assert_eq!(summary.num_models, quality.num_models);

    // Descriptor fields should match
    assert_eq!(summary.num_residues, descriptors.num_residues);
    assert!((summary.radius_of_gyration - descriptors.radius_of_gyration).abs() < 1e-10);
    assert!((summary.hydrophobic_ratio - descriptors.hydrophobic_ratio).abs() < 1e-10);
}

// ============================================================================
// Analysis Readiness Tests
// ============================================================================

#[test]
fn test_is_analysis_ready() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();

    // 1UBQ should be analysis-ready if no altlocs
    if !summary.has_altlocs {
        assert!(summary.is_analysis_ready(), "1UBQ should be analysis-ready");
    }
}

#[test]
fn test_is_analysis_ready_multi_model() {
    let path = get_test_file("multi_model.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse multi_model.pdb");

    let summary = structure.summary();

    // Multi-model structure should not be analysis-ready
    assert!(summary.has_multiple_models);
    assert!(!summary.is_analysis_ready());
}

#[test]
fn test_is_clean() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();

    // Check is_clean
    if !summary.has_ca_only && !summary.has_altlocs {
        assert!(summary.is_clean());
    }
}

// ============================================================================
// CSV Output Tests
// ============================================================================

#[test]
fn test_field_names() {
    let names = StructureSummary::field_names();

    assert!(!names.is_empty());
    assert!(names.contains(&"num_atoms"));
    assert!(names.contains(&"radius_of_gyration"));
    assert!(names.contains(&"has_altlocs"));
    assert!(names.contains(&"hydrophobic_ratio"));
}

#[test]
fn test_to_csv_values() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();
    let values = summary.to_csv_values();

    // Should have same number of values as field names
    assert_eq!(
        values.len(),
        StructureSummary::field_names().len(),
        "CSV values should match field names count"
    );

    // All values should be non-empty strings
    for (i, value) in values.iter().enumerate() {
        assert!(
            !value.is_empty(),
            "Field {} should have a value",
            StructureSummary::field_names()[i]
        );
    }
}

#[test]
fn test_summaries_to_csv_with_header() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summaries = vec![structure.summary()];
    let csv = summaries_to_csv(&summaries, true);

    // Should have header + data row
    let lines: Vec<&str> = csv.lines().collect();
    assert_eq!(lines.len(), 2, "Should have header and one data row");

    // Header should contain field names
    let header = lines[0];
    assert!(header.contains("num_atoms"));
    assert!(header.contains("radius_of_gyration"));
}

#[test]
fn test_summaries_to_csv_without_header() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summaries = vec![structure.summary()];
    let csv = summaries_to_csv(&summaries, false);

    // Should have only data row
    let lines: Vec<&str> = csv.lines().collect();
    assert_eq!(lines.len(), 1, "Should have only one data row");

    // Should not contain field names at the start
    assert!(!lines[0].starts_with("has_ca_only"));
}

// ============================================================================
// Batch Processing Tests
// ============================================================================

#[test]
fn test_batch_summarize() {
    let ubq = parse_pdb_file(get_test_file("1UBQ.pdb")).expect("Failed to parse 1UBQ.pdb");
    let hm2 = parse_pdb_file(get_test_file("8HM2.pdb")).expect("Failed to parse 8HM2.pdb");

    let structures = vec![ubq, hm2];
    let summaries = batch_summarize(&structures);

    assert_eq!(summaries.len(), 2);

    // Both should have valid summaries
    for summary in &summaries {
        assert!(summary.num_atoms > 0);
        assert!(summary.num_residues > 0);
    }
}

#[test]
fn test_batch_summarize_empty() {
    let structures: Vec<PdbStructure> = vec![];
    let summaries = batch_summarize(&structures);

    assert!(summaries.is_empty());
}

#[test]
fn test_batch_to_csv() {
    let ubq = parse_pdb_file(get_test_file("1UBQ.pdb")).expect("Failed to parse 1UBQ.pdb");
    let hm2 = parse_pdb_file(get_test_file("8HM2.pdb")).expect("Failed to parse 8HM2.pdb");

    let structures = vec![ubq, hm2];
    let summaries = batch_summarize(&structures);
    let csv = summaries_to_csv(&summaries, true);

    // Should have header + 2 data rows
    let lines: Vec<&str> = csv.lines().collect();
    assert_eq!(lines.len(), 3, "Should have header and two data rows");
}

// ============================================================================
// Edge Cases
// ============================================================================

#[test]
fn test_empty_structure_summary() {
    let structure = PdbStructure::new();
    let summary = structure.summary();

    assert_eq!(summary.num_atoms, 0);
    assert_eq!(summary.num_residues, 0);
    assert_eq!(summary.num_chains, 0);
    assert!(!summary.is_analysis_ready());
    assert!(!summary.is_clean());
}

#[test]
fn test_summary_default() {
    let summary = StructureSummary::default();

    assert_eq!(summary.num_atoms, 0);
    assert!(!summary.has_ca_only);
    assert!(summary.aa_composition.is_empty());
    assert_eq!(summary.radius_of_gyration, 0.0);
}

// ============================================================================
// from_parts Tests
// ============================================================================

#[test]
fn test_from_parts() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let quality = structure.quality_report();
    let descriptors = structure.structure_descriptors();

    let summary_from_parts = StructureSummary::from_parts(quality.clone(), descriptors.clone());
    let summary_direct = structure.summary();

    // Both methods should produce equivalent results
    assert_eq!(summary_from_parts.num_atoms, summary_direct.num_atoms);
    assert_eq!(summary_from_parts.has_ca_only, summary_direct.has_ca_only);
    assert!((summary_from_parts.radius_of_gyration - summary_direct.radius_of_gyration).abs() < 1e-10);
}

// ============================================================================
// Composition Access Tests
// ============================================================================

#[test]
fn test_aa_composition_in_summary() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let summary = structure.summary();

    // Should have AA composition
    assert!(!summary.aa_composition.is_empty());

    // All fractions should sum to approximately 1.0
    let total: f64 = summary.aa_composition.values().sum();
    assert!(
        (total - 1.0).abs() < 1e-10,
        "AA composition should sum to 1.0, got {}",
        total
    );
}
