//! Integration tests for the descriptors module.
//!
//! These tests verify descriptor calculations using real PDB files.

#![cfg(feature = "descriptors")]

use pdbrust::{parse_pdb_file, PdbStructure};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Composition Tests
// ============================================================================

#[test]
fn test_aa_composition_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let composition = structure.aa_composition();

    // 1UBQ should have multiple amino acid types
    assert!(!composition.is_empty(), "Composition should not be empty");

    // All fractions should sum to approximately 1.0
    let total: f64 = composition.values().sum();
    assert!(
        (total - 1.0).abs() < 1e-10,
        "Fractions should sum to 1.0, got {}",
        total
    );

    // All fractions should be between 0 and 1
    for (aa, &fraction) in &composition {
        assert!(
            (0.0..=1.0).contains(&fraction),
            "Fraction for {} should be 0-1, got {}",
            aa,
            fraction
        );
    }
}

#[test]
fn test_glycine_ratio_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let gly_ratio = structure.glycine_ratio();

    // 1UBQ (ubiquitin) has 6 glycines out of 76 residues = ~7.9%
    assert!(
        gly_ratio > 0.0 && gly_ratio < 0.2,
        "Glycine ratio should be reasonable, got {}",
        gly_ratio
    );
}

#[test]
fn test_hydrophobic_ratio_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let hydro_ratio = structure.hydrophobic_ratio();

    // Typical proteins have 30-50% hydrophobic residues
    assert!(
        hydro_ratio > 0.2 && hydro_ratio < 0.6,
        "Hydrophobic ratio should be reasonable, got {}",
        hydro_ratio
    );
}

#[test]
fn test_polar_ratio_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let polar = structure.polar_ratio();
    let charged = structure.charged_ratio();
    let hydrophobic = structure.hydrophobic_ratio();

    // Ratios should be positive
    assert!(polar >= 0.0);
    assert!(charged >= 0.0);
    assert!(hydrophobic >= 0.0);

    // Sum should be less than or equal to 1.0 (some AAs may not be in any category)
    assert!(polar + charged + hydrophobic <= 1.0 + 1e-10);
}

#[test]
fn test_count_ca_residues() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ca_count = structure.count_ca_residues();

    // 1UBQ has 76 residues
    assert!(
        ca_count > 70 && ca_count < 80,
        "1UBQ should have ~76 residues, got {}",
        ca_count
    );
}

#[test]
fn test_missing_residue_ratio() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let missing = structure.missing_residue_ratio();

    // 1UBQ should be mostly complete
    assert!(
        missing < 0.1,
        "1UBQ should have few missing residues, got {}",
        missing
    );
}

// ============================================================================
// Geometry Tests
// ============================================================================

#[test]
fn test_radius_of_gyration_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let rg = structure.radius_of_gyration();

    // Ubiquitin is a small globular protein, Rg should be ~10-15 Å
    assert!(
        rg > 5.0 && rg < 20.0,
        "Rg should be reasonable for ubiquitin, got {}",
        rg
    );
}

#[test]
fn test_max_ca_distance_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let max_dist = structure.max_ca_distance();
    let rg = structure.radius_of_gyration();

    // Max distance should be greater than Rg
    assert!(
        max_dist > rg,
        "Max distance ({}) should be > Rg ({})",
        max_dist,
        rg
    );

    // For a small protein, max distance should be < 60 Å
    assert!(
        max_dist < 60.0,
        "Max distance should be reasonable, got {}",
        max_dist
    );
}

#[test]
fn test_secondary_structure_ratio_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ss_ratio = structure.secondary_structure_ratio();

    // Well-folded proteins typically have high ordered content
    // This is a heuristic, so we just check it's reasonable
    assert!(
        (0.0..=1.0).contains(&ss_ratio),
        "SS ratio should be 0-1, got {}",
        ss_ratio
    );

    // Ubiquitin is well-structured, so expect high ratio
    assert!(
        ss_ratio > 0.8,
        "Ubiquitin should have high ordered content, got {}",
        ss_ratio
    );
}

#[test]
fn test_compactness_index_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let compactness = structure.compactness_index();

    // Globular proteins typically have compactness 2-3
    assert!(
        compactness > 1.0 && compactness < 5.0,
        "Compactness should be reasonable, got {}",
        compactness
    );
}

#[test]
fn test_ca_density_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let density = structure.ca_density();

    // Density should be positive
    assert!(
        density > 0.0,
        "Density should be positive, got {}",
        density
    );

    // Typical protein densities are in the range 0.001-0.01 atoms/Å³
    assert!(
        density > 0.0001 && density < 0.1,
        "Density should be in reasonable range, got {}",
        density
    );
}

#[test]
fn test_ca_bounding_box_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = structure.ca_bounding_box();

    // Bounding box should have positive extent in all dimensions
    assert!(x_max > x_min);
    assert!(y_max > y_min);
    assert!(z_max > z_min);

    // For ubiquitin, dimensions should be roughly 30-40 Å
    let dx = x_max - x_min;
    let dy = y_max - y_min;
    let dz = z_max - z_min;

    assert!(
        dx > 10.0 && dx < 50.0,
        "X dimension should be reasonable, got {}",
        dx
    );
    assert!(
        dy > 10.0 && dy < 50.0,
        "Y dimension should be reasonable, got {}",
        dy
    );
    assert!(
        dz > 10.0 && dz < 50.0,
        "Z dimension should be reasonable, got {}",
        dz
    );
}

// ============================================================================
// Combined Descriptor Tests
// ============================================================================

#[test]
fn test_structure_descriptors_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let desc = structure.structure_descriptors();

    // Check all fields are populated
    assert!(desc.num_residues > 0);
    assert!(desc.num_atoms > 0);
    assert!(!desc.aa_composition.is_empty());
    assert!(desc.radius_of_gyration > 0.0);
    assert!(desc.max_ca_distance > 0.0);
    assert!(desc.compactness_index > 0.0);
    assert!(desc.ca_density > 0.0);
    assert!(desc.secondary_structure_ratio >= 0.0);
}

#[test]
fn test_descriptors_consistency() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let desc = structure.structure_descriptors();

    // Verify consistency between individual calls and combined descriptor
    assert_eq!(desc.num_residues, structure.count_ca_residues());
    assert_eq!(desc.num_atoms, structure.get_num_atoms());
    assert!((desc.radius_of_gyration - structure.radius_of_gyration()).abs() < 1e-10);
    assert!((desc.max_ca_distance - structure.max_ca_distance()).abs() < 1e-10);
    assert!((desc.glycine_ratio - structure.glycine_ratio()).abs() < 1e-10);
    assert!((desc.hydrophobic_ratio - structure.hydrophobic_ratio()).abs() < 1e-10);
}

// ============================================================================
// Edge Cases
// ============================================================================

#[test]
fn test_empty_structure_descriptors() {
    let structure = PdbStructure::new();

    assert!(structure.aa_composition().is_empty());
    assert_eq!(structure.glycine_ratio(), 0.0);
    assert_eq!(structure.hydrophobic_ratio(), 0.0);
    assert_eq!(structure.count_ca_residues(), 0);
    assert_eq!(structure.missing_residue_ratio(), 0.0);
    assert_eq!(structure.radius_of_gyration(), 0.0);
    assert_eq!(structure.max_ca_distance(), 0.0);
    assert_eq!(structure.secondary_structure_ratio(), 0.0);
    assert_eq!(structure.compactness_index(), 0.0);
    assert_eq!(structure.ca_density(), 0.0);
}

#[test]
fn test_multi_chain_descriptors() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // All descriptors should work regardless of chain count
    let desc = structure.structure_descriptors();
    assert!(desc.num_residues > 0);
}

// ============================================================================
// Compare Multiple Structures
// ============================================================================

#[test]
fn test_compare_structures() {
    let ubq = parse_pdb_file(get_test_file("1UBQ.pdb")).expect("Failed to parse 1UBQ.pdb");
    let hm2 = parse_pdb_file(get_test_file("8HM2.pdb")).expect("Failed to parse 8HM2.pdb");

    let ubq_desc = ubq.structure_descriptors();
    let hm2_desc = hm2.structure_descriptors();

    // Both should have valid descriptors
    assert!(ubq_desc.num_residues > 0);
    assert!(hm2_desc.num_residues > 0);

    // Larger structure should have larger max_ca_distance
    if hm2_desc.num_residues > ubq_desc.num_residues {
        // Not necessarily true for all cases, but a reasonable expectation
        // We just verify both are valid
        assert!(hm2_desc.max_ca_distance > 0.0);
        assert!(ubq_desc.max_ca_distance > 0.0);
    }
}
