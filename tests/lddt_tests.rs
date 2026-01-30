//! Integration tests for LDDT (Local Distance Difference Test) functionality.
//!
//! These tests verify LDDT calculation using real PDB files.

#![cfg(feature = "geometry")]

use pdbrust::geometry::{AtomSelection, LddtOptions};
use pdbrust::{PdbStructure, parse_pdb_file};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Self-LDDT Tests (structure vs itself)
// ============================================================================

#[test]
fn test_lddt_self_1ubq() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let result = structure.lddt_to(&structure).unwrap();

    // Self-comparison should have perfect LDDT
    assert!(
        (result.score - 1.0).abs() < 1e-10,
        "Self-LDDT should be 1.0, got {}",
        result.score
    );
    assert!(result.num_pairs > 0, "Should have distance pairs");
    assert!(result.num_residues > 0, "Should have residues");
}

#[test]
fn test_lddt_self_with_selection() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Test with different selections
    let result_ca = structure
        .lddt_to_with_options(&structure, AtomSelection::CaOnly, LddtOptions::default())
        .unwrap();
    let result_bb = structure
        .lddt_to_with_options(&structure, AtomSelection::Backbone, LddtOptions::default())
        .unwrap();

    assert!(
        (result_ca.score - 1.0).abs() < 1e-10,
        "Self-LDDT (CA) should be 1.0"
    );
    assert!(
        (result_bb.score - 1.0).abs() < 1e-10,
        "Self-LDDT (backbone) should be 1.0"
    );

    // Backbone should have more pairs than CA only
    assert!(
        result_bb.num_pairs >= result_ca.num_pairs,
        "Backbone should have >= pairs than CA"
    );
}

// ============================================================================
// Translation Invariance Tests
// ============================================================================

#[test]
fn test_lddt_translation_invariance_1ubq() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create translated copy
    let mut model = reference.clone();
    for atom in &mut model.atoms {
        atom.x += 100.0;
        atom.y += 200.0;
        atom.z += 300.0;
    }

    let result = model.lddt_to(&reference).unwrap();

    // LDDT should be invariant to translation
    assert!(
        (result.score - 1.0).abs() < 1e-10,
        "LDDT should be translation invariant, got {}",
        result.score
    );
}

// ============================================================================
// Rotation Invariance Tests
// ============================================================================

#[test]
fn test_lddt_rotation_invariance_1ubq() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create rotated copy (90 degrees around z-axis)
    let mut model = reference.clone();
    for atom in &mut model.atoms {
        let x = atom.x;
        let y = atom.y;
        atom.x = -y;
        atom.y = x;
    }

    let result = model.lddt_to(&reference).unwrap();

    // LDDT should be invariant to rotation
    assert!(
        (result.score - 1.0).abs() < 1e-10,
        "LDDT should be rotation invariant, got {}",
        result.score
    );
}

#[test]
fn test_lddt_combined_transform_invariance() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create copy with rotation + translation
    let mut model = reference.clone();
    let cos_45 = std::f64::consts::FRAC_PI_4.cos();
    let sin_45 = std::f64::consts::FRAC_PI_4.sin();

    for atom in &mut model.atoms {
        // Rotate 45 degrees around z-axis
        let x = atom.x;
        let y = atom.y;
        atom.x = x * cos_45 - y * sin_45 + 50.0;
        atom.y = x * sin_45 + y * cos_45 + 100.0;
        atom.z += 75.0;
    }

    let result = model.lddt_to(&reference).unwrap();

    // LDDT should be invariant to rigid body transformation
    assert!(
        (result.score - 1.0).abs() < 1e-10,
        "LDDT should be invariant to rotation+translation, got {}",
        result.score
    );
}

// ============================================================================
// Perturbed Structure Tests
// ============================================================================

#[test]
fn test_lddt_small_perturbation() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Small perturbation (0.3 Angstrom) on CA atoms specifically
    let mut model = reference.clone();
    let mut ca_count = 0;
    for atom in model.atoms.iter_mut() {
        if atom.name.trim() == "CA" {
            if ca_count % 3 == 0 {
                atom.y += 0.3;
            }
            ca_count += 1;
        }
    }

    let result = model.lddt_to(&reference).unwrap();

    // Small perturbation should still have high LDDT
    assert!(
        result.score > 0.9,
        "Small perturbation should have LDDT > 0.9, got {}",
        result.score
    );
    // Note: With very small perturbations (<0.5A) and only affecting some atoms,
    // the LDDT might still be 1.0 as distances might not change enough.
    // So we just check it's >= 0.9, not strictly < 1.0
}

#[test]
fn test_lddt_large_perturbation() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Large perturbation (5 Angstroms)
    let mut model = reference.clone();
    for (i, atom) in model.atoms.iter_mut().enumerate() {
        if i % 3 == 0 {
            atom.y += 5.0;
        }
    }

    let result = model.lddt_to(&reference).unwrap();

    // Large perturbation should have lower LDDT
    assert!(
        result.score < 0.9,
        "Large perturbation should have LDDT < 0.9, got {}",
        result.score
    );
    assert!(
        result.score > 0.0,
        "LDDT should still be > 0.0, got {}",
        result.score
    );
}

// ============================================================================
// Custom Options Tests
// ============================================================================

#[test]
fn test_lddt_custom_thresholds() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Small perturbation
    let mut model = reference.clone();
    model.atoms[10].y += 0.7; // 0.7 Angstrom perturbation

    let options_default = LddtOptions::default();
    let options_strict = LddtOptions::default().with_thresholds(vec![0.1, 0.2, 0.3]);

    let result_default = model
        .lddt_to_with_options(&reference, AtomSelection::CaOnly, options_default)
        .unwrap();
    let result_strict = model
        .lddt_to_with_options(&reference, AtomSelection::CaOnly, options_strict)
        .unwrap();

    // Stricter thresholds should give lower or equal score
    assert!(
        result_strict.score <= result_default.score,
        "Stricter thresholds should not give higher LDDT: {} vs {}",
        result_strict.score,
        result_default.score
    );
}

#[test]
fn test_lddt_custom_inclusion_radius() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let options_small = LddtOptions::default().with_inclusion_radius(8.0);
    let options_large = LddtOptions::default().with_inclusion_radius(20.0);

    let result_small = structure
        .lddt_to_with_options(&structure, AtomSelection::CaOnly, options_small)
        .unwrap();
    let result_large = structure
        .lddt_to_with_options(&structure, AtomSelection::CaOnly, options_large)
        .unwrap();

    // Smaller radius should have fewer pairs
    assert!(
        result_small.num_pairs <= result_large.num_pairs,
        "Smaller radius should have <= pairs: {} vs {}",
        result_small.num_pairs,
        result_large.num_pairs
    );

    // Both should still be 1.0 for self-comparison
    assert!((result_small.score - 1.0).abs() < 1e-10);
    assert!((result_large.score - 1.0).abs() < 1e-10);
}

// ============================================================================
// Per-Residue LDDT Tests
// ============================================================================

#[test]
fn test_per_residue_lddt_1ubq() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let per_res = structure.per_residue_lddt_to(&structure).unwrap();

    // 1UBQ has 76 residues
    assert!(!per_res.is_empty(), "Should have residues");

    // All residues should have LDDT = 1.0 for self-comparison
    for r in &per_res {
        assert!(
            (r.score - 1.0).abs() < 1e-10,
            "Self per-residue LDDT should be 1.0, got {} for {}{}",
            r.score,
            r.residue_id.0,
            r.residue_id.1
        );
    }
}

#[test]
fn test_per_residue_lddt_perturbed() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Perturb a specific region
    let mut model = reference.clone();
    let target_residue = 40; // Perturb residue 40

    for atom in &mut model.atoms {
        if atom.residue_seq == target_residue {
            atom.y += 5.0;
        }
    }

    let per_res = model.per_residue_lddt_to(&reference).unwrap();

    // Find the perturbed residue
    let perturbed = per_res.iter().find(|r| r.residue_id.1 == target_residue);

    if let Some(p) = perturbed {
        // The perturbed residue should have lower LDDT than average
        let avg_others: f64 = per_res
            .iter()
            .filter(|r| r.residue_id.1 != target_residue)
            .map(|r| r.score)
            .sum::<f64>()
            / (per_res.len() - 1) as f64;

        assert!(
            p.score < avg_others,
            "Perturbed residue should have lower LDDT: {} vs avg {}",
            p.score,
            avg_others
        );
    }
}

#[test]
fn test_per_residue_lddt_sorted() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let per_res = structure.per_residue_lddt_to(&structure).unwrap();

    // Check that results are sorted by chain_id then residue_seq
    for i in 1..per_res.len() {
        let prev = &per_res[i - 1];
        let curr = &per_res[i];

        let is_sorted = prev.residue_id.0 < curr.residue_id.0
            || (prev.residue_id.0 == curr.residue_id.0 && prev.residue_id.1 <= curr.residue_id.1);

        assert!(is_sorted, "Per-residue results should be sorted");
    }
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_lddt_empty_structure() {
    let structure = PdbStructure::new();

    let result = structure.lddt_to(&structure);

    assert!(result.is_err(), "Empty structure should return error");
}

#[test]
fn test_lddt_mismatched_structures() {
    let path = get_test_file("1UBQ.pdb");
    let structure1 = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create structure with different number of CA atoms
    // by removing all atoms from a specific residue
    let mut structure2 = structure1.clone();
    let target_residue = structure2
        .atoms
        .iter()
        .find(|a| a.name.trim() == "CA")
        .map(|a| a.residue_seq)
        .unwrap();
    structure2.atoms.retain(|a| a.residue_seq != target_residue);

    let result = structure1.lddt_to(&structure2);

    assert!(result.is_err(), "Mismatched structures should return error");
}

// ============================================================================
// Per-Threshold Score Tests
// ============================================================================

#[test]
fn test_lddt_per_threshold_scores() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create perturbed model with varying deviations
    let mut model = reference.clone();
    for (i, atom) in model.atoms.iter_mut().enumerate() {
        // Add perturbations of different sizes
        if i % 10 == 0 {
            atom.y += 0.7; // Within 1.0 threshold
        }
        if i % 15 == 0 {
            atom.z += 1.5; // Within 2.0 threshold
        }
    }

    let result = model.lddt_to(&reference).unwrap();

    // Per-threshold scores should be monotonically increasing
    // (larger thresholds are more lenient)
    assert_eq!(
        result.per_threshold_scores.len(),
        4,
        "Should have 4 threshold scores"
    );

    for i in 1..result.per_threshold_scores.len() {
        assert!(
            result.per_threshold_scores[i] >= result.per_threshold_scores[i - 1],
            "Per-threshold scores should be monotonically increasing: {:?}",
            result.per_threshold_scores
        );
    }
}

// ============================================================================
// Comparison with RMSD Tests
// ============================================================================

#[test]
fn test_lddt_vs_rmsd_translated() {
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Translate model
    let mut model = reference.clone();
    for atom in &mut model.atoms {
        atom.x += 100.0;
        atom.y += 100.0;
        atom.z += 100.0;
    }

    // RMSD without alignment should be large
    let rmsd = model.rmsd_to(&reference).unwrap();
    assert!(rmsd > 100.0, "RMSD without alignment should be large");

    // LDDT should still be 1.0 (superposition-free)
    let lddt = model.lddt_to(&reference).unwrap();
    assert!(
        (lddt.score - 1.0).abs() < 1e-10,
        "LDDT should be 1.0 for translated structure"
    );
}
