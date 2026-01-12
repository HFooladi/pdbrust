//! Integration tests for the geometry module.
//!
//! These tests verify RMSD calculation and structure alignment
//! using real PDB files.

#![cfg(feature = "geometry")]

use pdbrust::geometry::{calculate_alignment, calculate_rmsd, rmsd_from_coords, AtomSelection};
use pdbrust::{parse_pdb_file, PdbStructure};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Basic RMSD Tests
// ============================================================================

#[test]
fn test_rmsd_from_coords_identical() {
    let coords = vec![
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    ];
    let rmsd = rmsd_from_coords(&coords, &coords).unwrap();
    assert!(rmsd < 1e-10, "RMSD of identical coords should be 0, got {}", rmsd);
}

#[test]
fn test_rmsd_from_coords_known_displacement() {
    // All atoms displaced by 2.0 in x direction
    let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)];
    let coords2 = vec![(2.0, 0.0, 0.0), (3.0, 0.0, 0.0), (4.0, 0.0, 0.0)];
    let rmsd = rmsd_from_coords(&coords1, &coords2).unwrap();
    assert!(
        (rmsd - 2.0).abs() < 1e-10,
        "RMSD should be 2.0, got {}",
        rmsd
    );
}

// ============================================================================
// Self-RMSD Tests (structure vs itself)
// ============================================================================

#[test]
fn test_self_rmsd_1ubq() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let rmsd = structure.rmsd_to(&structure).unwrap();
    assert!(rmsd < 1e-10, "Self-RMSD should be 0, got {}", rmsd);
}

#[test]
fn test_self_rmsd_with_selection() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Test with different selections
    let rmsd_ca = structure
        .rmsd_to_with_selection(&structure, AtomSelection::CaOnly)
        .unwrap();
    let rmsd_bb = structure
        .rmsd_to_with_selection(&structure, AtomSelection::Backbone)
        .unwrap();
    let rmsd_all = structure
        .rmsd_to_with_selection(&structure, AtomSelection::AllAtoms)
        .unwrap();

    assert!(rmsd_ca < 1e-10, "Self-RMSD (CA) should be 0");
    assert!(rmsd_bb < 1e-10, "Self-RMSD (backbone) should be 0");
    assert!(rmsd_all < 1e-10, "Self-RMSD (all) should be 0");
}

// ============================================================================
// Alignment Tests
// ============================================================================

#[test]
fn test_self_alignment_1ubq() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let (aligned, result) = structure.align_to(&structure).unwrap();

    assert!(result.rmsd < 1e-10, "Self-alignment RMSD should be 0");
    assert!(result.num_atoms > 0, "Should have aligned some atoms");
    assert_eq!(
        aligned.atoms.len(),
        structure.atoms.len(),
        "Aligned structure should have same number of atoms"
    );
}

#[test]
fn test_alignment_reduces_rmsd() {
    let path = get_test_file("1UBQ.pdb");
    let target = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create a translated copy
    let mut mobile = target.clone();
    for atom in &mut mobile.atoms {
        atom.x += 50.0;
        atom.y += 50.0;
        atom.z += 50.0;
    }

    // RMSD before alignment should be ~86.6 (sqrt(50^2 + 50^2 + 50^2))
    let rmsd_before = mobile.rmsd_to(&target).unwrap();
    assert!(rmsd_before > 80.0, "Pre-alignment RMSD should be large");

    // Align
    let (aligned, result) = mobile.align_to(&target).unwrap();

    // RMSD after alignment should be ~0
    assert!(
        result.rmsd < 1e-6,
        "Post-alignment RMSD should be ~0, got {}",
        result.rmsd
    );

    // Verify aligned coordinates match target
    for (aligned_atom, target_atom) in aligned.atoms.iter().zip(target.atoms.iter()) {
        let dx = (aligned_atom.x - target_atom.x).abs();
        let dy = (aligned_atom.y - target_atom.y).abs();
        let dz = (aligned_atom.z - target_atom.z).abs();
        assert!(
            dx < 1e-4 && dy < 1e-4 && dz < 1e-4,
            "Aligned atoms should match target: ({}, {}, {})",
            dx,
            dy,
            dz
        );
    }
}

#[test]
fn test_alignment_with_rotation() {
    let path = get_test_file("1UBQ.pdb");
    let target = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create a rotated copy (90 degree rotation around z-axis)
    let mut mobile = target.clone();
    for atom in &mut mobile.atoms {
        let x = atom.x;
        let y = atom.y;
        atom.x = -y;
        atom.y = x;
    }

    // Align
    let (_aligned, result) = mobile.align_to(&target).unwrap();

    // RMSD after alignment should be ~0 (rotation is reversible)
    assert!(
        result.rmsd < 1e-4,
        "Post-alignment RMSD should be ~0, got {}",
        result.rmsd
    );

    // Rotation matrix should have determinant 1 (proper rotation)
    let det = result.rotation[0][0]
        * (result.rotation[1][1] * result.rotation[2][2] - result.rotation[1][2] * result.rotation[2][1])
        - result.rotation[0][1]
            * (result.rotation[1][0] * result.rotation[2][2] - result.rotation[1][2] * result.rotation[2][0])
        + result.rotation[0][2]
            * (result.rotation[1][0] * result.rotation[2][1] - result.rotation[1][1] * result.rotation[2][0]);
    assert!(
        (det - 1.0).abs() < 1e-6,
        "Rotation matrix determinant should be 1, got {}",
        det
    );
}

#[test]
fn test_calculate_alignment() {
    let path = get_test_file("1UBQ.pdb");
    let target = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let mut mobile = target.clone();
    for atom in &mut mobile.atoms {
        atom.x += 10.0;
        atom.y += 20.0;
        atom.z += 30.0;
    }

    let result = calculate_alignment(&mobile, &target, AtomSelection::CaOnly).unwrap();

    assert!(result.rmsd < 1e-6, "Aligned RMSD should be ~0");
    assert!(result.num_atoms > 0, "Should have used some atoms");
}

// ============================================================================
// Per-Residue RMSD Tests
// ============================================================================

#[test]
fn test_per_residue_rmsd_self() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let per_res = structure.per_residue_rmsd_to(&structure).unwrap();

    assert!(!per_res.is_empty(), "Should have per-residue RMSDs");

    // All per-residue RMSDs should be 0 for self-comparison
    for r in &per_res {
        assert!(
            r.rmsd < 1e-10,
            "Per-residue self-RMSD should be 0, got {} for {}{}",
            r.rmsd,
            r.residue_id.0,
            r.residue_id.1
        );
    }
}

#[test]
fn test_per_residue_rmsd_translated() {
    let path = get_test_file("1UBQ.pdb");
    let target = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Create translated copy
    let mut mobile = target.clone();
    for atom in &mut mobile.atoms {
        atom.x += 100.0;
    }

    let per_res = mobile.per_residue_rmsd_to(&target).unwrap();

    // After alignment, all per-residue RMSDs should be ~0
    for r in &per_res {
        assert!(
            r.rmsd < 1e-4,
            "Per-residue RMSD after alignment should be ~0, got {} for {}{}",
            r.rmsd,
            r.residue_id.0,
            r.residue_id.1
        );
    }
}

#[test]
fn test_per_residue_rmsd_properties() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let per_res = structure.per_residue_rmsd_to(&structure).unwrap();

    // Verify each result has valid properties
    for r in &per_res {
        // RMSD should be non-negative
        assert!(r.rmsd >= 0.0, "RMSD should be non-negative");

        // Should have at least one atom
        assert!(r.num_atoms > 0, "Should have at least one atom per residue");

        // Residue name should not be empty
        assert!(!r.residue_name.is_empty(), "Residue name should not be empty");
    }
}

// ============================================================================
// Atom Selection Tests
// ============================================================================

#[test]
fn test_atom_selection_matches() {
    assert!(AtomSelection::CaOnly.matches("CA"));
    assert!(AtomSelection::CaOnly.matches(" CA "));
    assert!(!AtomSelection::CaOnly.matches("CB"));

    assert!(AtomSelection::Backbone.matches("N"));
    assert!(AtomSelection::Backbone.matches("CA"));
    assert!(AtomSelection::Backbone.matches("C"));
    assert!(AtomSelection::Backbone.matches("O"));
    assert!(!AtomSelection::Backbone.matches("CB"));

    assert!(AtomSelection::AllAtoms.matches("anything"));

    let custom = AtomSelection::Custom(vec!["CA".to_string(), "CB".to_string()]);
    assert!(custom.matches("CA"));
    assert!(custom.matches("CB"));
    assert!(!custom.matches("N"));
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_error_on_empty_structure() {
    let empty = PdbStructure::new();
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // RMSD with empty structure should fail
    let result = empty.rmsd_to(&structure);
    assert!(result.is_err(), "RMSD with empty structure should fail");

    let result = structure.rmsd_to(&empty);
    assert!(result.is_err(), "RMSD with empty structure should fail");
}

#[test]
fn test_error_on_mismatched_structures() {
    let path1 = get_test_file("1UBQ.pdb");
    let path2 = get_test_file("8HM2.pdb");

    let structure1 = parse_pdb_file(&path1).expect("Failed to parse 1UBQ.pdb");
    let structure2 = parse_pdb_file(&path2).expect("Failed to parse 8HM2.pdb");

    // These structures have different numbers of residues
    let result = calculate_rmsd(&structure1, &structure2, AtomSelection::CaOnly);
    assert!(
        result.is_err(),
        "RMSD between different structures should fail"
    );
}

#[test]
fn test_error_on_insufficient_atoms() {
    // Create structure with only 2 atoms
    let mut structure = PdbStructure::new();
    structure.atoms.push(test_helpers::create_test_atom(0.0, 0.0, 0.0, 1));
    structure.atoms.push(test_helpers::create_test_atom(1.0, 0.0, 0.0, 2));

    let result = structure.align_to(&structure);
    assert!(
        result.is_err(),
        "Alignment with < 3 atoms should fail"
    );
}

// ============================================================================
// AlignmentResult Tests
// ============================================================================

#[test]
fn test_alignment_result_properties() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let (_, result) = structure.align_to(&structure).unwrap();

    // For self-alignment, rotation should be identity matrix
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (result.rotation[i][j] - expected).abs() < 1e-6,
                "Identity rotation expected at [{},{}]",
                i,
                j
            );
        }
    }

    // Translation should be ~0 for self-alignment
    for (i, t) in result.translation.iter().enumerate() {
        assert!(
            t.abs() < 1e-6,
            "Translation[{}] should be ~0, got {}",
            i,
            t
        );
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

mod test_helpers {
    use pdbrust::records::Atom;

    pub fn create_test_atom(x: f64, y: f64, z: f64, residue_seq: i32) -> Atom {
        Atom {
            serial: residue_seq,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        }
    }
}
