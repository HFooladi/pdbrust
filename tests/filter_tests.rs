//! Integration tests for the filter module.
//!
//! These tests verify the filtering and cleaning functionality
//! using real PDB files.

#![cfg(feature = "filter")]

use pdbrust::{parse_pdb_file, PdbStructure};
use std::path::PathBuf;
use tempfile::NamedTempFile;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

// ============================================================================
// Extraction Tests
// ============================================================================

#[test]
fn test_get_ca_coords_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let ca_coords = structure.get_ca_coords(None);

    // 1UBQ has 76 residues, so should have 76 CA atoms
    assert!(!ca_coords.is_empty(), "Should have CA coordinates");
    assert!(ca_coords.len() > 70, "1UBQ should have ~76 CA atoms, got {}", ca_coords.len());

    // All coordinates should be valid (not NaN or Inf)
    for (x, y, z) in &ca_coords {
        assert!(x.is_finite(), "x should be finite");
        assert!(y.is_finite(), "y should be finite");
        assert!(z.is_finite(), "z should be finite");
    }
}

#[test]
fn test_get_ca_coords_by_chain() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let chains = structure.get_chain_ids();

    // Get CA coords for first chain
    if let Some(first_chain) = chains.first() {
        let chain_ca = structure.get_ca_coords(Some(first_chain));
        let all_ca = structure.get_ca_coords(None);

        // Chain-specific should be <= total
        assert!(chain_ca.len() <= all_ca.len());
    }
}

#[test]
fn test_remove_ligands() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let original_atoms = structure.get_num_atoms();
    let cleaned = structure.remove_ligands();

    // Cleaned structure should have same or fewer atoms
    assert!(cleaned.get_num_atoms() <= original_atoms);

    // All remaining atoms should be from standard residues
    for atom in &cleaned.atoms {
        assert!(
            pdbrust::filter::is_standard_residue(&atom.residue_name),
            "Found non-standard residue after remove_ligands: {}",
            atom.residue_name
        );
    }
}

#[test]
fn test_keep_only_chain() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let chains = structure.get_chain_ids();
    if let Some(first_chain) = chains.first() {
        let single_chain = structure.keep_only_chain(first_chain);

        // Should only have one chain
        assert_eq!(single_chain.get_num_chains(), 1);

        // All atoms should be from that chain
        for atom in &single_chain.atoms {
            assert_eq!(&atom.chain_id, first_chain);
        }
    }
}

#[test]
fn test_keep_only_ca() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let original_atoms = structure.get_num_atoms();
    let ca_only = structure.keep_only_ca();

    // CA-only should have fewer atoms
    assert!(ca_only.get_num_atoms() < original_atoms);

    // All atoms should be CA
    for atom in &ca_only.atoms {
        assert_eq!(atom.name.trim(), "CA", "Expected CA, got {}", atom.name);
    }
}

#[test]
fn test_keep_only_backbone() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let original_atoms = structure.get_num_atoms();
    let backbone = structure.keep_only_backbone();

    // Backbone should have fewer atoms than full structure
    assert!(backbone.get_num_atoms() < original_atoms);

    // Backbone should have more atoms than CA-only
    let ca_only = structure.keep_only_ca();
    assert!(backbone.get_num_atoms() > ca_only.get_num_atoms());

    // All atoms should be backbone atoms
    for atom in &backbone.atoms {
        let name = atom.name.trim();
        assert!(
            ["N", "CA", "C", "O"].contains(&name),
            "Expected backbone atom, got {}",
            name
        );
    }
}

#[test]
fn test_filter_atoms_custom() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Filter by residue name
    let alanines = structure.filter_atoms(|atom| atom.residue_name == "ALA");

    for atom in &alanines.atoms {
        assert_eq!(atom.residue_name, "ALA");
    }
}

// ============================================================================
// Cleaning Tests
// ============================================================================

#[test]
fn test_normalize_chain_ids() {
    let path = get_test_file("1UBQ.pdb");
    let mut structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let original_chains = structure.get_chain_ids();
    structure.normalize_chain_ids();
    let new_chains = structure.get_chain_ids();

    // Should have same number of chains
    assert_eq!(original_chains.len(), new_chains.len());

    // First chain should be 'A'
    if !new_chains.is_empty() {
        assert_eq!(new_chains[0], "A");
    }
}

#[test]
fn test_reindex_residues() {
    let path = get_test_file("1UBQ.pdb");
    let mut structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    structure.reindex_residues();

    // Check that residues start at 1 for each chain
    for chain_id in structure.get_chain_ids() {
        let residues = structure.get_residues_for_chain(&chain_id);
        if !residues.is_empty() {
            assert_eq!(residues[0].0, 1, "First residue should be 1 after reindexing");
        }
    }

    // Check that residues are continuous
    for chain_id in structure.get_chain_ids() {
        let residues = structure.get_residues_for_chain(&chain_id);
        for (i, (seq, _)) in residues.iter().enumerate() {
            assert_eq!(
                *seq as usize,
                i + 1,
                "Residue {} should have seq {} after reindexing",
                i,
                i + 1
            );
        }
    }
}

#[test]
fn test_center_structure() {
    let path = get_test_file("1UBQ.pdb");
    let mut structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Get original centroid
    let (orig_cx, orig_cy, orig_cz) = structure.get_ca_centroid();
    assert!(
        orig_cx.abs() > 1.0 || orig_cy.abs() > 1.0 || orig_cz.abs() > 1.0,
        "Original centroid should not be at origin"
    );

    structure.center_structure();

    // After centering, centroid should be at origin
    let (cx, cy, cz) = structure.get_ca_centroid();
    assert!(
        cx.abs() < 1e-10,
        "Centered cx should be ~0, got {}",
        cx
    );
    assert!(
        cy.abs() < 1e-10,
        "Centered cy should be ~0, got {}",
        cy
    );
    assert!(
        cz.abs() < 1e-10,
        "Centered cz should be ~0, got {}",
        cz
    );
}

#[test]
fn test_renumber_atoms() {
    let path = get_test_file("1UBQ.pdb");
    let mut structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    structure.renumber_atoms();

    // Check that atoms are numbered 1, 2, 3, ...
    for (i, atom) in structure.atoms.iter().enumerate() {
        assert_eq!(
            atom.serial as usize,
            i + 1,
            "Atom {} should have serial {}",
            i,
            i + 1
        );
    }
}

// ============================================================================
// Chained Operations Tests
// ============================================================================

#[test]
fn test_chained_operations() {
    let path = get_test_file("1UBQ.pdb");
    let mut structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Chain multiple operations
    structure
        .normalize_chain_ids()
        .reindex_residues()
        .renumber_atoms()
        .center_structure();

    // Verify all operations were applied
    assert_eq!(structure.get_chain_ids()[0], "A");
    assert_eq!(structure.atoms[0].serial, 1);

    let (cx, cy, cz) = structure.get_ca_centroid();
    assert!(cx.abs() < 1e-10);
    assert!(cy.abs() < 1e-10);
    assert!(cz.abs() < 1e-10);
}

#[test]
fn test_filter_then_clean() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // First filter, then clean
    let mut cleaned = structure.remove_ligands().keep_only_ca();
    cleaned.center_structure();

    // All atoms should be CA atoms from standard residues
    for atom in &cleaned.atoms {
        assert_eq!(atom.name.trim(), "CA");
        assert!(pdbrust::filter::is_standard_residue(&atom.residue_name));
    }

    // Should be centered
    let (cx, cy, cz) = cleaned.get_ca_centroid();
    assert!(cx.abs() < 1e-10);
    assert!(cy.abs() < 1e-10);
    assert!(cz.abs() < 1e-10);
}

// ============================================================================
// Write and Re-read Tests
// ============================================================================

#[test]
fn test_write_filtered_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Filter to CA only
    let ca_only = structure.keep_only_ca();
    let original_ca_count = ca_only.get_num_atoms();

    // Write to temp file
    let temp_file = NamedTempFile::new().expect("Failed to create temp file");
    ca_only
        .to_file(temp_file.path())
        .expect("Failed to write filtered structure");

    // Re-read and verify
    let reread = parse_pdb_file(temp_file.path()).expect("Failed to re-read filtered structure");
    assert_eq!(
        reread.get_num_atoms(),
        original_ca_count,
        "Re-read structure should have same number of atoms"
    );

    // All atoms should still be CA
    for atom in &reread.atoms {
        assert_eq!(atom.name.trim(), "CA");
    }
}

#[test]
fn test_write_cleaned_structure() {
    let path = get_test_file("1UBQ.pdb");
    let mut structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Clean the structure
    structure
        .normalize_chain_ids()
        .reindex_residues()
        .renumber_atoms();

    // Write to temp file
    let temp_file = NamedTempFile::new().expect("Failed to create temp file");
    structure
        .to_file(temp_file.path())
        .expect("Failed to write cleaned structure");

    // Re-read and verify
    let reread = parse_pdb_file(temp_file.path()).expect("Failed to re-read cleaned structure");

    assert_eq!(reread.get_chain_ids()[0], "A");
    assert_eq!(reread.atoms[0].serial, 1);
}

// ============================================================================
// Edge Cases
// ============================================================================

#[test]
fn test_empty_structure_operations() {
    let mut structure = PdbStructure::new();

    // All operations should work on empty structure without panic
    let ca_coords = structure.get_ca_coords(None);
    assert!(ca_coords.is_empty());

    let cleaned = structure.remove_ligands();
    assert!(cleaned.atoms.is_empty());

    let chain_a = structure.keep_only_chain("A");
    assert!(chain_a.atoms.is_empty());

    let ca_only = structure.keep_only_ca();
    assert!(ca_only.atoms.is_empty());

    structure.normalize_chain_ids();
    structure.reindex_residues();
    structure.center_structure();
    structure.renumber_atoms();
}

#[test]
fn test_nonexistent_chain() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Extract nonexistent chain
    let chain_z = structure.keep_only_chain("Z");
    assert!(chain_z.atoms.is_empty());

    // Get CA coords for nonexistent chain
    let ca_coords = structure.get_ca_coords(Some("Z"));
    assert!(ca_coords.is_empty());
}

// ============================================================================
// Multi-model Tests
// ============================================================================

#[test]
fn test_filter_multi_model() {
    let path = get_test_file("multi_model.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse multi_model.pdb");

    // Filtering should work on multi-model structures
    let ca_only = structure.keep_only_ca();

    for atom in &ca_only.atoms {
        assert_eq!(atom.name.trim(), "CA");
    }
}
