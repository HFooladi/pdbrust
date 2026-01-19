//! Integration tests for the selection language.

#![cfg(feature = "filter")]

use pdbrust::{PdbStructure, SelectionError, parse_pdb_file};
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

fn load_test_structure() -> PdbStructure {
    parse_pdb_file(get_test_file("1UBQ.pdb")).expect("Failed to load test structure")
}

// ============================================================================
// Basic Selector Tests
// ============================================================================

#[test]
fn test_select_chain() {
    let structure = load_test_structure();
    let selected = structure.select("chain A").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(atom.chain_id, "A", "All atoms should be from chain A");
    }
}

#[test]
fn test_select_name_ca() {
    let structure = load_test_structure();
    let selected = structure.select("name CA").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(atom.name.trim(), "CA", "All atoms should be CA");
    }
}

#[test]
fn test_select_resname() {
    let structure = load_test_structure();
    let selected = structure.select("resname MET").unwrap();

    // 1UBQ has at least one MET residue (Met1)
    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(
            atom.residue_name.trim(),
            "MET",
            "All atoms should be from MET residues"
        );
    }
}

#[test]
fn test_select_resid() {
    let structure = load_test_structure();
    let selected = structure.select("resid 1").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(atom.residue_seq, 1, "All atoms should be from residue 1");
    }
}

#[test]
fn test_select_resid_range() {
    let structure = load_test_structure();
    let selected = structure.select("resid 1:10").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert!(
            atom.residue_seq >= 1 && atom.residue_seq <= 10,
            "All atoms should be from residues 1-10, got {}",
            atom.residue_seq
        );
    }
}

#[test]
fn test_select_element() {
    let structure = load_test_structure();
    let selected = structure.select("element N").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(
            atom.element.trim().to_uppercase(),
            "N",
            "All atoms should be nitrogen"
        );
    }
}

// ============================================================================
// Keyword Tests
// ============================================================================

#[test]
fn test_select_backbone() {
    let structure = load_test_structure();
    let selected = structure.select("backbone").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert!(
            atom.is_backbone(),
            "All atoms should be backbone atoms, got {}",
            atom.name
        );
    }
}

#[test]
fn test_select_protein() {
    let structure = load_test_structure();
    let selected = structure.select("protein").unwrap();

    // 1UBQ is a protein, so most atoms should match
    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert!(
            pdbrust::filter::is_standard_amino_acid(&atom.residue_name),
            "All atoms should be from standard amino acids, got {}",
            atom.residue_name
        );
    }
}

#[test]
fn test_select_all() {
    let structure = load_test_structure();
    let selected_all = structure.select("all").unwrap();
    let selected_star = structure.select("*").unwrap();

    assert_eq!(
        selected_all.atoms.len(),
        structure.atoms.len(),
        "all should select all atoms"
    );
    assert_eq!(
        selected_star.atoms.len(),
        structure.atoms.len(),
        "* should select all atoms"
    );
}

// ============================================================================
// Boolean Operator Tests
// ============================================================================

#[test]
fn test_select_and() {
    let structure = load_test_structure();
    let selected = structure.select("chain A and name CA").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(atom.chain_id, "A", "All atoms should be from chain A");
        assert_eq!(atom.name.trim(), "CA", "All atoms should be CA");
    }
}

#[test]
fn test_select_or() {
    let structure = load_test_structure();
    let selected = structure.select("resname MET or resname GLN").unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        let resname = atom.residue_name.trim();
        assert!(
            resname == "MET" || resname == "GLN",
            "All atoms should be from MET or GLN residues, got {}",
            resname
        );
    }
}

#[test]
fn test_select_not() {
    let structure = load_test_structure();
    let selected = structure.select("not hydrogen").unwrap();

    for atom in &selected.atoms {
        assert!(!atom.is_hydrogen(), "No atoms should be hydrogen");
    }
}

#[test]
fn test_select_parentheses() {
    let structure = load_test_structure();
    // (resid 1 or resid 2) and name CA
    let selected = structure
        .select("(resid 1 or resid 2) and name CA")
        .unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert!(
            atom.residue_seq == 1 || atom.residue_seq == 2,
            "Residue should be 1 or 2, got {}",
            atom.residue_seq
        );
        assert_eq!(atom.name.trim(), "CA", "All atoms should be CA");
    }
}

#[test]
fn test_select_complex() {
    let structure = load_test_structure();
    let selected = structure
        .select("chain A and backbone and not hydrogen")
        .unwrap();

    assert!(!selected.atoms.is_empty(), "Selection should not be empty");
    for atom in &selected.atoms {
        assert_eq!(atom.chain_id, "A", "All atoms should be from chain A");
        assert!(atom.is_backbone(), "All atoms should be backbone");
        assert!(!atom.is_hydrogen(), "No atoms should be hydrogen");
    }
}

// ============================================================================
// Numeric Comparison Tests
// ============================================================================

#[test]
fn test_select_bfactor_lt() {
    let structure = load_test_structure();
    let selected = structure.select("bfactor < 20.0").unwrap();

    for atom in &selected.atoms {
        assert!(
            atom.temp_factor < 20.0,
            "B-factor should be < 20.0, got {}",
            atom.temp_factor
        );
    }
}

#[test]
fn test_select_bfactor_ge() {
    let structure = load_test_structure();
    let selected = structure.select("bfactor >= 10.0").unwrap();

    for atom in &selected.atoms {
        assert!(
            atom.temp_factor >= 10.0,
            "B-factor should be >= 10.0, got {}",
            atom.temp_factor
        );
    }
}

#[test]
fn test_select_occupancy() {
    let structure = load_test_structure();
    let selected = structure.select("occupancy > 0.5").unwrap();

    for atom in &selected.atoms {
        assert!(
            atom.occupancy > 0.5,
            "Occupancy should be > 0.5, got {}",
            atom.occupancy
        );
    }
}

// ============================================================================
// Alternative Keyword Tests
// ============================================================================

#[test]
fn test_select_alternative_keywords() {
    let structure = load_test_structure();

    // Test various alternative forms
    let _ = structure.select("chainid A").unwrap();
    let _ = structure.select("atomname CA").unwrap();
    let _ = structure.select("resn MET").unwrap();
    let _ = structure.select("resi 1").unwrap();
    let _ = structure.select("elem N").unwrap();
    let _ = structure.select("bb").unwrap(); // backbone
}

#[test]
fn test_select_case_insensitive() {
    let structure = load_test_structure();

    // Keywords should be case-insensitive
    let upper = structure.select("CHAIN A AND NAME CA").unwrap();
    let lower = structure.select("chain a and name ca").unwrap();

    assert_eq!(
        upper.atoms.len(),
        lower.atoms.len(),
        "Case should not matter for keywords"
    );
}

// ============================================================================
// Equivalence Tests (compare with existing filter methods)
// ============================================================================

#[test]
fn test_select_equivalence_chain() {
    let structure = load_test_structure();

    let via_select = structure.select("chain A").unwrap();
    let via_filter = structure.keep_only_chain("A");

    assert_eq!(
        via_select.atoms.len(),
        via_filter.atoms.len(),
        "select('chain A') should be equivalent to keep_only_chain('A')"
    );
}

#[test]
fn test_select_equivalence_ca() {
    let structure = load_test_structure();

    let via_select = structure.select("name CA").unwrap();
    let via_filter = structure.keep_only_ca();

    assert_eq!(
        via_select.atoms.len(),
        via_filter.atoms.len(),
        "select('name CA') should be equivalent to keep_only_ca()"
    );
}

#[test]
fn test_select_equivalence_backbone() {
    let structure = load_test_structure();

    let via_select = structure.select("backbone").unwrap();
    let via_filter = structure.keep_only_backbone();

    assert_eq!(
        via_select.atoms.len(),
        via_filter.atoms.len(),
        "select('backbone') should be equivalent to keep_only_backbone()"
    );
}

#[test]
fn test_select_equivalence_protein() {
    let structure = load_test_structure();

    let via_select = structure.select("protein").unwrap();
    let via_filter = structure.remove_ligands();

    assert_eq!(
        via_select.atoms.len(),
        via_filter.atoms.len(),
        "select('protein') should be equivalent to remove_ligands() for a pure protein"
    );
}

#[test]
fn test_select_equivalence_chain_and_ca() {
    let structure = load_test_structure();

    let via_select = structure.select("chain A and name CA").unwrap();
    let via_filter = structure.keep_only_chain("A").keep_only_ca();

    assert_eq!(
        via_select.atoms.len(),
        via_filter.atoms.len(),
        "select('chain A and name CA') should be equivalent to keep_only_chain('A').keep_only_ca()"
    );
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_select_error_empty() {
    let structure = load_test_structure();
    let result = structure.select("");

    assert!(matches!(result, Err(SelectionError::EmptySelection)));
}

#[test]
fn test_select_error_unclosed_paren() {
    let structure = load_test_structure();
    let result = structure.select("(chain A and name CA");

    assert!(matches!(
        result,
        Err(SelectionError::UnclosedParenthesis { .. })
    ));
}

#[test]
fn test_select_error_unknown_keyword() {
    let structure = load_test_structure();
    let result = structure.select("unknown_keyword");

    assert!(matches!(result, Err(SelectionError::UnknownKeyword { .. })));
}

#[test]
fn test_select_error_missing_value() {
    let structure = load_test_structure();
    let result = structure.select("chain and name CA");

    // "chain and" should fail because "and" is parsed as the chain value
    // then "name" becomes an unknown keyword (we need an operator)
    assert!(result.is_err());
}

// ============================================================================
// Validation Tests
// ============================================================================

#[test]
fn test_validate_selection_valid() {
    assert!(PdbStructure::validate_selection("chain A").is_ok());
    assert!(PdbStructure::validate_selection("chain A and name CA").is_ok());
    assert!(PdbStructure::validate_selection("(chain A or chain B) and backbone").is_ok());
    assert!(PdbStructure::validate_selection("resid 1:100 and protein").is_ok());
    assert!(PdbStructure::validate_selection("bfactor < 30.0").is_ok());
}

#[test]
fn test_validate_selection_invalid() {
    assert!(PdbStructure::validate_selection("").is_err());
    assert!(PdbStructure::validate_selection("(chain A").is_err());
    assert!(PdbStructure::validate_selection("unknown_keyword").is_err());
    assert!(PdbStructure::validate_selection("chain and").is_err());
}
