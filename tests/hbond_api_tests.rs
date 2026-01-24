//! Integration tests for hydrogen bond network analysis.

use pdbrust::{PdbStructure, parse_pdb_file};

#[cfg(all(feature = "descriptors", feature = "dssp"))]
use pdbrust::HBondType;

/// Create a test structure with residues that should form H-bonds
fn create_test_structure() -> PdbStructure {
    // This creates a minimal structure - real H-bonds need proper geometry
    let pdb_content = r#"
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00 20.00           O
ATOM      5  N   ALA A   2       3.320   1.567   0.000  1.00 20.00           N
ATOM      6  CA  ALA A   2       3.954   2.881   0.000  1.00 20.00           C
ATOM      7  C   ALA A   2       5.464   2.771   0.000  1.00 20.00           C
ATOM      8  O   ALA A   2       6.108   1.726   0.000  1.00 20.00           O
ATOM      9  N   ALA A   3       6.012   3.973   0.000  1.00 20.00           N
ATOM     10  CA  ALA A   3       7.445   4.145   0.000  1.00 20.00           C
ATOM     11  C   ALA A   3       8.115   2.811   0.000  1.00 20.00           C
ATOM     12  O   ALA A   3       7.418   1.800   0.000  1.00 20.00           O
END
"#;
    pdbrust::parse_pdb_string(pdb_content).unwrap()
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_mainchain_hbonds_basic() {
    let structure = create_test_structure();
    let hbonds = structure.mainchain_hbonds();

    // The simple test structure may not have H-bonds due to geometry
    // Just verify it doesn't crash
    println!("Found {} H-bonds in test structure", hbonds.len());

    for hb in &hbonds {
        // Verify all fields are populated
        assert!(!hb.donor_chain.is_empty());
        assert!(!hb.acceptor_chain.is_empty());
        assert!(hb.energy < 0.0); // H-bond energy should be negative
        assert!(hb.n_o_distance > 0.0); // Distance should be positive
    }
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_mainchain_hbonds_empty_structure() {
    let structure = PdbStructure::new();
    let hbonds = structure.mainchain_hbonds();
    assert!(hbonds.is_empty());
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_hbonds_for_residue() {
    let structure = create_test_structure();

    // Query H-bonds for residue 2
    let res_hbonds = structure.hbonds_for_residue("A", 2);

    // Check that the result is valid
    assert!(res_hbonds.total() <= structure.mainchain_hbonds().len() * 2);

    // Test methods
    if res_hbonds.total() > 0 {
        assert!(res_hbonds.has_hbonds());
    }
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_hbond_statistics() {
    let structure = create_test_structure();
    let stats = structure.hbond_statistics();

    // Verify statistics are consistent
    assert!(
        stats.intra_helical + stats.beta_sheet + stats.turn + stats.long_range + stats.inter_chain
            <= stats.total_hbonds
    );

    if stats.total_hbonds > 0 {
        // Mean energy should be negative for H-bonds
        assert!(stats.mean_energy < 0.0);
    }
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_hbond_type_classification() {
    // Test H-bond type methods
    assert!(HBondType::IntraHelical.is_secondary_structure());
    assert!(HBondType::BetaSheet.is_secondary_structure());
    assert!(HBondType::Turn.is_secondary_structure());
    assert!(!HBondType::LongRange.is_secondary_structure());
    assert!(!HBondType::InterChain.is_secondary_structure());
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_mainchain_hbond_methods() {
    let structure = create_test_structure();
    let hbonds = structure.mainchain_hbonds();

    for hb in &hbonds {
        // Test is_strong (energy < -1.0)
        if hb.energy < -1.0 {
            assert!(hb.is_strong());
        }

        // Test is_helical
        if matches!(hb.hbond_type, HBondType::IntraHelical) {
            assert!(hb.is_helical());
            assert!(!hb.is_beta_sheet());
        }

        // Test is_beta_sheet
        if matches!(hb.hbond_type, HBondType::BetaSheet) {
            assert!(hb.is_beta_sheet());
            assert!(!hb.is_helical());
        }
    }
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_with_real_pdb_file() {
    // Try to load a real PDB file if available
    let path = "examples/pdb_files/1UBQ.pdb";
    if let Ok(structure) = parse_pdb_file(path) {
        let _hbonds = structure.mainchain_hbonds();
        let stats = structure.hbond_statistics();

        println!("H-bond analysis for 1UBQ:");
        println!("  Total H-bonds: {}", stats.total_hbonds);
        println!("  Intra-helical: {}", stats.intra_helical);
        println!("  Beta-sheet: {}", stats.beta_sheet);
        println!("  Turns: {}", stats.turn);
        println!("  Long-range: {}", stats.long_range);
        println!("  Mean energy: {:.2} kcal/mol", stats.mean_energy);

        // 1UBQ has a helix and beta sheet, should have both types
        assert!(stats.total_hbonds > 0);

        // Test getting H-bonds for a specific residue
        let res_hbonds = structure.hbonds_for_residue("A", 10);
        println!(
            "  Residue A10 H-bonds: donated={}, accepted={}",
            res_hbonds.donated.len(),
            res_hbonds.accepted.len()
        );
    }
}
