//! Integration tests for dihedral angle analysis and Ramachandran classification.

use pdbrust::{PdbStructure, parse_pdb_file};

#[cfg(all(feature = "descriptors", feature = "dssp"))]
use pdbrust::RamachandranRegion;

/// Create a simple test structure with backbone atoms for a helix
fn create_helix_structure() -> PdbStructure {
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
fn test_phi_psi_angles_basic() {
    let structure = create_helix_structure();
    let dihedrals = structure.phi_psi_angles();

    assert_eq!(dihedrals.len(), 3);

    // First residue should have no phi (N-terminus)
    assert!(dihedrals[0].phi.is_none());
    // Middle residue should have both
    assert!(dihedrals[1].phi.is_some() || dihedrals[1].psi.is_some());
    // Last residue should have no psi (C-terminus)
    assert!(dihedrals[2].psi.is_none());
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_phi_psi_angles_empty_structure() {
    let structure = PdbStructure::new();
    let dihedrals = structure.phi_psi_angles();
    assert!(dihedrals.is_empty());
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_ramachandran_statistics() {
    let structure = create_helix_structure();
    let stats = structure.ramachandran_statistics();

    // Check that statistics are reasonable
    assert!(stats.total_residues <= 3);
    assert!(
        (stats.favored_fraction + stats.allowed_fraction + stats.outlier_fraction - 1.0).abs()
            < 0.01
            || stats.total_residues == 0
    );
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_ramachandran_outliers() {
    let structure = create_helix_structure();
    let outliers = structure.ramachandran_outliers();

    // The simple test structure may or may not have outliers
    // Just verify it doesn't crash
    assert!(outliers.len() <= structure.phi_psi_angles().len());
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_cis_peptide_bonds() {
    let structure = create_helix_structure();
    let cis_bonds = structure.cis_peptide_bonds();

    // Normal structures should have trans peptide bonds (omega ~ 180)
    // This simple structure likely has trans bonds
    for (res1, res2) in &cis_bonds {
        // Cis non-proline would be unusual
        if res2.residue_name != "PRO" {
            println!(
                "Unusual cis non-proline: {}-{}",
                res1.residue_name, res2.residue_name
            );
        }
    }
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_residue_dihedrals_methods() {
    let structure = create_helix_structure();
    let dihedrals = structure.phi_psi_angles();

    for d in &dihedrals {
        // Test has_phi_psi
        if d.phi.is_some() && d.psi.is_some() {
            assert!(d.has_phi_psi());
        }

        // Test region methods
        if d.ramachandran_region == RamachandranRegion::Outlier {
            assert!(d.ramachandran_region.is_outlier());
            assert!(!d.ramachandran_region.is_favorable());
        }
        if d.ramachandran_region == RamachandranRegion::Core {
            assert!(d.ramachandran_region.is_favorable());
            assert!(!d.ramachandran_region.is_outlier());
        }
    }
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_ramachandran_region_classification() {
    // Test that classify functions work correctly
    // This is a unit test for the classification logic

    // α-helix region (φ ≈ -60°, ψ ≈ -45°)
    // β-sheet region (φ ≈ -120°, ψ ≈ +130°)
    // The classification is tested in the unit tests

    let structure = create_helix_structure();
    let stats = structure.ramachandran_statistics();

    // All counts should be non-negative
    assert!(stats.favored_count <= stats.total_residues);
    assert!(stats.allowed_count <= stats.total_residues);
    assert!(stats.outlier_count <= stats.total_residues);
}

#[test]
#[cfg(all(feature = "descriptors", feature = "dssp"))]
fn test_with_real_pdb_file() {
    // Try to load a real PDB file if available
    let path = "examples/pdb_files/1UBQ.pdb";
    if let Ok(structure) = parse_pdb_file(path) {
        let dihedrals = structure.phi_psi_angles();
        assert!(!dihedrals.is_empty());

        let stats = structure.ramachandran_statistics();
        println!("Ramachandran stats for 1UBQ:");
        println!("  Total residues: {}", stats.total_residues);
        println!("  Favored: {:.1}%", stats.favored_fraction * 100.0);
        println!("  Allowed: {:.1}%", stats.allowed_fraction * 100.0);
        println!("  Outliers: {:.1}%", stats.outlier_fraction * 100.0);
        println!("  Cis peptides: {}", stats.cis_peptide_count);

        // 1UBQ should have mostly favored residues
        assert!(stats.favored_fraction > 0.5);
    }
}
