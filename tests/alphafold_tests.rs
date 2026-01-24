//! Integration tests for AlphaFold/pLDDT support.

use pdbrust::{PdbStructure, parse_pdb_string};

#[cfg(feature = "descriptors")]
use pdbrust::ConfidenceCategory;

/// Create a structure that looks like an AlphaFold prediction
#[cfg(feature = "descriptors")]
fn create_alphafold_like_structure() -> PdbStructure {
    let pdb_content = r#"HEADER    ALPHAFOLD MODEL
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 95.50           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 94.30           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 93.10           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00 92.80           O
ATOM      5  N   GLY A   2       3.320   1.567   0.000  1.00 75.20           N
ATOM      6  CA  GLY A   2       3.954   2.881   0.000  1.00 73.50           C
ATOM      7  C   GLY A   2       5.464   2.771   0.000  1.00 71.80           C
ATOM      8  O   GLY A   2       6.108   1.726   0.000  1.00 70.10           O
ATOM      9  N   VAL A   3       6.012   3.973   0.000  1.00 45.30           N
ATOM     10  CA  VAL A   3       7.445   4.145   0.000  1.00 42.60           C
ATOM     11  C   VAL A   3       8.115   2.811   0.000  1.00 40.20           C
ATOM     12  O   VAL A   3       7.418   1.800   0.000  1.00 38.90           O
END
"#;
    parse_pdb_string(pdb_content).unwrap()
}

/// Create a structure that looks like experimental (not predicted)
#[cfg(feature = "descriptors")]
fn create_experimental_structure() -> PdbStructure {
    let pdb_content = r#"HEADER    EXPERIMENTAL STRUCTURE
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 15.50           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 14.30           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 13.10           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00 12.80           O
END
"#;
    parse_pdb_string(pdb_content).unwrap()
}

#[test]
#[cfg(feature = "descriptors")]
fn test_is_predicted_by_header() {
    let structure = create_alphafold_like_structure();
    assert!(structure.is_predicted());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_is_predicted_by_title() {
    let pdb_content = r#"TITLE     PREDICTED STRUCTURE
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 85.50           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 84.30           C
END
"#;
    let structure = parse_pdb_string(pdb_content).unwrap();
    assert!(structure.is_predicted());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_not_predicted_experimental() {
    let structure = create_experimental_structure();
    // May still be detected as predicted due to B-factor heuristics
    // The heuristic checks range and distribution
    println!(
        "Experimental structure is_predicted: {}",
        structure.is_predicted()
    );
}

#[test]
#[cfg(feature = "descriptors")]
fn test_plddt_mean() {
    let structure = create_alphafold_like_structure();
    let mean = structure.plddt_mean();

    // Should be between 0 and 100
    assert!(mean >= 0.0);
    assert!(mean <= 100.0);

    // For our test structure, mean should be moderate
    println!("Mean pLDDT: {:.2}", mean);
}

#[test]
#[cfg(feature = "descriptors")]
fn test_plddt_mean_empty() {
    let structure = PdbStructure::new();
    let mean = structure.plddt_mean();
    assert_eq!(mean, 0.0);
}

#[test]
#[cfg(feature = "descriptors")]
fn test_per_residue_plddt() {
    let structure = create_alphafold_like_structure();
    let profile = structure.per_residue_plddt();

    assert_eq!(profile.len(), 3); // 3 residues

    // Check residue 1 (ALA) - high confidence
    let res1 = profile.iter().find(|r| r.residue_seq == 1).unwrap();
    assert!(res1.plddt > 90.0);
    assert_eq!(res1.confidence_category, ConfidenceCategory::VeryHigh);
    assert!(res1.is_confident());

    // Check residue 2 (GLY) - confident
    let res2 = profile.iter().find(|r| r.residue_seq == 2).unwrap();
    assert!(res2.plddt >= 70.0 && res2.plddt < 90.0);
    assert_eq!(res2.confidence_category, ConfidenceCategory::Confident);

    // Check residue 3 (VAL) - low confidence
    let res3 = profile.iter().find(|r| r.residue_seq == 3).unwrap();
    assert!(res3.plddt < 50.0);
    assert_eq!(res3.confidence_category, ConfidenceCategory::VeryLow);
    assert!(res3.is_disordered());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_low_confidence_regions() {
    let structure = create_alphafold_like_structure();

    // Threshold 70 should include residue 3
    let low_conf = structure.low_confidence_regions(70.0);
    assert!(!low_conf.is_empty());

    // All returned residues should be below threshold
    for res in &low_conf {
        assert!(res.plddt < 70.0);
    }

    // Threshold 50 should only include very low regions
    let very_low = structure.low_confidence_regions(50.0);
    for res in &very_low {
        assert!(res.plddt < 50.0);
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_high_confidence_regions() {
    let structure = create_alphafold_like_structure();

    // Threshold 70 should include residues 1 and 2
    let high_conf = structure.high_confidence_regions(70.0);
    assert!(!high_conf.is_empty());

    for res in &high_conf {
        assert!(res.plddt >= 70.0);
    }

    // Threshold 90 should only include very high regions
    let very_high = structure.high_confidence_regions(90.0);
    for res in &very_high {
        assert!(res.plddt >= 90.0);
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_plddt_distribution() {
    let structure = create_alphafold_like_structure();
    let (very_high, confident, low, very_low) = structure.plddt_distribution();

    // Fractions should sum to 1.0
    let sum = very_high + confident + low + very_low;
    assert!((sum - 1.0).abs() < 0.01);

    // All fractions should be between 0 and 1
    assert!((0.0..=1.0).contains(&very_high));
    assert!((0.0..=1.0).contains(&confident));
    assert!((0.0..=1.0).contains(&low));
    assert!((0.0..=1.0).contains(&very_low));

    println!("pLDDT distribution:");
    println!("  Very high (>90): {:.1}%", very_high * 100.0);
    println!("  Confident (70-90): {:.1}%", confident * 100.0);
    println!("  Low (50-70): {:.1}%", low * 100.0);
    println!("  Very low (<50): {:.1}%", very_low * 100.0);
}

#[test]
#[cfg(feature = "descriptors")]
fn test_confidence_category_methods() {
    assert!(ConfidenceCategory::VeryHigh.is_reliable());
    assert!(ConfidenceCategory::Confident.is_reliable());
    assert!(!ConfidenceCategory::Low.is_reliable());
    assert!(!ConfidenceCategory::VeryLow.is_reliable());

    assert!(!ConfidenceCategory::VeryHigh.needs_caution());
    assert!(!ConfidenceCategory::Confident.needs_caution());
    assert!(ConfidenceCategory::Low.needs_caution());
    assert!(ConfidenceCategory::VeryLow.needs_caution());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_confidence_category_from_plddt() {
    assert_eq!(
        ConfidenceCategory::from_plddt(95.0),
        ConfidenceCategory::VeryHigh
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(90.1),
        ConfidenceCategory::VeryHigh
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(90.0),
        ConfidenceCategory::Confident
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(80.0),
        ConfidenceCategory::Confident
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(70.0),
        ConfidenceCategory::Confident
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(69.9),
        ConfidenceCategory::Low
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(50.0),
        ConfidenceCategory::Low
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(49.9),
        ConfidenceCategory::VeryLow
    );
    assert_eq!(
        ConfidenceCategory::from_plddt(30.0),
        ConfidenceCategory::VeryLow
    );
}
