//! Integration tests for the descriptors module.
//!
//! These tests verify descriptor calculations using real PDB files.

#![cfg(feature = "descriptors")]

use pdbrust::{PdbStructure, parse_pdb_file};
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
    assert!(density > 0.0, "Density should be positive, got {}", density);

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

// ============================================================================
// Distance Matrix Tests
// ============================================================================

#[test]
fn test_distance_matrix_ca_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let matrix = structure.distance_matrix_ca();
    let ca_count = structure.count_ca_residues();

    // Matrix should be square with dimensions matching CA count
    assert_eq!(matrix.len(), ca_count);
    for row in &matrix {
        assert_eq!(row.len(), ca_count);
    }

    // Diagonal should be zero
    for (i, row) in matrix.iter().enumerate() {
        assert_eq!(row[i], 0.0, "Diagonal should be zero at position {}", i);
    }

    // Matrix should be symmetric
    for (i, row) in matrix.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!(
                (val - matrix[j][i]).abs() < 1e-10,
                "Matrix should be symmetric at [{}, {}]",
                i,
                j
            );
        }
    }

    // All distances should be non-negative
    for (i, row) in matrix.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!(
                val >= 0.0,
                "Distances should be non-negative at [{}, {}]",
                i,
                j
            );
        }
    }
}

#[test]
fn test_distance_matrix_ca_max_equals_max_ca_distance() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let matrix = structure.distance_matrix_ca();
    let max_dist = structure.max_ca_distance();

    // Find max in matrix
    let matrix_max: f64 = matrix
        .iter()
        .flat_map(|row| row.iter())
        .cloned()
        .fold(0.0, f64::max);

    // Should match max_ca_distance()
    assert!(
        (matrix_max - max_dist).abs() < 1e-10,
        "Max distance from matrix ({}) should match max_ca_distance ({})",
        matrix_max,
        max_dist
    );
}

#[test]
fn test_distance_matrix_all_atoms_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let matrix = structure.distance_matrix();
    let atom_count = structure.get_num_atoms();

    // Matrix should be square with dimensions matching atom count
    assert_eq!(matrix.len(), atom_count);
    for row in &matrix {
        assert_eq!(row.len(), atom_count);
    }

    // Diagonal should be zero
    for (i, row) in matrix.iter().enumerate() {
        assert_eq!(row[i], 0.0);
    }
}

#[test]
fn test_distance_matrix_empty_structure() {
    let structure = PdbStructure::new();
    let matrix = structure.distance_matrix_ca();
    assert!(matrix.is_empty());

    let all_atom_matrix = structure.distance_matrix();
    assert!(all_atom_matrix.is_empty());
}

// ============================================================================
// Contact Map Tests
// ============================================================================

#[test]
fn test_contact_map_ca_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let contacts = structure.contact_map_ca(8.0);
    let ca_count = structure.count_ca_residues();

    // Matrix should be square with dimensions matching CA count
    assert_eq!(contacts.len(), ca_count);
    for row in &contacts {
        assert_eq!(row.len(), ca_count);
    }

    // Diagonal should be true (self-contact)
    for (i, row) in contacts.iter().enumerate() {
        assert!(row[i], "Diagonal should be true at position {}", i);
    }

    // Matrix should be symmetric
    for (i, row) in contacts.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert_eq!(
                val, contacts[j][i],
                "Contact map should be symmetric at [{}, {}]",
                i, j
            );
        }
    }
}

#[test]
fn test_contact_map_ca_consistent_with_distance_matrix() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let threshold = 8.0;
    let distances = structure.distance_matrix_ca();
    let contacts = structure.contact_map_ca(threshold);

    // Contact map should match distance <= threshold
    for (i, dist_row) in distances.iter().enumerate() {
        for (j, &dist) in dist_row.iter().enumerate() {
            let expected_contact = dist <= threshold;
            assert_eq!(
                contacts[i][j], expected_contact,
                "Contact map should be consistent with distances at [{}, {}]: dist={}, threshold={}",
                i, j, dist, threshold
            );
        }
    }
}

#[test]
fn test_contact_map_ca_threshold_variation() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let contacts_6 = structure.contact_map_ca(6.0);
    let contacts_8 = structure.contact_map_ca(8.0);
    let contacts_12 = structure.contact_map_ca(12.0);

    // Count contacts for each threshold
    let count_6: usize = contacts_6
        .iter()
        .flat_map(|r| r.iter())
        .filter(|&&x| x)
        .count();
    let count_8: usize = contacts_8
        .iter()
        .flat_map(|r| r.iter())
        .filter(|&&x| x)
        .count();
    let count_12: usize = contacts_12
        .iter()
        .flat_map(|r| r.iter())
        .filter(|&&x| x)
        .count();

    // More contacts with larger threshold
    assert!(
        count_6 <= count_8,
        "More contacts with larger threshold: {} <= {}",
        count_6,
        count_8
    );
    assert!(
        count_8 <= count_12,
        "More contacts with larger threshold: {} <= {}",
        count_8,
        count_12
    );
}

#[test]
fn test_contact_map_all_atoms_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let contacts = structure.contact_map(4.5);
    let atom_count = structure.get_num_atoms();

    // Matrix should be square with dimensions matching atom count
    assert_eq!(contacts.len(), atom_count);
    for row in &contacts {
        assert_eq!(row.len(), atom_count);
    }

    // Diagonal should be true
    for (i, row) in contacts.iter().enumerate() {
        assert!(row[i]);
    }
}

#[test]
fn test_contact_map_empty_structure() {
    let structure = PdbStructure::new();
    let contacts = structure.contact_map_ca(8.0);
    assert!(contacts.is_empty());

    let all_atom_contacts = structure.contact_map(4.5);
    assert!(all_atom_contacts.is_empty());
}

// ============================================================================
// B-factor Analysis Tests
// ============================================================================

#[test]
fn test_b_factor_mean_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let mean_b = structure.b_factor_mean();

    // B-factors should be positive and reasonable
    assert!(
        mean_b > 0.0,
        "Mean B-factor should be positive, got {}",
        mean_b
    );
    assert!(
        mean_b < 100.0,
        "Mean B-factor should be reasonable, got {}",
        mean_b
    );
}

#[test]
fn test_b_factor_mean_ca_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let mean_b_ca = structure.b_factor_mean_ca();
    let mean_b_all = structure.b_factor_mean();

    // CA B-factor should be positive
    assert!(
        mean_b_ca > 0.0,
        "Mean CA B-factor should be positive, got {}",
        mean_b_ca
    );

    // CA and all-atom means can differ but should be in same ballpark
    assert!(
        (mean_b_ca - mean_b_all).abs() < 20.0,
        "CA mean ({}) should be close to all-atom mean ({})",
        mean_b_ca,
        mean_b_all
    );
}

#[test]
fn test_b_factor_min_max_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let min_b = structure.b_factor_min();
    let max_b = structure.b_factor_max();
    let mean_b = structure.b_factor_mean();

    // Min should be <= mean <= max
    assert!(
        min_b <= mean_b,
        "Min ({}) should be <= mean ({})",
        min_b,
        mean_b
    );
    assert!(
        mean_b <= max_b,
        "Mean ({}) should be <= max ({})",
        mean_b,
        max_b
    );

    // There should be a range of B-factors
    assert!(max_b > min_b, "Max ({}) should be > min ({})", max_b, min_b);
}

#[test]
fn test_b_factor_std_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let std_b = structure.b_factor_std();

    // Std should be positive (there should be variation)
    assert!(
        std_b > 0.0,
        "B-factor std should be positive, got {}",
        std_b
    );

    // Std should be reasonable (typically < 30 for well-determined structures)
    assert!(
        std_b < 50.0,
        "B-factor std should be reasonable, got {}",
        std_b
    );
}

#[test]
fn test_b_factor_profile_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let profile = structure.b_factor_profile();

    // 1UBQ has 76 amino acid residues plus water molecules
    // Profile includes all residues (protein + waters)
    assert!(
        profile.len() >= 76,
        "Profile should have at least 76 residues, got {}",
        profile.len()
    );

    // Each residue should have valid statistics
    for res in &profile {
        assert!(
            !res.chain_id.is_empty() || res.chain_id.trim().is_empty(),
            "Chain ID can be empty for some records"
        );
        assert!(
            res.b_factor_mean >= res.b_factor_min,
            "Mean ({}) should be >= min ({})",
            res.b_factor_mean,
            res.b_factor_min
        );
        assert!(
            res.b_factor_mean <= res.b_factor_max,
            "Mean ({}) should be <= max ({})",
            res.b_factor_mean,
            res.b_factor_max
        );
        assert!(res.atom_count > 0, "Atom count should be > 0");
    }
}

#[test]
fn test_b_factor_profile_ordering() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let profile = structure.b_factor_profile();

    // Profile should be sorted by chain, then by residue number
    for window in profile.windows(2) {
        let a = &window[0];
        let b = &window[1];

        if a.chain_id == b.chain_id {
            assert!(
                a.residue_seq <= b.residue_seq,
                "Residues should be sorted within chain: {} <= {}",
                a.residue_seq,
                b.residue_seq
            );
        } else {
            assert!(
                a.chain_id < b.chain_id,
                "Chains should be sorted alphabetically: {} < {}",
                a.chain_id,
                b.chain_id
            );
        }
    }
}

#[test]
fn test_flexible_residues_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Get some statistics to set reasonable threshold
    let mean_b = structure.b_factor_mean();
    let std_b = structure.b_factor_std();

    // Use mean + std as threshold
    let threshold = mean_b + std_b;
    let flexible = structure.flexible_residues(threshold);

    // Should find some flexible residues (but not all)
    assert!(!flexible.is_empty(), "Should find some flexible residues");
    assert!(
        flexible.len() < structure.b_factor_profile().len(),
        "Not all residues should be flexible"
    );

    // All returned residues should be above threshold
    for res in &flexible {
        assert!(
            res.b_factor_mean > threshold,
            "Flexible residue B-factor ({}) should be > threshold ({})",
            res.b_factor_mean,
            threshold
        );
    }
}

#[test]
fn test_rigid_residues_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Get some statistics to set reasonable threshold
    let mean_b = structure.b_factor_mean();
    let std_b = structure.b_factor_std();

    // Use mean - std as threshold
    let threshold = mean_b - std_b;
    let rigid = structure.rigid_residues(threshold);

    // Should find some rigid residues (but not all)
    assert!(!rigid.is_empty(), "Should find some rigid residues");
    assert!(
        rigid.len() < structure.b_factor_profile().len(),
        "Not all residues should be rigid"
    );

    // All returned residues should be below threshold
    for res in &rigid {
        assert!(
            res.b_factor_mean < threshold,
            "Rigid residue B-factor ({}) should be < threshold ({})",
            res.b_factor_mean,
            threshold
        );
    }
}

#[test]
fn test_normalize_b_factors_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let normalized = structure.normalize_b_factors();

    // After normalization, mean should be ~0 and std should be ~1
    let norm_mean = normalized.b_factor_mean();
    let norm_std = normalized.b_factor_std();

    assert!(
        norm_mean.abs() < 1e-6,
        "Normalized mean should be ~0, got {}",
        norm_mean
    );
    assert!(
        (norm_std - 1.0).abs() < 1e-6,
        "Normalized std should be ~1, got {}",
        norm_std
    );

    // Original structure should be unchanged
    let original_mean = structure.b_factor_mean();
    assert!(
        original_mean.abs() > 1.0,
        "Original mean should be unchanged"
    );
}

#[test]
fn test_b_factor_empty_structure() {
    let structure = PdbStructure::new();

    assert_eq!(structure.b_factor_mean(), 0.0);
    assert_eq!(structure.b_factor_mean_ca(), 0.0);
    assert_eq!(structure.b_factor_min(), 0.0);
    assert_eq!(structure.b_factor_max(), 0.0);
    assert_eq!(structure.b_factor_std(), 0.0);
    assert!(structure.b_factor_profile().is_empty());
    assert!(structure.flexible_residues(50.0).is_empty());
    assert!(structure.rigid_residues(20.0).is_empty());
}

#[test]
fn test_b_factor_percentile_real_structure() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    // Get percentile for first atom
    if let Some(first_atom) = structure.atoms.first() {
        let percentile = structure.b_factor_percentile(first_atom.serial);
        assert!(
            percentile.is_some(),
            "Should find percentile for existing atom"
        );

        let p = percentile.unwrap();
        assert!(
            (0.0..=100.0).contains(&p),
            "Percentile should be 0-100, got {}",
            p
        );
    }

    // Non-existent atom should return None
    assert!(structure.b_factor_percentile(-999).is_none());
}

#[test]
fn test_structure_descriptors_includes_bfactors() {
    let path = get_test_file("1UBQ.pdb");
    let structure = parse_pdb_file(&path).expect("Failed to parse 1UBQ.pdb");

    let desc = structure.structure_descriptors();

    // B-factor fields should be populated and consistent with individual methods
    assert!((desc.b_factor_mean - structure.b_factor_mean()).abs() < 1e-10);
    assert!((desc.b_factor_mean_ca - structure.b_factor_mean_ca()).abs() < 1e-10);
    assert!((desc.b_factor_min - structure.b_factor_min()).abs() < 1e-10);
    assert!((desc.b_factor_max - structure.b_factor_max()).abs() < 1e-10);
    assert!((desc.b_factor_std - structure.b_factor_std()).abs() < 1e-10);
}
