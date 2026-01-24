//! Integration tests for protein-ligand interaction analysis.

use pdbrust::{PdbStructure, parse_pdb_string};

/// Create a structure with a protein and a ligand
fn create_protein_ligand_structure() -> PdbStructure {
    let pdb_content = r#"
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00 20.00           O
ATOM      5  N   LYS A   2       3.320   1.567   0.000  1.00 20.00           N
ATOM      6  CA  LYS A   2       3.954   2.881   0.000  1.00 20.00           C
ATOM      7  C   LYS A   2       5.464   2.771   0.000  1.00 20.00           C
ATOM      8  O   LYS A   2       6.108   1.726   0.000  1.00 20.00           O
ATOM      9  NZ  LYS A   2       4.500   4.000   0.000  1.00 20.00           N
ATOM     10  N   SER A   3       6.012   3.973   0.000  1.00 20.00           N
ATOM     11  CA  SER A   3       7.445   4.145   0.000  1.00 20.00           C
ATOM     12  C   SER A   3       8.115   2.811   0.000  1.00 20.00           C
ATOM     13  O   SER A   3       7.418   1.800   0.000  1.00 20.00           O
ATOM     14  OG  SER A   3       7.500   5.500   0.000  1.00 20.00           O
HETATM   15  C1  ATP A 100       5.000   5.000   0.000  1.00 20.00           C
HETATM   16  O1  ATP A 100       4.000   5.500   0.500  1.00 20.00           O
HETATM   17  N1  ATP A 100       5.500   5.500  -0.500  1.00 20.00           N
END
"#;
    parse_pdb_string(pdb_content).unwrap()
}

/// Create a structure without ligands
fn create_protein_only_structure() -> PdbStructure {
    let pdb_content = r#"
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00 20.00           O
END
"#;
    parse_pdb_string(pdb_content).unwrap()
}

#[test]
#[cfg(feature = "descriptors")]
fn test_binding_site_basic() {
    let structure = create_protein_ligand_structure();

    let site = structure.binding_site("ATP", 6.0);
    assert!(site.is_some());

    let site = site.unwrap();
    assert_eq!(site.ligand_name, "ATP");
    assert!(!site.contact_residues.is_empty());

    println!("Binding site for ATP:");
    println!("  Contact residues: {}", site.num_residues());
    for res in &site.contact_residues {
        println!(
            "    {}{} {}: {:.2} Å",
            res.chain_id, res.residue_seq, res.residue_name, res.min_distance
        );
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_binding_site_not_found() {
    let structure = create_protein_ligand_structure();

    let site = structure.binding_site("XXX", 5.0);
    assert!(site.is_none());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_binding_site_no_ligands() {
    let structure = create_protein_only_structure();

    let site = structure.binding_site("ATP", 5.0);
    assert!(site.is_none());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_binding_site_distance_cutoff() {
    let structure = create_protein_ligand_structure();

    // Smaller cutoff = fewer contacts
    let site_small = structure.binding_site("ATP", 3.0);
    let site_large = structure.binding_site("ATP", 10.0);

    if let (Some(s1), Some(s2)) = (site_small, site_large) {
        assert!(s1.num_residues() <= s2.num_residues());
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_ligand_interactions_basic() {
    let structure = create_protein_ligand_structure();

    let profile = structure.ligand_interactions("ATP");
    assert!(profile.is_some());

    let profile = profile.unwrap();
    assert_eq!(profile.ligand_name, "ATP");
    assert!(!profile.contact_residues.is_empty());

    println!("Ligand interactions for ATP:");
    println!("  Contact residues: {}", profile.contact_residues.len());
    println!("  H-bonds: {}", profile.hydrogen_bonds.len());
    println!("  Salt bridges: {}", profile.salt_bridges.len());
    println!(
        "  Hydrophobic contacts: {}",
        profile.hydrophobic_contacts.len()
    );
    println!("  Total interactions: {}", profile.total_interactions());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_ligand_interactions_not_found() {
    let structure = create_protein_ligand_structure();

    let profile = structure.ligand_interactions("XXX");
    assert!(profile.is_none());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_all_ligand_interactions() {
    let structure = create_protein_ligand_structure();

    let profiles = structure.all_ligand_interactions();

    // Should find ATP
    assert!(!profiles.is_empty());

    for profile in &profiles {
        println!(
            "Ligand {}: {} contacts, {} total interactions",
            profile.ligand_name,
            profile.contact_residues.len(),
            profile.total_interactions()
        );
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_all_ligand_interactions_no_ligands() {
    let structure = create_protein_only_structure();

    let profiles = structure.all_ligand_interactions();
    assert!(profiles.is_empty());
}

#[test]
#[cfg(feature = "descriptors")]
fn test_interaction_profile_methods() {
    let structure = create_protein_ligand_structure();

    if let Some(profile) = structure.ligand_interactions("ATP") {
        // Test total_interactions
        let total = profile.total_interactions();
        let expected = profile.hydrogen_bonds.len()
            + profile.salt_bridges.len()
            + profile.hydrophobic_contacts.len();
        assert_eq!(total, expected);

        // Test has_interactions
        if total > 0 {
            assert!(profile.has_interactions());
        }
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_binding_site_residues_by_distance() {
    let structure = create_protein_ligand_structure();

    if let Some(site) = structure.binding_site("ATP", 10.0) {
        let sorted = site.residues_by_distance();

        // Verify sorted by distance
        for i in 1..sorted.len() {
            assert!(sorted[i - 1].min_distance <= sorted[i].min_distance);
        }
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_contact_residue_fields() {
    let structure = create_protein_ligand_structure();

    if let Some(site) = structure.binding_site("ATP", 10.0) {
        for res in &site.contact_residues {
            // All fields should be populated
            assert!(!res.chain_id.is_empty());
            assert!(!res.residue_name.is_empty());
            assert!(res.min_distance > 0.0);
            assert!(res.num_contacts > 0);
        }
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_hydrogen_bond_detection() {
    let structure = create_protein_ligand_structure();

    if let Some(profile) = structure.ligand_interactions("ATP") {
        for hb in &profile.hydrogen_bonds {
            // Distance should be within H-bond range
            assert!(hb.distance <= 3.5);

            // Should have protein atom and ligand atom names
            assert!(!hb.protein_atom.is_empty());
            assert!(!hb.ligand_atom.is_empty());
        }
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_salt_bridge_detection() {
    let structure = create_protein_ligand_structure();

    if let Some(profile) = structure.ligand_interactions("ATP") {
        for sb in &profile.salt_bridges {
            // Distance should be within salt bridge range
            assert!(sb.distance <= 4.0);

            // Should have valid residue info
            assert!(!sb.protein_resname.is_empty());
            assert!(!sb.protein_atom.is_empty());
        }
    }
}

#[test]
#[cfg(feature = "descriptors")]
fn test_hydrophobic_contact_detection() {
    let structure = create_protein_ligand_structure();

    if let Some(profile) = structure.ligand_interactions("ATP") {
        for hc in &profile.hydrophobic_contacts {
            // Distance should be within hydrophobic range
            assert!(hc.distance <= 4.0);

            // Should have valid info
            assert!(!hc.protein_resname.is_empty());
        }
    }
}
