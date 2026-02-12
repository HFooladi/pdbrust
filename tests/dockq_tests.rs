#![cfg(feature = "dockq")]

use pdbrust::PdbStructure;
use pdbrust::dockq::{ChainMappingStrategy, DockQOptions, DockQQuality, find_chain_mapping};
use pdbrust::records::Atom;

// ============================================================================
// Helper Functions
// ============================================================================

#[allow(clippy::too_many_arguments)]
fn create_atom(
    serial: i32,
    name: &str,
    residue_name: &str,
    chain_id: &str,
    residue_seq: i32,
    x: f64,
    y: f64,
    z: f64,
    element: &str,
) -> Atom {
    Atom {
        serial,
        name: name.to_string(),
        alt_loc: None,
        residue_name: residue_name.to_string(),
        chain_id: chain_id.to_string(),
        residue_seq,
        ins_code: None,
        is_hetatm: false,
        x,
        y,
        z,
        occupancy: 1.0,
        temp_factor: 20.0,
        element: element.to_string(),
    }
}

/// Create a two-chain complex with an interface.
/// Chain A: 5 residues along x-axis (y=0)
/// Chain B: 5 residues along x-axis (y=4.0, within contact distance)
fn create_dimer() -> PdbStructure {
    let residue_names_a = ["ALA", "GLY", "VAL", "LEU", "ILE"];
    let residue_names_b = ["PHE", "TRP", "TYR", "SER", "THR"];

    let mut structure = PdbStructure::new();
    let mut serial = 1;

    for (i, resname) in residue_names_a.iter().enumerate() {
        let x = i as f64 * 3.8;
        let seq = (i + 1) as i32;
        structure.atoms.push(create_atom(
            serial,
            "N",
            resname,
            "A",
            seq,
            x - 0.5,
            -0.5,
            0.0,
            "N",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial, "CA", resname, "A", seq, x, 0.0, 0.0, "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "C",
            resname,
            "A",
            seq,
            x + 0.5,
            0.0,
            0.0,
            "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "O",
            resname,
            "A",
            seq,
            x + 0.5,
            0.5,
            0.0,
            "O",
        ));
        serial += 1;
    }

    for (i, resname) in residue_names_b.iter().enumerate() {
        let x = i as f64 * 3.8;
        let seq = (i + 1) as i32;
        structure.atoms.push(create_atom(
            serial,
            "N",
            resname,
            "B",
            seq,
            x - 0.5,
            3.5,
            0.0,
            "N",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial, "CA", resname, "B", seq, x, 4.0, 0.0, "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "C",
            resname,
            "B",
            seq,
            x + 0.5,
            4.0,
            0.0,
            "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "O",
            resname,
            "B",
            seq,
            x + 0.5,
            4.5,
            0.0,
            "O",
        ));
        serial += 1;
    }

    structure
}

/// Create a three-chain complex with multiple interfaces.
fn create_trimer() -> PdbStructure {
    let mut structure = create_dimer();
    let mut serial = structure.atoms.len() as i32 + 1;

    let residue_names_c = ["ASP", "GLU", "LYS", "ARG", "HIS"];

    for (i, resname) in residue_names_c.iter().enumerate() {
        let x = i as f64 * 3.8;
        let seq = (i + 1) as i32;
        // Chain C offset in z direction, close to chain A
        structure.atoms.push(create_atom(
            serial,
            "N",
            resname,
            "C",
            seq,
            x - 0.5,
            0.0,
            3.5,
            "N",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial, "CA", resname, "C", seq, x, 0.0, 4.0, "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "C",
            resname,
            "C",
            seq,
            x + 0.5,
            0.0,
            4.0,
            "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "O",
            resname,
            "C",
            seq,
            x + 0.5,
            0.5,
            4.0,
            "O",
        ));
        serial += 1;
    }

    structure
}

// ============================================================================
// Self-Comparison Tests
// ============================================================================

#[test]
fn test_dockq_self_comparison_perfect_score() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    assert!(
        (result.total_dockq - 1.0).abs() < 1e-6,
        "Self-comparison should give DockQ = 1.0, got {}",
        result.total_dockq
    );
    assert_eq!(result.num_interfaces, 1);
}

#[test]
fn test_dockq_self_comparison_fnat_one() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    for iface in &result.interfaces {
        assert!(
            (iface.fnat - 1.0).abs() < 1e-6,
            "Self-comparison fnat should be 1.0, got {}",
            iface.fnat
        );
    }
}

#[test]
fn test_dockq_self_comparison_irmsd_zero() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    for iface in &result.interfaces {
        assert!(
            iface.irmsd < 1e-6,
            "Self-comparison iRMSD should be ~0, got {}",
            iface.irmsd
        );
    }
}

#[test]
fn test_dockq_self_comparison_lrmsd_zero() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    for iface in &result.interfaces {
        assert!(
            iface.lrmsd < 1e-6,
            "Self-comparison LRMSD should be ~0, got {}",
            iface.lrmsd
        );
    }
}

#[test]
fn test_dockq_self_comparison_high_quality() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    for iface in &result.interfaces {
        assert_eq!(iface.quality, DockQQuality::High);
    }
}

// ============================================================================
// Translated Model Tests
// ============================================================================

#[test]
fn test_dockq_rigid_translation_preserves_contacts() {
    let native = create_dimer();
    let mut model = native.clone();

    // Translate the entire model rigidly
    for atom in &mut model.atoms {
        atom.x += 5.0;
        atom.y += 3.0;
        atom.z += 1.0;
    }

    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };

    let result = model.dockq_to_with_options(&native, options).unwrap();

    // Rigid translation preserves all internal distances → fnat = 1.0
    for iface in &result.interfaces {
        assert!(
            (iface.fnat - 1.0).abs() < 1e-6,
            "Rigid translation should preserve contacts, got fnat={}",
            iface.fnat
        );
    }
}

// ============================================================================
// Chain Mapping Tests
// ============================================================================

#[test]
fn test_chain_mapping_auto_identical() {
    let structure = create_dimer();
    let mapping = find_chain_mapping(&structure, &structure).unwrap();

    assert_eq!(mapping.len(), 2);
    assert!(mapping.contains(&("A".to_string(), "A".to_string())));
    assert!(mapping.contains(&("B".to_string(), "B".to_string())));
}

#[test]
fn test_chain_mapping_auto_swapped() {
    let native = create_dimer();

    // Create model with swapped chain labels
    let mut model = PdbStructure::new();
    for atom in &native.atoms {
        let mut new_atom = atom.clone();
        new_atom.chain_id = if atom.chain_id == "A" {
            "B".to_string()
        } else {
            "A".to_string()
        };
        model.atoms.push(new_atom);
    }

    let mapping = find_chain_mapping(&model, &native).unwrap();
    assert_eq!(mapping.len(), 2);

    // Model B should map to Native A (same sequence) and vice versa
    assert!(mapping.contains(&("B".to_string(), "A".to_string())));
    assert!(mapping.contains(&("A".to_string(), "B".to_string())));
}

#[test]
fn test_dockq_with_swapped_chains() {
    let native = create_dimer();

    // Create model with swapped chain labels
    let mut model = PdbStructure::new();
    for atom in &native.atoms {
        let mut new_atom = atom.clone();
        new_atom.chain_id = if atom.chain_id == "A" {
            "B".to_string()
        } else {
            "A".to_string()
        };
        model.atoms.push(new_atom);
    }

    // Auto chain mapping should handle this
    let result = model.dockq_to(&native).unwrap();
    assert!(
        (result.total_dockq - 1.0).abs() < 1e-6,
        "Swapped chains should still give DockQ = 1.0 with auto-mapping, got {}",
        result.total_dockq
    );
}

#[test]
fn test_explicit_chain_mapping() {
    let structure = create_dimer();
    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };

    let result = structure
        .dockq_to_with_options(&structure, options)
        .unwrap();
    assert!((result.total_dockq - 1.0).abs() < 1e-6);
}

// ============================================================================
// Quality Classification Tests
// ============================================================================

#[test]
fn test_quality_classification_boundaries() {
    assert_eq!(DockQQuality::from_score(0.0), DockQQuality::Incorrect);
    assert_eq!(DockQQuality::from_score(0.22), DockQQuality::Incorrect);
    assert_eq!(DockQQuality::from_score(0.23), DockQQuality::Acceptable);
    assert_eq!(DockQQuality::from_score(0.48), DockQQuality::Acceptable);
    assert_eq!(DockQQuality::from_score(0.49), DockQQuality::Medium);
    assert_eq!(DockQQuality::from_score(0.79), DockQQuality::Medium);
    assert_eq!(DockQQuality::from_score(0.80), DockQQuality::High);
    assert_eq!(DockQQuality::from_score(1.0), DockQQuality::High);
}

#[test]
fn test_quality_display() {
    assert_eq!(format!("{}", DockQQuality::Incorrect), "Incorrect");
    assert_eq!(format!("{}", DockQQuality::Acceptable), "Acceptable");
    assert_eq!(format!("{}", DockQQuality::Medium), "Medium");
    assert_eq!(format!("{}", DockQQuality::High), "High");
}

// ============================================================================
// Multi-Interface Tests
// ============================================================================

#[test]
fn test_trimer_multiple_interfaces() {
    let structure = create_trimer();
    let result = structure.dockq_to(&structure).unwrap();

    // A trimer should have up to 3 interfaces: A-B, A-C, B-C
    // At least A-B should have contacts (from the dimer base)
    assert!(
        result.num_interfaces >= 1,
        "Trimer should have at least 1 interface, got {}",
        result.num_interfaces
    );

    // Self-comparison should still be perfect
    assert!(
        (result.total_dockq - 1.0).abs() < 1e-6,
        "Trimer self-comparison should give DockQ = 1.0, got {}",
        result.total_dockq
    );
}

// ============================================================================
// DockQ Options Tests
// ============================================================================

#[test]
fn test_default_options() {
    let opts = DockQOptions::default();
    assert!((opts.contact_threshold - 5.0).abs() < 1e-10);
    assert!((opts.interface_threshold - 10.0).abs() < 1e-10);
}

#[test]
fn test_custom_contact_threshold() {
    let structure = create_dimer();

    // Very small threshold should find fewer contacts
    let options_tight = DockQOptions {
        contact_threshold: 2.0,
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };

    // Very large threshold should find more contacts
    let options_wide = DockQOptions {
        contact_threshold: 10.0,
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };

    // Tight threshold might fail with no contacts, which is expected
    let result_tight = structure.dockq_to_with_options(&structure, options_tight);
    let result_wide = structure
        .dockq_to_with_options(&structure, options_wide)
        .unwrap();

    // Wide threshold should find contacts and give perfect score
    assert!((result_wide.total_dockq - 1.0).abs() < 1e-6);

    // If tight threshold succeeds, it should have fewer contacts
    if let Ok(tight) = result_tight {
        for (iface_t, iface_w) in tight.interfaces.iter().zip(result_wide.interfaces.iter()) {
            assert!(iface_t.num_native_contacts <= iface_w.num_native_contacts);
        }
    }
}

// ============================================================================
// Perturbed Interface Tests
// ============================================================================

#[test]
fn test_dockq_one_chain_moved() {
    let native = create_dimer();
    let mut model = native.clone();

    // Move chain B away from chain A in the model
    for atom in &mut model.atoms {
        if atom.chain_id == "B" {
            atom.y += 20.0; // Move 20 A away
        }
    }

    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };

    let result = model.dockq_to_with_options(&native, options).unwrap();

    // Moving one chain away should degrade fnat (no contacts in model)
    for iface in &result.interfaces {
        assert!(
            iface.fnat < 0.5,
            "Moving chain should reduce fnat, got {}",
            iface.fnat
        );
        assert!(
            iface.dockq < 0.5,
            "Moving chain should reduce DockQ, got {}",
            iface.dockq
        );
    }
}

// ============================================================================
// Edge Case Tests
// ============================================================================

#[test]
fn test_dockq_no_interface_contacts() {
    // Create two chains very far apart
    let mut structure = PdbStructure::new();
    let mut serial = 1;
    for i in 0..4_i32 {
        let x = i as f64 * 3.8;
        let seq = i + 1;
        structure
            .atoms
            .push(create_atom(serial, "N", "ALA", "A", seq, x, 0.0, 0.0, "N"));
        serial += 1;
        structure
            .atoms
            .push(create_atom(serial, "CA", "ALA", "A", seq, x, 0.0, 0.0, "C"));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "C",
            "ALA",
            "A",
            seq,
            x + 0.5,
            0.0,
            0.0,
            "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "O",
            "ALA",
            "A",
            seq,
            x + 0.5,
            0.5,
            0.0,
            "O",
        ));
        serial += 1;
    }
    for i in 0..4_i32 {
        let x = i as f64 * 3.8;
        let seq = i + 1;
        // Chain B is 100 A away
        structure.atoms.push(create_atom(
            serial, "N", "GLY", "B", seq, x, 100.0, 0.0, "N",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial, "CA", "GLY", "B", seq, x, 100.0, 0.0, "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "C",
            "GLY",
            "B",
            seq,
            x + 0.5,
            100.0,
            0.0,
            "C",
        ));
        serial += 1;
        structure.atoms.push(create_atom(
            serial,
            "O",
            "GLY",
            "B",
            seq,
            x + 0.5,
            100.5,
            0.0,
            "O",
        ));
        serial += 1;
    }

    let result = structure.dockq_to(&structure);
    assert!(
        result.is_err(),
        "Should fail when chains have no interface contacts"
    );
}

// ============================================================================
// Sequence Alignment Tests
// ============================================================================

#[test]
fn test_sequence_alignment_identity() {
    use pdbrust::dockq::sequence_align::{AlignmentParams, align_sequences};

    let seq: Vec<String> = ["ALA", "GLY", "VAL"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let result = align_sequences(&seq, &seq, &AlignmentParams::default());

    assert!((result.identity - 1.0).abs() < 1e-10);
    assert_eq!(result.num_aligned, 3);
}

#[test]
fn test_sequence_identity_function() {
    use pdbrust::dockq::sequence_identity;

    let seq1: Vec<String> = ["ALA", "GLY", "VAL"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let seq2: Vec<String> = ["ALA", "GLY", "VAL"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let seq3: Vec<String> = ["LEU", "ILE", "PHE"]
        .iter()
        .map(|s| s.to_string())
        .collect();

    assert!((sequence_identity(&seq1, &seq2) - 1.0).abs() < 1e-10);
    assert!((sequence_identity(&seq1, &seq3) - 0.0).abs() < 1e-10);
}

// ============================================================================
// Interface Result Validation Tests
// ============================================================================

#[test]
fn test_interface_result_fields() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    for iface in &result.interfaces {
        // All scores should be in valid ranges
        assert!(iface.fnat >= 0.0 && iface.fnat <= 1.0);
        assert!(iface.fnonnat >= 0.0 && iface.fnonnat <= 1.0);
        assert!(iface.f1 >= 0.0 && iface.f1 <= 1.0);
        assert!(iface.irmsd >= 0.0);
        assert!(iface.lrmsd >= 0.0);
        assert!(iface.dockq >= 0.0 && iface.dockq <= 1.0);
        assert!(iface.num_native_contacts > 0);
    }
}

#[test]
fn test_dockq_result_chain_mapping_populated() {
    let structure = create_dimer();
    let result = structure.dockq_to(&structure).unwrap();

    assert!(!result.chain_mapping.is_empty());
    assert_eq!(result.chain_mapping.len(), 2);
}
