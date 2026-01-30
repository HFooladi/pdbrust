//! LDDT (Local Distance Difference Test) demonstration.
//!
//! This example demonstrates how to calculate LDDT scores between
//! protein structures. LDDT is a superposition-free metric widely
//! used in AlphaFold (pLDDT) and CASP evaluations.
//!
//! Run with:
//!   cargo run --example lddt_demo --features geometry

use pdbrust::geometry::{AtomSelection, LddtOptions};
use pdbrust::parse_pdb_file;
use std::path::PathBuf;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== LDDT (Local Distance Difference Test) Demo ===\n");

    // Load the structure
    let path = get_test_file("1UBQ.pdb");
    let reference = parse_pdb_file(&path)?;
    println!(
        "Loaded reference: 1UBQ.pdb ({} atoms)",
        reference.atoms.len()
    );

    // =========================================================================
    // Example 1: Self-comparison (should be perfect)
    // =========================================================================
    println!("\n--- Example 1: Self-LDDT ---");

    let result = reference.lddt_to(&reference)?;
    println!("LDDT Score: {:.4}", result.score);
    println!("Distance pairs evaluated: {}", result.num_pairs);
    println!("Residues evaluated: {}", result.num_residues);
    println!(
        "Per-threshold scores: {:?}",
        result
            .per_threshold_scores
            .iter()
            .map(|s| format!("{:.4}", s))
            .collect::<Vec<_>>()
    );

    // =========================================================================
    // Example 2: Translated structure (LDDT should still be 1.0)
    // =========================================================================
    println!("\n--- Example 2: Translation Invariance ---");

    let mut translated = reference.clone();
    for atom in &mut translated.atoms {
        atom.x += 100.0;
        atom.y += 50.0;
        atom.z += 25.0;
    }

    let result = translated.lddt_to(&reference)?;
    println!("After translation of (100, 50, 25) Angstroms:");
    println!("  LDDT Score: {:.4} (should be 1.0)", result.score);

    // Compare with RMSD (which requires alignment)
    let rmsd = translated.rmsd_to(&reference)?;
    println!(
        "  Direct RMSD: {:.2} Angstroms (large without alignment)",
        rmsd
    );

    // =========================================================================
    // Example 3: Rotated structure (LDDT should still be 1.0)
    // =========================================================================
    println!("\n--- Example 3: Rotation Invariance ---");

    let mut rotated = reference.clone();
    let cos_45 = std::f64::consts::FRAC_PI_4.cos();
    let sin_45 = std::f64::consts::FRAC_PI_4.sin();
    for atom in &mut rotated.atoms {
        let x = atom.x;
        let y = atom.y;
        atom.x = x * cos_45 - y * sin_45;
        atom.y = x * sin_45 + y * cos_45;
    }

    let result = rotated.lddt_to(&reference)?;
    println!("After 45-degree rotation around Z-axis:");
    println!("  LDDT Score: {:.4} (should be 1.0)", result.score);

    // =========================================================================
    // Example 4: Perturbed structure
    // =========================================================================
    println!("\n--- Example 4: Perturbed Structure ---");

    let mut perturbed = reference.clone();
    // Perturb every 5th residue by 2 Angstroms
    for (i, atom) in perturbed.atoms.iter_mut().enumerate() {
        if i % 5 == 0 {
            atom.y += 2.0;
        }
    }

    let result = perturbed.lddt_to(&reference)?;
    println!("After perturbing every 5th residue by 2 Angstroms:");
    println!("  LDDT Score: {:.4}", result.score);
    println!("  Per-threshold scores:");
    for (i, score) in result.per_threshold_scores.iter().enumerate() {
        let threshold = match i {
            0 => 0.5,
            1 => 1.0,
            2 => 2.0,
            3 => 4.0,
            _ => 0.0,
        };
        println!("    {:.1} Angstrom threshold: {:.4}", threshold, score);
    }

    // =========================================================================
    // Example 5: Custom options
    // =========================================================================
    println!("\n--- Example 5: Custom Options ---");

    // Stricter thresholds
    let strict_options = LddtOptions::default().with_thresholds(vec![0.25, 0.5, 1.0]);
    let result_strict =
        perturbed.lddt_to_with_options(&reference, AtomSelection::CaOnly, strict_options)?;

    println!("With stricter thresholds [0.25, 0.5, 1.0]:");
    println!("  LDDT Score: {:.4}", result_strict.score);

    // Smaller inclusion radius
    let small_radius = LddtOptions::default().with_inclusion_radius(8.0);
    let result_small =
        perturbed.lddt_to_with_options(&reference, AtomSelection::CaOnly, small_radius)?;

    println!("With smaller inclusion radius (8.0 Angstroms):");
    println!(
        "  LDDT Score: {:.4} ({} pairs)",
        result_small.score, result_small.num_pairs
    );

    // =========================================================================
    // Example 6: Per-residue LDDT
    // =========================================================================
    println!("\n--- Example 6: Per-Residue LDDT ---");

    let per_res = perturbed.per_residue_lddt_to(&reference)?;

    println!("Per-residue LDDT (showing residues with LDDT < 0.9):");
    let low_lddt: Vec<_> = per_res.iter().filter(|r| r.score < 0.9).collect();

    if low_lddt.is_empty() {
        println!("  All residues have LDDT >= 0.9");
    } else {
        for r in low_lddt.iter().take(10) {
            println!(
                "  {}{} {}: LDDT = {:.3} ({} pairs)",
                r.residue_id.0, r.residue_id.1, r.residue_name, r.score, r.num_pairs
            );
        }
        if low_lddt.len() > 10 {
            println!("  ... and {} more residues", low_lddt.len() - 10);
        }
    }

    // =========================================================================
    // Example 7: LDDT vs RMSD comparison
    // =========================================================================
    println!("\n--- Example 7: LDDT vs RMSD Comparison ---");

    // Create a heavily perturbed structure
    let mut heavily_perturbed = reference.clone();
    for (i, atom) in heavily_perturbed.atoms.iter_mut().enumerate() {
        if i % 2 == 0 {
            atom.z += 5.0;
        }
    }

    let lddt_result = heavily_perturbed.lddt_to(&reference)?;
    let rmsd_direct = heavily_perturbed.rmsd_to(&reference)?;
    let (_, aligned_result) = heavily_perturbed.align_to(&reference)?;

    println!("For heavily perturbed structure:");
    println!("  LDDT Score: {:.4}", lddt_result.score);
    println!("  Direct RMSD: {:.2} Angstroms", rmsd_direct);
    println!("  Aligned RMSD: {:.2} Angstroms", aligned_result.rmsd);
    println!("\nNote: LDDT is superposition-free, while RMSD depends on alignment.");

    // =========================================================================
    // Example 8: Backbone atoms selection
    // =========================================================================
    println!("\n--- Example 8: Different Atom Selections ---");

    let lddt_ca = reference.lddt_to_with_options(
        &reference,
        AtomSelection::CaOnly,
        LddtOptions::default(),
    )?;
    let lddt_bb = reference.lddt_to_with_options(
        &reference,
        AtomSelection::Backbone,
        LddtOptions::default(),
    )?;

    println!("Self-LDDT with different selections:");
    println!(
        "  CA atoms only: {:.4} ({} pairs)",
        lddt_ca.score, lddt_ca.num_pairs
    );
    println!(
        "  Backbone atoms: {:.4} ({} pairs)",
        lddt_bb.score, lddt_bb.num_pairs
    );

    println!("\n=== Demo Complete ===");

    Ok(())
}
