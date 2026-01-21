//! B-factor (temperature factor) analysis example
//!
//! This example demonstrates:
//! - Computing B-factor statistics (mean, std, min, max)
//! - Generating per-residue B-factor profiles
//! - Identifying flexible and rigid regions
//! - Normalizing B-factors for cross-structure comparison
//!
//! Run with: cargo run --example b_factor_demo --features descriptors

use pdbrust::parse_pdb_file;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("================================================");
    println!("PDBRust B-factor Analysis Demo");
    println!("================================================\n");

    // Load structure
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join("1UBQ.pdb");

    println!("Loading structure from {:?}", pdb_path);
    let structure = parse_pdb_file(&pdb_path)?;
    println!("Total atoms: {}", structure.atoms.len());
    println!();

    // ====================
    // B-factor Statistics
    // ====================
    println!("=== B-factor Statistics ===");

    let mean_b = structure.b_factor_mean();
    let mean_b_ca = structure.b_factor_mean_ca();
    let min_b = structure.b_factor_min();
    let max_b = structure.b_factor_max();
    let std_b = structure.b_factor_std();

    println!("Mean B-factor (all atoms): {:.2} Å²", mean_b);
    println!("Mean B-factor (CA only):   {:.2} Å²", mean_b_ca);
    println!("Min B-factor:              {:.2} Å²", min_b);
    println!("Max B-factor:              {:.2} Å²", max_b);
    println!("Std deviation:             {:.2} Å²", std_b);
    println!();

    // ====================
    // Per-Residue Profile
    // ====================
    println!("=== Per-Residue B-factor Profile ===");

    let profile = structure.b_factor_profile();
    println!("Residues analyzed: {}", profile.len());

    println!("\nFirst 10 residues:");
    println!(
        "{:<6} {:<8} {:<5} {:<10} {:<10} {:<10}",
        "Chain", "ResSeq", "Name", "Mean", "Min", "Max"
    );
    println!("{}", "-".repeat(55));

    for res in profile.iter().take(10) {
        println!(
            "{:<6} {:<8} {:<5} {:<10.2} {:<10.2} {:<10.2}",
            res.chain_id,
            res.residue_seq,
            res.residue_name,
            res.b_factor_mean,
            res.b_factor_min,
            res.b_factor_max
        );
    }
    println!();

    // ====================
    // Flexible Regions
    // ====================
    println!("=== Flexible Regions (High B-factor) ===");

    let threshold_high = 30.0;
    let flexible = structure.flexible_residues(threshold_high);
    println!(
        "Residues with mean B-factor > {:.1} Å²: {}",
        threshold_high,
        flexible.len()
    );

    if !flexible.is_empty() {
        println!("\nTop 5 most flexible residues:");
        for res in flexible.iter().take(5) {
            println!(
                "  {}{} {}: mean={:.2} Å², max={:.2} Å²",
                res.chain_id,
                res.residue_seq,
                res.residue_name,
                res.b_factor_mean,
                res.b_factor_max
            );
        }
    }
    println!();

    // ====================
    // Rigid Regions
    // ====================
    println!("=== Rigid Regions (Low B-factor) ===");

    let threshold_low = 15.0;
    let rigid = structure.rigid_residues(threshold_low);
    println!(
        "Residues with mean B-factor < {:.1} Å²: {}",
        threshold_low,
        rigid.len()
    );

    if !rigid.is_empty() {
        println!("\nTop 5 most rigid residues:");
        for res in rigid.iter().take(5) {
            println!(
                "  {}{} {}: mean={:.2} Å², min={:.2} Å²",
                res.chain_id,
                res.residue_seq,
                res.residue_name,
                res.b_factor_mean,
                res.b_factor_min
            );
        }
    }
    println!();

    // ====================
    // B-factor Normalization
    // ====================
    println!("=== B-factor Normalization ===");

    let normalized = structure.normalize_b_factors();
    let norm_mean = normalized.b_factor_mean();
    let norm_std = normalized.b_factor_std();

    println!("After Z-score normalization:");
    println!("  Mean: {:.4} (should be ~0)", norm_mean);
    println!("  Std:  {:.4} (should be ~1)", norm_std);
    println!();

    println!("Normalization is useful for:");
    println!("  - Comparing B-factors across different structures");
    println!("  - Identifying relative flexibility within a structure");
    println!("  - Detecting unusually mobile or rigid regions");
    println!();

    // ====================
    // Percentile Analysis
    // ====================
    println!("=== Percentile Analysis ===");

    // Get a few sample atoms and their percentiles
    if structure.atoms.len() >= 3 {
        let atom1 = &structure.atoms[0];
        let atom2 = &structure.atoms[structure.atoms.len() / 2];
        let atom3 = &structure.atoms[structure.atoms.len() - 1];

        if let Some(p1) = structure.b_factor_percentile(atom1.serial) {
            println!(
                "Atom {} ({}): B={:.2} Å², percentile={:.1}%",
                atom1.serial,
                atom1.name.trim(),
                atom1.temp_factor,
                p1 * 100.0
            );
        }
        if let Some(p2) = structure.b_factor_percentile(atom2.serial) {
            println!(
                "Atom {} ({}): B={:.2} Å², percentile={:.1}%",
                atom2.serial,
                atom2.name.trim(),
                atom2.temp_factor,
                p2 * 100.0
            );
        }
        if let Some(p3) = structure.b_factor_percentile(atom3.serial) {
            println!(
                "Atom {} ({}): B={:.2} Å², percentile={:.1}%",
                atom3.serial,
                atom3.name.trim(),
                atom3.temp_factor,
                p3 * 100.0
            );
        }
    }
    println!();

    // ====================
    // Summary
    // ====================
    println!("=== Summary ===");
    println!(
        "
B-factor analysis for {}:
  - Total atoms: {}
  - Mean B-factor: {:.2} Å²
  - Range: {:.2} - {:.2} Å²
  - Flexible residues (B > {:.0} Å²): {}
  - Rigid residues (B < {:.0} Å²): {}
",
        pdb_path.file_name().unwrap().to_string_lossy(),
        structure.atoms.len(),
        mean_b,
        min_b,
        max_b,
        threshold_high,
        flexible.len(),
        threshold_low,
        rigid.len()
    );

    println!("================================================");
    println!("B-factor analysis demo completed successfully!");
    println!("================================================");

    Ok(())
}
