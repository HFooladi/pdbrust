//! Ligand Pose Quality Demo
//!
//! Demonstrates PoseBusters-style geometry checks for protein-ligand complexes.
//!
//! Run with:
//! ```sh
//! cargo run --example ligand_quality_demo --features "ligand-quality"
//! ```
//!
//! Or with a specific PDB file:
//! ```sh
//! cargo run --example ligand_quality_demo --features "ligand-quality" -- path/to/complex.pdb
//! ```

use pdbrust::records::Atom;
use pdbrust::{PdbStructure, parse_pdb_file};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Check if a file was provided as argument
    let args: Vec<String> = std::env::args().collect();

    let structure = if args.len() > 1 {
        // Load user-provided PDB file
        let path = &args[1];
        println!("Loading structure from: {}", path);
        parse_pdb_file(path)?
    } else {
        // Create a demo structure
        println!("Creating demo protein-ligand complex...\n");
        create_demo_complex()
    };

    // Get list of ligands
    let ligand_names = structure.get_ligand_names();
    println!(
        "Structure contains {} ligand(s): {:?}\n",
        ligand_names.len(),
        ligand_names
    );

    if ligand_names.is_empty() {
        println!("No ligands found in structure.");
        return Ok(());
    }

    // Analyze each ligand
    println!("=== Ligand Pose Quality Analysis ===\n");

    for ligand_name in &ligand_names {
        if let Some(report) = structure.ligand_pose_quality(ligand_name) {
            // Print header
            let status = if report.is_geometry_valid {
                "✓ PASS"
            } else {
                "✗ FAIL"
            };
            println!(
                "Ligand: {} ({}{}) - {}",
                report.ligand_name, report.ligand_chain_id, report.ligand_residue_seq, status
            );
            println!("  Atoms: {}", report.ligand_atom_count);
            println!();

            // Distance check
            println!("  Distance Check:");
            if report.min_protein_ligand_distance.is_finite() {
                println!(
                    "    Min protein-ligand distance: {:.2} Å",
                    report.min_protein_ligand_distance
                );
            } else {
                println!("    Min protein-ligand distance: N/A (no protein atoms)");
            }
            println!("    Protein clashes: {}", report.num_clashes);
            if !report.clashes.is_empty() {
                println!(
                    "    Worst clash severity: {:.2}x",
                    report.worst_clash_severity
                );
                println!("    Top clashes:");
                for (i, clash) in report.clashes.iter().take(3).enumerate() {
                    println!(
                        "      {}. {} {} {} {} - {} {}: {:.2}Å (expected >{:.2}Å)",
                        i + 1,
                        clash.protein_chain_id,
                        clash.protein_residue_name,
                        clash.protein_residue_seq,
                        clash.protein_atom_name.trim(),
                        clash.ligand_atom_name.trim(),
                        clash.ligand_element,
                        clash.distance,
                        clash.expected_min_distance
                    );
                }
            }
            let dist_status = if report.passes_distance_check {
                "PASS"
            } else {
                "FAIL"
            };
            println!("    Status: {}", dist_status);
            println!();

            // Volume overlap check
            println!("  Volume Overlap Check:");
            println!(
                "    Protein volume overlap: {:.1}%",
                report.protein_volume_overlap_pct
            );
            println!("    Threshold: <7.5%");
            let overlap_status = if report.passes_overlap_check {
                "PASS"
            } else {
                "FAIL"
            };
            println!("    Status: {}", overlap_status);
            println!();

            // Cofactor clashes
            if !report.cofactor_clashes.is_empty() {
                println!("  Cofactor Clashes: {}", report.num_cofactor_clashes);
            }

            // Overall result
            println!("  Overall: {}", status);
            println!();
            println!("{}", "-".repeat(50));
            println!();
        }
    }

    // Summary
    let reports = structure.all_ligand_pose_quality();
    let passed = reports.iter().filter(|r| r.is_geometry_valid).count();
    let failed = reports.len() - passed;

    println!("\n=== Summary ===");
    println!("Total ligands: {}", reports.len());
    println!("Passed: {}", passed);
    println!("Failed: {}", failed);

    if failed > 0 {
        println!("\nFailed ligands:");
        for report in reports.iter().filter(|r| !r.is_geometry_valid) {
            println!(
                "  - {} ({}{}) - {} clashes, {:.1}% overlap",
                report.ligand_name,
                report.ligand_chain_id,
                report.ligand_residue_seq,
                report.num_clashes,
                report.protein_volume_overlap_pct
            );
        }
    }

    Ok(())
}

/// Create a demo protein-ligand complex with both good and bad poses
fn create_demo_complex() -> PdbStructure {
    let mut structure = PdbStructure::new();

    // Create a small protein (simplified backbone)
    let protein_atoms = [
        // Residue 1 - ALA
        (1, "N", "ALA", 1, 0.0, 0.0, 0.0, "N"),
        (2, "CA", "ALA", 1, 1.458, 0.0, 0.0, "C"),
        (3, "C", "ALA", 1, 2.009, 1.420, 0.0, "C"),
        (4, "O", "ALA", 1, 1.251, 2.390, 0.0, "O"),
        (5, "CB", "ALA", 1, 1.988, -0.767, -1.199, "C"),
        // Residue 2 - VAL
        (6, "N", "VAL", 2, 3.303, 1.618, 0.0, "N"),
        (7, "CA", "VAL", 2, 3.920, 2.940, 0.0, "C"),
        (8, "C", "VAL", 2, 5.440, 2.840, 0.0, "C"),
        (9, "O", "VAL", 2, 6.040, 1.760, 0.0, "O"),
        (10, "CB", "VAL", 2, 3.450, 3.810, 1.170, "C"),
        // Residue 3 - LEU
        (11, "N", "LEU", 3, 6.030, 3.990, 0.0, "N"),
        (12, "CA", "LEU", 3, 7.480, 4.130, 0.0, "C"),
        (13, "C", "LEU", 3, 8.040, 5.530, 0.0, "C"),
        (14, "O", "LEU", 3, 7.250, 6.480, 0.0, "O"),
        (15, "CB", "LEU", 3, 8.010, 3.350, 1.200, "C"),
    ];

    for (serial, name, res_name, res_seq, x, y, z, element) in protein_atoms {
        structure.atoms.push(Atom::new(
            serial,
            name.to_string(),
            None,
            res_name.to_string(),
            "A".to_string(),
            res_seq,
            x,
            y,
            z,
            1.0,
            20.0,
            element.to_string(),
            None,
        ));
    }

    // Good ligand - well-positioned in binding site
    let good_ligand_atoms = [
        (100, "C1", 15.0, 5.0, 0.0, "C"),
        (101, "C2", 16.4, 5.0, 0.0, "C"),
        (102, "C3", 17.1, 6.2, 0.0, "C"),
        (103, "C4", 16.4, 7.4, 0.0, "C"),
        (104, "C5", 15.0, 7.4, 0.0, "C"),
        (105, "C6", 14.3, 6.2, 0.0, "C"),
        (106, "O1", 17.8, 6.2, 0.0, "O"),
        (107, "N1", 13.0, 6.2, 0.0, "N"),
    ];

    for (serial, name, x, y, z, element) in good_ligand_atoms {
        structure.atoms.push(Atom::new_hetatm(
            serial,
            name.to_string(),
            None,
            "LIG".to_string(),
            "A".to_string(),
            100,
            x,
            y,
            z,
            1.0,
            25.0,
            element.to_string(),
            None,
        ));
    }

    // Bad ligand - clashing with protein
    let bad_ligand_atoms = [
        (200, "C1", 1.8, 0.0, 0.0, "C"), // Very close to protein CA
        (201, "C2", 2.5, 1.0, 0.0, "C"), // Close to protein backbone
        (202, "C3", 3.5, 1.5, 0.0, "C"),
        (203, "O1", 4.5, 2.0, 0.0, "O"),
    ];

    for (serial, name, x, y, z, element) in bad_ligand_atoms {
        structure.atoms.push(Atom::new_hetatm(
            serial,
            name.to_string(),
            None,
            "BAD".to_string(),
            "A".to_string(),
            200,
            x,
            y,
            z,
            1.0,
            30.0,
            element.to_string(),
            None,
        ));
    }

    // Add some water (should be excluded from analysis)
    structure.atoms.push(Atom::new_hetatm(
        300,
        "O".to_string(),
        None,
        "HOH".to_string(),
        "A".to_string(),
        300,
        20.0,
        20.0,
        20.0,
        1.0,
        40.0,
        "O".to_string(),
        None,
    ));
    structure.atoms.push(Atom::new_hetatm(
        301,
        "O".to_string(),
        None,
        "HOH".to_string(),
        "A".to_string(),
        301,
        22.0,
        20.0,
        20.0,
        1.0,
        40.0,
        "O".to_string(),
        None,
    ));

    structure
}
