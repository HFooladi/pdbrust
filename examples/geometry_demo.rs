//! Geometry operations example: RMSD and structure alignment
//!
//! This example demonstrates:
//! - RMSD calculation between structures
//! - Structure alignment using Kabsch algorithm
//! - Different atom selections
//! - Per-residue RMSD for flexibility analysis
//!
//! Run with: cargo run --example geometry_demo --features geometry

use pdbrust::geometry::AtomSelection;
use pdbrust::parse_pdb_file;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("================================================");
    println!("PDBRust Geometry Demo: RMSD and Alignment");
    println!("================================================\n");

    // Load structure
    let structure1 = parse_pdb_file("examples/pdb_files/1UBQ.pdb")?;
    println!("Loaded structure: {} atoms", structure1.atoms.len());

    // Create a modified copy for comparison by loading again and translating
    let mut structure2 = parse_pdb_file("examples/pdb_files/1UBQ.pdb")?;

    // Translate structure2 to create a measurable difference
    for atom in &mut structure2.atoms {
        atom.x += 1.0;
        atom.y += 0.5;
        atom.z += 0.2;
    }
    println!("Created second structure by translating (1.0, 0.5, 0.2) Å\n");

    // --- Basic RMSD ---
    println!("1. BASIC RMSD CALCULATION");
    println!("-----------------------------------------");

    let rmsd = structure1.rmsd_to(&structure2)?;
    println!("RMSD (CA atoms, default): {:.4} Å", rmsd);

    // --- Different Atom Selections ---
    println!("\n2. RMSD WITH DIFFERENT ATOM SELECTIONS");
    println!("-----------------------------------------");

    // CA only (default)
    let rmsd_ca = structure1.rmsd_to_with_selection(&structure2, AtomSelection::CaOnly)?;
    println!("RMSD (CA only): {:.4} Å", rmsd_ca);

    // Backbone atoms
    let rmsd_bb = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;
    println!("RMSD (backbone): {:.4} Å", rmsd_bb);

    // All atoms
    let rmsd_all = structure1.rmsd_to_with_selection(&structure2, AtomSelection::AllAtoms)?;
    println!("RMSD (all atoms): {:.4} Å", rmsd_all);

    // Custom selection
    let rmsd_custom = structure1.rmsd_to_with_selection(
        &structure2,
        AtomSelection::Custom(vec!["CA".to_string(), "CB".to_string()]),
    )?;
    println!("RMSD (CA + CB): {:.4} Å", rmsd_custom);

    // --- Structure Alignment ---
    println!("\n3. STRUCTURE ALIGNMENT (KABSCH ALGORITHM)");
    println!("-----------------------------------------");

    // Align structure1 to structure2
    let (aligned, result) = structure1.align_to(&structure2)?;
    println!("Alignment RMSD: {:.4} Å", result.rmsd);
    println!("Atoms used: {}", result.num_atoms);

    // Verify alignment by computing RMSD after
    let rmsd_after = aligned.rmsd_to(&structure2)?;
    println!("RMSD after alignment: {:.6} Å (should be ~0)", rmsd_after);

    // Backbone alignment
    let (_aligned_bb, result_bb) =
        structure1.align_to_with_selection(&structure2, AtomSelection::Backbone)?;
    println!(
        "\nBackbone alignment: {:.4} Å ({} atoms)",
        result_bb.rmsd, result_bb.num_atoms
    );

    // --- Per-Residue RMSD ---
    println!("\n4. PER-RESIDUE RMSD (FLEXIBILITY ANALYSIS)");
    println!("-----------------------------------------");

    let per_res = structure1.per_residue_rmsd_to(&structure2)?;
    println!("Residues analyzed: {}", per_res.len());

    println!("\nFirst 10 residues:");
    println!(
        "{:<6} {:<8} {:<5} {:<10} {:<6}",
        "Chain", "ResSeq", "Name", "RMSD (Å)", "Atoms"
    );
    println!("{}", "-".repeat(40));

    for r in per_res.iter().take(10) {
        let (chain_id, residue_seq) = &r.residue_id;
        println!(
            "{:<6} {:<8} {:<5} {:<10.4} {:<6}",
            chain_id, residue_seq, r.residue_name, r.rmsd, r.num_atoms
        );
    }

    // Find flexible regions
    let threshold = 0.5;
    let flexible: Vec<_> = per_res.iter().filter(|r| r.rmsd > threshold).collect();
    if !flexible.is_empty() {
        println!(
            "\nFlexible regions (RMSD > {} Å): {}",
            threshold,
            flexible.len()
        );
        for r in flexible.iter().take(5) {
            let (chain_id, residue_seq) = &r.residue_id;
            println!(
                "  {}{} {}: {:.4} Å",
                chain_id, residue_seq, r.residue_name, r.rmsd
            );
        }
    } else {
        println!(
            "\nNo residues with RMSD > {} Å (structures are very similar)",
            threshold
        );
    }

    // --- Self Comparison ---
    println!("\n5. SELF-COMPARISON (VALIDATION)");
    println!("-----------------------------------------");

    let self_rmsd = structure1.rmsd_to(&structure1)?;
    println!("Self-RMSD: {:.10} Å (should be 0)", self_rmsd);

    let (_, self_result) = structure1.align_to(&structure1)?;
    println!(
        "Self-alignment RMSD: {:.10} Å (should be 0)",
        self_result.rmsd
    );

    // --- Summary ---
    println!("\n6. SUMMARY");
    println!("-----------------------------------------");
    println!(
        "
Structure comparison summary:
  - Structure 1: {} atoms
  - Structure 2: {} atoms (translated)
  - RMSD before alignment: {:.4} Å
  - RMSD after alignment: {:.6} Å
  - Atoms used in alignment: {}
",
        structure1.atoms.len(),
        structure2.atoms.len(),
        rmsd,
        rmsd_after,
        result.num_atoms
    );

    println!("================================================");
    println!("Geometry demo completed successfully!");
    println!("================================================");

    Ok(())
}
