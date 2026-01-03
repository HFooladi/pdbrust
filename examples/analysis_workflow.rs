//! Complete Analysis Workflow Example
//!
//! This example demonstrates a typical cheminformatics workflow:
//! 1. Load a structure (auto-detect PDB/mmCIF)
//! 2. Assess quality
//! 3. Clean the structure
//! 4. Compute structural descriptors
//! 5. Generate unified summary
//! 6. Export cleaned structure
//!
//! Run with:
//! ```bash
//! cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"
//! cargo run --example analysis_workflow --features "filter,descriptors,quality,summary" -- path/to/file.pdb
//! ```

use pdbrust::{parse_structure_file, write_pdb_file};
use std::env;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Get input file from command line or use default
    let args: Vec<String> = env::args().collect();
    let pdb_file = args
        .get(1)
        .map(|s| s.as_str())
        .unwrap_or("examples/pdb_files/1UBQ.pdb");

    println!("=== PDBRust Analysis Workflow ===\n");

    // ========== Step 1: Load Structure ==========
    println!("Step 1: Loading structure from '{}'", pdb_file);
    let structure = parse_structure_file(pdb_file)?;

    println!(
        "  Loaded: {} atoms, {} chains ({:?})",
        structure.atoms.len(),
        structure.get_chain_ids().len(),
        structure.get_chain_ids()
    );
    if let Some(title) = &structure.title {
        println!("  Title: {}", title);
    }

    // ========== Step 2: Quality Assessment ==========
    println!("\nStep 2: Quality Assessment");
    let report = structure.quality_report();

    println!("  Chains: {}", report.num_chains);
    println!("  Models: {}", report.num_models);
    println!("  Atoms: {}", report.num_atoms);
    println!("  Has HETATM (ligands/waters): {}", report.has_hetatm);
    println!("  Has hydrogens: {}", report.has_hydrogens);
    println!("  Has alternate locations: {}", report.has_altlocs);
    println!("  Has disulfide bonds: {}", report.has_ssbonds);
    println!(
        "  Analysis ready: {}",
        if report.is_analysis_ready() {
            "Yes"
        } else {
            "No (see flags above)"
        }
    );

    // ========== Step 3: Clean Structure ==========
    println!("\nStep 3: Cleaning Structure");

    // Remove ligands and waters (common preprocessing step)
    let cleaned = structure.remove_ligands();
    println!(
        "  After removing ligands/waters: {} atoms (removed {})",
        cleaned.atoms.len(),
        structure.atoms.len() - cleaned.atoms.len()
    );

    // Optionally remove hydrogens for lighter representation
    let cleaned = cleaned.remove_hydrogens();
    println!("  After removing hydrogens: {} atoms", cleaned.atoms.len());

    // ========== Step 4: Compute Structural Descriptors ==========
    println!("\nStep 4: Structural Descriptors");

    let rg = cleaned.radius_of_gyration();
    println!("  Radius of gyration: {:.2} A", rg);

    let max_dist = cleaned.max_ca_distance();
    println!("  Max CA-CA distance: {:.2} A", max_dist);

    let n_residues = cleaned.count_ca_residues();
    println!("  Number of residues: {}", n_residues);

    // Amino acid composition
    let composition = cleaned.aa_composition();
    println!("\n  Amino acid composition (top 5):");
    let mut sorted: Vec<_> = composition.iter().collect();
    sorted.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
    for (aa, fraction) in sorted.iter().take(5) {
        println!("    {}: {:.1}%", aa, *fraction * 100.0);
    }

    // Specific ratios
    let hydrophobic = cleaned.hydrophobic_ratio();
    let polar = cleaned.polar_ratio();
    let charged = cleaned.charged_ratio();
    println!("\n  Composition ratios:");
    println!("    Hydrophobic: {:.1}%", hydrophobic * 100.0);
    println!("    Polar: {:.1}%", polar * 100.0);
    println!("    Charged: {:.1}%", charged * 100.0);

    // ========== Step 5: Unified Summary ==========
    println!("\nStep 5: Unified Summary");
    let summary = cleaned.summary();

    println!("  Residues: {}", summary.num_residues);
    println!("  Atoms: {}", summary.num_atoms);
    println!("  Rg: {:.2} A", summary.radius_of_gyration);
    println!("  Compactness index: {:.3}", summary.compactness_index);
    println!("  CA density: {:.6} atoms/A^3", summary.ca_density);
    println!(
        "  Analysis ready: {}",
        if summary.is_analysis_ready() {
            "Yes"
        } else {
            "No"
        }
    );

    // Export summary as CSV row
    let csv_header = pdbrust::summary::StructureSummary::field_names().join(",");
    let csv_row = summary.to_csv_values().join(",");
    println!("\n  CSV export:");
    println!("  Header: {}", csv_header);
    println!("  Values: {}", csv_row);

    // ========== Step 6: Export Cleaned Structure ==========
    let output_file = "cleaned_output.pdb";
    println!("\nStep 6: Exporting cleaned structure to '{}'", output_file);
    write_pdb_file(&cleaned, output_file)?;
    println!("  Done! Cleaned structure saved.");

    // ========== Summary ==========
    println!("\n=== Workflow Complete ===");
    println!("Input:  {} atoms", structure.atoms.len());
    println!("Output: {} atoms (cleaned)", cleaned.atoms.len());
    println!(
        "Reduction: {:.1}%",
        (1.0 - cleaned.atoms.len() as f64 / structure.atoms.len() as f64) * 100.0
    );

    Ok(())
}
