//! Batch Processing Example
//!
//! This example demonstrates how to process multiple PDB files
//! in batch and export results to CSV. Common use cases:
//! - Dataset characterization
//! - Quality filtering
//! - Feature extraction for ML
//!
//! Run with:
//! ```bash
//! cargo run --example batch_processing --features "descriptors,summary"
//! ```

use pdbrust::parse_structure_file;
use pdbrust::summary::{batch_summarize, summaries_to_csv, StructureSummary};
use std::error::Error;
use std::fs;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== PDBRust Batch Processing ===\n");

    let pdb_dir = "examples/pdb_files";

    // ========== Step 1: Find PDB Files ==========
    println!("Step 1: Finding PDB files in '{}'", pdb_dir);

    let pdb_files: Vec<_> = fs::read_dir(pdb_dir)?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension()
                .map_or(false, |ext| ext == "pdb" || ext == "cif")
        })
        .collect();

    println!("Found {} structure files:", pdb_files.len());
    for path in &pdb_files {
        println!("  - {}", path.file_name().unwrap().to_string_lossy());
    }

    // ========== Step 2: Parse All Structures ==========
    println!("\nStep 2: Parsing structures...");

    let mut structures = Vec::new();
    let mut filenames = Vec::new();

    for path in &pdb_files {
        let filename = path.file_name().unwrap().to_string_lossy().to_string();
        print!("  Parsing {}... ", filename);

        match parse_structure_file(path) {
            Ok(structure) => {
                println!("{} atoms", structure.atoms.len());
                structures.push(structure);
                filenames.push(filename);
            }
            Err(e) => {
                println!("ERROR: {}", e);
            }
        }
    }

    println!("Successfully parsed {} structures", structures.len());

    if structures.is_empty() {
        println!("\nNo structures to process. Exiting.");
        return Ok(());
    }

    // ========== Step 3: Compute Summaries ==========
    println!("\nStep 3: Computing summaries...");

    let summaries = batch_summarize(&structures);

    // Print summary table
    println!("\n{:<15} {:>8} {:>8} {:>8} {:>10}", "File", "Atoms", "Residues", "Chains", "Rg (A)");
    println!("{}", "-".repeat(55));

    for (filename, summary) in filenames.iter().zip(summaries.iter()) {
        println!(
            "{:<15} {:>8} {:>8} {:>8} {:>10.2}",
            truncate_filename(filename, 15),
            summary.num_atoms,
            summary.num_residues,
            summary.num_chains,
            summary.radius_of_gyration
        );
    }

    // ========== Step 4: Quality Filtering ==========
    println!("\nStep 4: Quality filtering...");

    let analysis_ready: Vec<_> = summaries
        .iter()
        .zip(filenames.iter())
        .filter(|(s, _)| s.is_analysis_ready())
        .collect();

    println!(
        "Structures ready for analysis: {}/{}",
        analysis_ready.len(),
        summaries.len()
    );

    for (_, filename) in &analysis_ready {
        println!("  - {}", filename);
    }

    // Filter by specific criteria
    let high_quality: Vec<_> = summaries
        .iter()
        .zip(filenames.iter())
        .filter(|(s, _)| s.num_residues >= 50 && s.radius_of_gyration > 10.0)
        .collect();

    println!(
        "\nStructures with >= 50 residues and Rg > 10A: {}",
        high_quality.len()
    );

    // ========== Step 5: Export to CSV ==========
    println!("\nStep 5: Exporting to CSV...");

    // Add filename column manually
    let csv_output = create_csv_with_filenames(&filenames, &summaries);

    let output_file = "batch_results.csv";
    fs::write(output_file, &csv_output)?;
    println!("Saved to: {}", output_file);

    // Also create standard CSV (without filenames)
    let standard_csv = summaries_to_csv(&summaries, true);
    fs::write("batch_summaries.csv", &standard_csv)?;
    println!("Saved standard format to: batch_summaries.csv");

    // ========== Step 6: Statistics ==========
    println!("\nStep 6: Dataset statistics");

    if !summaries.is_empty() {
        let total_atoms: usize = summaries.iter().map(|s| s.num_atoms).sum();
        let total_residues: usize = summaries.iter().map(|s| s.num_residues).sum();
        let avg_rg: f64 =
            summaries.iter().map(|s| s.radius_of_gyration).sum::<f64>() / summaries.len() as f64;
        let avg_hydrophobic: f64 =
            summaries.iter().map(|s| s.hydrophobic_ratio).sum::<f64>() / summaries.len() as f64;

        println!("  Total atoms: {}", total_atoms);
        println!("  Total residues: {}", total_residues);
        println!("  Average Rg: {:.2} A", avg_rg);
        println!("  Average hydrophobic ratio: {:.1}%", avg_hydrophobic * 100.0);

        // Size distribution
        let sizes: Vec<usize> = summaries.iter().map(|s| s.num_residues).collect();
        let min_size = sizes.iter().min().unwrap_or(&0);
        let max_size = sizes.iter().max().unwrap_or(&0);
        println!("  Size range: {} - {} residues", min_size, max_size);
    }

    println!("\n=== Batch Processing Complete ===");

    Ok(())
}

/// Truncate filename for display
fn truncate_filename(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len - 3])
    }
}

/// Create CSV with filename as first column
fn create_csv_with_filenames(filenames: &[String], summaries: &[StructureSummary]) -> String {
    let mut output = String::new();

    // Header
    output.push_str("filename,");
    output.push_str(&StructureSummary::field_names().join(","));
    output.push('\n');

    // Data rows
    for (filename, summary) in filenames.iter().zip(summaries.iter()) {
        output.push_str(filename);
        output.push(',');
        output.push_str(&summary.to_csv_values().join(","));
        output.push('\n');
    }

    output
}
