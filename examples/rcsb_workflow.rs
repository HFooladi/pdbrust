//! RCSB PDB Search and Download Workflow
//!
//! This example demonstrates how to:
//! 1. Build search queries for RCSB PDB
//! 2. Search for structures matching criteria
//! 3. Download and analyze matching structures
//!
//! **Note**: This example requires network access to RCSB PDB.
//!
//! Run with:
//! ```bash
//! cargo run --example rcsb_workflow --features "rcsb,descriptors"
//! ```

use pdbrust::rcsb::{
    ExperimentalMethod, FileFormat, PolymerType, SearchQuery, download_structure, download_to_file,
    rcsb_search,
};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== RCSB PDB Workflow ===\n");

    // ========== Example 1: Simple Text Search ==========
    println!("--- Example 1: Simple Text Search ---");

    let query = SearchQuery::new().with_text("ubiquitin");

    println!("Searching for: 'ubiquitin'");
    match rcsb_search(&query, 5) {
        Ok(result) => {
            println!("Found {} total structures", result.total_count);
            println!("Top 5 PDB IDs: {:?}", result.pdb_ids);
        }
        Err(e) => println!("Search failed: {}", e),
    }

    // ========== Example 2: Complex Query with Multiple Filters ==========
    println!("\n--- Example 2: Complex Query ---");

    let query = SearchQuery::new()
        .with_text("kinase")
        .with_organism("Homo sapiens")
        .with_resolution_max(2.0)
        .with_experimental_method(ExperimentalMethod::XRay)
        .with_polymer_type(PolymerType::Protein);

    println!("Searching for:");
    println!("  - Text: 'kinase'");
    println!("  - Organism: Homo sapiens");
    println!("  - Resolution: <= 2.0 A");
    println!("  - Method: X-ray crystallography");
    println!("  - Type: Protein");

    match rcsb_search(&query, 5) {
        Ok(result) => {
            println!("\nFound {} total structures", result.total_count);
            println!("Top 5: {:?}", result.pdb_ids);
        }
        Err(e) => println!("Search failed: {}", e),
    }

    // ========== Example 3: Date Range Search ==========
    println!("\n--- Example 3: Recent Structures ---");

    let query = SearchQuery::new()
        .with_release_date_min("2023-01-01")
        .with_resolution_max(1.5)
        .with_experimental_method(ExperimentalMethod::XRay);

    println!("Searching for:");
    println!("  - Released after: 2023-01-01");
    println!("  - Resolution: <= 1.5 A (high resolution)");
    println!("  - Method: X-ray");

    match rcsb_search(&query, 5) {
        Ok(result) => {
            println!("\nFound {} structures", result.total_count);
            println!("Top 5: {:?}", result.pdb_ids);
        }
        Err(e) => println!("Search failed: {}", e),
    }

    // ========== Example 4: Enzyme Classification Search ==========
    println!("\n--- Example 4: Enzyme Classification ---");

    let query = SearchQuery::new()
        .with_ec_number("2.7.11.1") // Protein kinases
        .with_resolution_max(2.5);

    println!("Searching for:");
    println!("  - EC number: 2.7.11.1 (protein kinases)");
    println!("  - Resolution: <= 2.5 A");

    match rcsb_search(&query, 5) {
        Ok(result) => {
            println!("\nFound {} protein kinases", result.total_count);
            println!("Top 5: {:?}", result.pdb_ids);
        }
        Err(e) => println!("Search failed: {}", e),
    }

    // ========== Example 5: Download and Analyze ==========
    println!("\n--- Example 5: Download and Analyze ---");

    let pdb_id = "1UBQ"; // Ubiquitin - a small, well-known protein

    println!("Downloading {} in PDB format...", pdb_id);
    match download_structure(pdb_id, FileFormat::Pdb) {
        Ok(structure) => {
            println!("Successfully downloaded!");
            println!("  Atoms: {}", structure.atoms.len());
            println!("  Chains: {:?}", structure.get_chain_ids());

            if let Some(title) = &structure.title {
                println!("  Title: {}", title);
            }

            // Compute descriptors (requires 'descriptors' feature)
            #[cfg(feature = "descriptors")]
            {
                let rg = structure.radius_of_gyration();
                let max_dist = structure.max_ca_distance();
                println!("\n  Structural descriptors:");
                println!("    Rg: {:.2} A", rg);
                println!("    Max CA distance: {:.2} A", max_dist);
            }
        }
        Err(e) => println!("Download failed: {}", e),
    }

    // ========== Example 6: Download to File ==========
    println!("\n--- Example 6: Download to File ---");

    println!("Downloading 1UBQ to file...");
    match download_to_file("1UBQ", "downloaded_1UBQ.pdb", FileFormat::Pdb) {
        Ok(()) => println!("Saved to: downloaded_1UBQ.pdb"),
        Err(e) => println!("Download failed: {}", e),
    }

    // ========== Example 7: Download mmCIF Format ==========
    println!("\n--- Example 7: mmCIF Format ---");

    println!("Downloading 1UBQ in mmCIF format...");
    match download_structure("1UBQ", FileFormat::Cif) {
        Ok(structure) => {
            println!("Successfully downloaded mmCIF!");
            println!("  Atoms: {}", structure.atoms.len());
        }
        Err(e) => println!("Download failed: {}", e),
    }

    // ========== Example 8: Search and Download Pipeline ==========
    println!("\n--- Example 8: Search -> Download -> Analyze Pipeline ---");

    // Search for small, high-resolution structures
    let query = SearchQuery::new()
        .with_text("insulin")
        .with_resolution_max(1.5)
        .with_sequence_length_max(100);

    println!("Searching for small insulin structures...");

    match rcsb_search(&query, 3) {
        Ok(result) => {
            println!(
                "Found {} structures, analyzing top 3:\n",
                result.total_count
            );

            for pdb_id in result.pdb_ids.iter().take(3) {
                print!("  {}: ", pdb_id);
                match download_structure(pdb_id, FileFormat::Pdb) {
                    Ok(structure) => {
                        #[cfg(feature = "descriptors")]
                        {
                            let rg = structure.radius_of_gyration();
                            println!(
                                "{} atoms, {} residues, Rg={:.1} A",
                                structure.atoms.len(),
                                structure.count_ca_residues(),
                                rg
                            );
                        }
                        #[cfg(not(feature = "descriptors"))]
                        {
                            println!("{} atoms", structure.atoms.len());
                        }
                    }
                    Err(e) => println!("failed: {}", e),
                }
            }
        }
        Err(e) => println!("Search failed: {}", e),
    }

    println!("\n=== Workflow Complete ===");
    println!("\nTip: Use search queries to find structures, then download");
    println!("and analyze them programmatically for large-scale studies.");

    Ok(())
}
