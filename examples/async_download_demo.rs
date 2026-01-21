//! Async RCSB PDB Download Demo
//!
//! This example demonstrates how to efficiently download multiple structures
//! from RCSB PDB using async/concurrent downloads with rate limiting.
//!
//! **Note**: This example requires network access to RCSB PDB.
//!
//! Run with:
//! ```bash
//! cargo run --example async_download_demo --features "rcsb-async,descriptors"
//! ```

use pdbrust::rcsb::{
    AsyncDownloadOptions, FileFormat, download_multiple_async, download_structure_async,
};
use std::error::Error;
use std::time::Instant;

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    println!("=== Async RCSB Download Demo ===\n");

    // ========== Example 1: Single Async Download ==========
    println!("--- Example 1: Single Async Download ---");

    let start = Instant::now();
    match download_structure_async("1UBQ", FileFormat::Pdb).await {
        Ok(structure) => {
            println!(
                "Downloaded 1UBQ: {} atoms in {:?}",
                structure.atoms.len(),
                start.elapsed()
            );
        }
        Err(e) => println!("Failed to download 1UBQ: {}", e),
    }

    // ========== Example 2: Multiple Downloads with Default Options ==========
    println!("\n--- Example 2: Multiple Downloads (Default Options) ---");

    let pdb_ids = vec!["1UBQ", "8HM2", "4INS", "1HHB", "2MBP"];
    println!(
        "Downloading {} structures with default options (5 concurrent, 100ms delay)...",
        pdb_ids.len()
    );

    let start = Instant::now();
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;
    let elapsed = start.elapsed();

    let successful: Vec<_> = results.iter().filter(|(_, r)| r.is_ok()).collect();
    let failed: Vec<_> = results.iter().filter(|(_, r)| r.is_err()).collect();

    println!("\nResults ({:?}):", elapsed);
    for (pdb_id, result) in &results {
        match result {
            Ok(structure) => {
                #[cfg(feature = "descriptors")]
                {
                    let rg = structure.radius_of_gyration();
                    println!(
                        "  {}: {} atoms, Rg={:.1} A",
                        pdb_id,
                        structure.atoms.len(),
                        rg
                    );
                }
                #[cfg(not(feature = "descriptors"))]
                {
                    println!("  {}: {} atoms", pdb_id, structure.atoms.len());
                }
            }
            Err(e) => println!("  {}: FAILED - {}", pdb_id, e),
        }
    }
    println!(
        "\nSummary: {}/{} succeeded, {}/{} failed",
        successful.len(),
        pdb_ids.len(),
        failed.len(),
        pdb_ids.len()
    );

    // ========== Example 3: Conservative Options ==========
    println!("\n--- Example 3: Conservative Options ---");

    let pdb_ids = vec!["1CRN", "1L2Y", "1VII"];
    let options = AsyncDownloadOptions::conservative();
    println!(
        "Using conservative options: {} concurrent, {}ms delay",
        options.max_concurrent, options.rate_limit_ms
    );

    let start = Instant::now();
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, Some(options)).await;
    let elapsed = start.elapsed();

    println!(
        "\nDownloaded {} structures in {:?}:",
        pdb_ids.len(),
        elapsed
    );
    for (pdb_id, result) in &results {
        match result {
            Ok(structure) => println!("  {}: {} atoms", pdb_id, structure.atoms.len()),
            Err(e) => println!("  {}: FAILED - {}", pdb_id, e),
        }
    }

    // ========== Example 4: Fast Options (use responsibly) ==========
    println!("\n--- Example 4: Fast Options ---");

    let pdb_ids = vec!["3PQR", "1TIM", "1AKE"];
    let options = AsyncDownloadOptions::fast();
    println!(
        "Using fast options: {} concurrent, {}ms delay",
        options.max_concurrent, options.rate_limit_ms
    );
    println!("(Use responsibly with appropriate API access)");

    let start = Instant::now();
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, Some(options)).await;
    let elapsed = start.elapsed();

    println!(
        "\nDownloaded {} structures in {:?}:",
        pdb_ids.len(),
        elapsed
    );
    for (pdb_id, result) in &results {
        match result {
            Ok(structure) => println!("  {}: {} atoms", pdb_id, structure.atoms.len()),
            Err(e) => println!("  {}: FAILED - {}", pdb_id, e),
        }
    }

    // ========== Example 5: Custom Options ==========
    println!("\n--- Example 5: Custom Options ---");

    let pdb_ids = vec!["1BNA", "1D66", "5DNA"];
    let options = AsyncDownloadOptions::default()
        .with_max_concurrent(3)
        .with_rate_limit_ms(200)
        .with_timeout_secs(60)
        .with_retries(3);

    println!(
        "Custom options: {} concurrent, {}ms delay, {}s timeout, {} retries",
        options.max_concurrent, options.rate_limit_ms, options.timeout_secs, options.retries
    );

    let start = Instant::now();
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, Some(options)).await;
    let elapsed = start.elapsed();

    println!(
        "\nDownloaded {} structures in {:?}:",
        pdb_ids.len(),
        elapsed
    );
    for (pdb_id, result) in &results {
        match result {
            Ok(structure) => println!("  {}: {} atoms", pdb_id, structure.atoms.len()),
            Err(e) => println!("  {}: FAILED - {}", pdb_id, e),
        }
    }

    // ========== Example 6: Mixed Success/Failure ==========
    println!("\n--- Example 6: Handling Mixed Results ---");

    let pdb_ids = vec!["1UBQ", "XXXX", "8HM2", "YYYY", "4INS"];
    println!("Downloading mix of valid and invalid PDB IDs...");

    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;

    let successful: Vec<_> = results
        .iter()
        .filter_map(|(id, r)| r.as_ref().ok().map(|s| (id, s)))
        .collect();
    let failed: Vec<_> = results
        .iter()
        .filter_map(|(id, r)| r.as_ref().err().map(|e| (id, e)))
        .collect();

    println!("\nSuccessful downloads:");
    for (pdb_id, structure) in &successful {
        println!("  {}: {} atoms", pdb_id, structure.atoms.len());
    }

    println!("\nFailed downloads:");
    for (pdb_id, error) in &failed {
        println!("  {}: {}", pdb_id, error);
    }

    // ========== Example 7: mmCIF Format ==========
    println!("\n--- Example 7: mmCIF Format ---");

    let pdb_ids = vec!["1UBQ", "4INS"];
    println!("Downloading in mmCIF format...");

    let results = download_multiple_async(&pdb_ids, FileFormat::Cif, None).await;

    for (pdb_id, result) in &results {
        match result {
            Ok(structure) => println!("  {}.cif: {} atoms", pdb_id, structure.atoms.len()),
            Err(e) => println!("  {}.cif: FAILED - {}", pdb_id, e),
        }
    }

    // ========== Example 8: Larger Batch ==========
    println!("\n--- Example 8: Larger Batch Download ---");

    // A diverse set of well-known proteins
    let pdb_ids = vec![
        "1UBQ", // Ubiquitin
        "1CRN", // Crambin
        "1L2Y", // Trp-cage
        "1VII", // Villin headpiece
        "2MBP", // Maltose binding protein
        "1HHB", // Hemoglobin
        "4INS", // Insulin
        "1TIM", // Triose phosphate isomerase
        "3PQR", // Cytochrome c
        "1AKE", // Adenylate kinase
    ];

    let options = AsyncDownloadOptions::default()
        .with_max_concurrent(5)
        .with_rate_limit_ms(100);

    println!(
        "Downloading {} structures ({} concurrent)...",
        pdb_ids.len(),
        options.max_concurrent
    );

    let start = Instant::now();
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, Some(options)).await;
    let elapsed = start.elapsed();

    let total_atoms: usize = results
        .iter()
        .filter_map(|(_, r)| r.as_ref().ok())
        .map(|s| s.atoms.len())
        .sum();

    let success_count = results.iter().filter(|(_, r)| r.is_ok()).count();

    println!(
        "\nCompleted in {:?}: {}/{} successful",
        elapsed,
        success_count,
        pdb_ids.len()
    );
    println!(
        "Total atoms downloaded: {} ({:.0} atoms/sec)",
        total_atoms,
        total_atoms as f64 / elapsed.as_secs_f64()
    );

    println!("\n=== Demo Complete ===");
    println!("\nTips:");
    println!("  - Use conservative() for rate-limited scenarios");
    println!("  - Use fast() only with appropriate API access");
    println!("  - Always handle both success and failure cases");
    println!("  - Consider your network and RCSB's rate limits");

    Ok(())
}
