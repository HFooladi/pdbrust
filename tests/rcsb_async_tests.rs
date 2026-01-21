//! Integration tests for the async RCSB module.
//!
//! Note: These tests require network access to the RCSB PDB servers.
//! Tests that require network are marked with #[ignore] and can be run
//! with `cargo test --features rcsb-async -- --ignored`

#![cfg(feature = "rcsb-async")]

use pdbrust::rcsb::{
    AsyncDownloadOptions, FileFormat, download_multiple_async, download_pdb_string_async,
    download_structure_async, download_to_file_async,
};
use tempfile::tempdir;

// ============================================================================
// AsyncDownloadOptions Tests (no network required)
// ============================================================================

#[test]
fn test_async_options_default() {
    let options = AsyncDownloadOptions::default();
    assert_eq!(options.max_concurrent, 5);
    assert_eq!(options.rate_limit_ms, 100);
    assert_eq!(options.timeout_secs, 30);
    assert_eq!(options.retries, 2);
}

#[test]
fn test_async_options_new() {
    let options = AsyncDownloadOptions::new();
    assert_eq!(options.max_concurrent, 5);
    assert_eq!(options.rate_limit_ms, 100);
}

#[test]
fn test_async_options_conservative() {
    let options = AsyncDownloadOptions::conservative();
    assert_eq!(options.max_concurrent, 2);
    assert_eq!(options.rate_limit_ms, 500);
    assert_eq!(options.timeout_secs, 60);
    assert_eq!(options.retries, 3);
}

#[test]
fn test_async_options_fast() {
    let options = AsyncDownloadOptions::fast();
    assert_eq!(options.max_concurrent, 20);
    assert_eq!(options.rate_limit_ms, 25);
    assert_eq!(options.timeout_secs, 30);
    assert_eq!(options.retries, 1);
}

#[test]
fn test_async_options_builder_max_concurrent() {
    let options = AsyncDownloadOptions::new().with_max_concurrent(10);
    assert_eq!(options.max_concurrent, 10);
    assert_eq!(options.rate_limit_ms, 100); // Other fields unchanged
}

#[test]
fn test_async_options_builder_rate_limit() {
    let options = AsyncDownloadOptions::new().with_rate_limit_ms(50);
    assert_eq!(options.rate_limit_ms, 50);
    assert_eq!(options.max_concurrent, 5); // Other fields unchanged
}

#[test]
fn test_async_options_builder_timeout() {
    let options = AsyncDownloadOptions::new().with_timeout_secs(60);
    assert_eq!(options.timeout_secs, 60);
}

#[test]
fn test_async_options_builder_retries() {
    let options = AsyncDownloadOptions::new().with_retries(5);
    assert_eq!(options.retries, 5);
}

#[test]
fn test_async_options_builder_chained() {
    let options = AsyncDownloadOptions::new()
        .with_max_concurrent(10)
        .with_rate_limit_ms(50)
        .with_timeout_secs(45)
        .with_retries(3);

    assert_eq!(options.max_concurrent, 10);
    assert_eq!(options.rate_limit_ms, 50);
    assert_eq!(options.timeout_secs, 45);
    assert_eq!(options.retries, 3);
}

#[test]
fn test_async_options_clone() {
    let options = AsyncDownloadOptions::conservative();
    let cloned = options.clone();
    assert_eq!(cloned.max_concurrent, options.max_concurrent);
    assert_eq!(cloned.rate_limit_ms, options.rate_limit_ms);
}

#[test]
fn test_async_options_debug() {
    let options = AsyncDownloadOptions::default();
    let debug = format!("{:?}", options);
    assert!(debug.contains("max_concurrent"));
    assert!(debug.contains("rate_limit_ms"));
}

// ============================================================================
// Network Tests (require network access)
// These are marked with #[ignore] and can be run with:
// cargo test --features rcsb-async -- --ignored
// ============================================================================

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_pdb_string_async_pdb_format() {
    let content = download_pdb_string_async("1UBQ", FileFormat::Pdb)
        .await
        .expect("Failed to download 1UBQ.pdb");

    // Check that it looks like a PDB file
    assert!(content.contains("HEADER"));
    assert!(content.contains("ATOM"));
    assert!(content.len() > 1000);
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_pdb_string_async_cif_format() {
    let content = download_pdb_string_async("1UBQ", FileFormat::Cif)
        .await
        .expect("Failed to download 1UBQ.cif");

    // Check that it looks like an mmCIF file
    assert!(content.contains("_atom_site"));
    assert!(content.len() > 1000);
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_structure_async_pdb_format() {
    let structure = download_structure_async("1UBQ", FileFormat::Pdb)
        .await
        .expect("Failed to download and parse 1UBQ.pdb");

    // 1UBQ (ubiquitin) has 76 residues and ~600 atoms
    assert!(structure.atoms.len() > 500);
    assert!(!structure.get_chain_ids().is_empty());
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_structure_async_cif_format() {
    let structure = download_structure_async("1UBQ", FileFormat::Cif)
        .await
        .expect("Failed to download and parse 1UBQ.cif");

    assert!(structure.atoms.len() > 500);
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_to_file_async() {
    let dir = tempdir().expect("Failed to create temp dir");
    let path = dir.path().join("1UBQ.pdb");

    download_to_file_async("1UBQ", &path, FileFormat::Pdb)
        .await
        .expect("Failed to download 1UBQ to file");

    // Verify file exists and has content
    assert!(path.exists());
    let content = std::fs::read_to_string(&path).expect("Failed to read file");
    assert!(content.contains("ATOM"));
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_async_nonexistent_pdb() {
    let result = download_pdb_string_async("XXXX", FileFormat::Pdb).await;
    assert!(result.is_err());
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_simple() {
    let pdb_ids = vec!["1UBQ", "8HM2"];
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;

    assert_eq!(results.len(), 2);

    // Check that results are in order
    assert_eq!(results[0].0, "1UBQ");
    assert_eq!(results[1].0, "8HM2");

    // Both should succeed
    for (pdb_id, result) in &results {
        assert!(
            result.is_ok(),
            "Failed to download {}: {:?}",
            pdb_id,
            result
        );
        let structure = result.as_ref().unwrap();
        assert!(structure.atoms.len() > 100, "{} has too few atoms", pdb_id);
    }
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_with_options() {
    let pdb_ids = vec!["1UBQ", "4INS"];
    let options = AsyncDownloadOptions::conservative();
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, Some(options)).await;

    assert_eq!(results.len(), 2);
    for (pdb_id, result) in &results {
        assert!(result.is_ok(), "Failed to download {}", pdb_id);
    }
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_partial_failure() {
    // Mix of valid and invalid PDB IDs
    let pdb_ids = vec!["1UBQ", "XXXX", "8HM2"];
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;

    assert_eq!(results.len(), 3);

    // 1UBQ should succeed
    assert!(results[0].1.is_ok());

    // XXXX should fail (not found)
    assert!(results[1].1.is_err());

    // 8HM2 should succeed
    assert!(results[2].1.is_ok());
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_larger_batch() {
    let pdb_ids = vec!["1UBQ", "8HM2", "4INS", "1HHB", "2MBP"];
    let options = AsyncDownloadOptions::default().with_max_concurrent(3);
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, Some(options)).await;

    assert_eq!(results.len(), 5);

    let successful: Vec<_> = results.iter().filter(|(_, r)| r.is_ok()).collect();
    // At least most should succeed (allow for transient failures)
    assert!(
        successful.len() >= 4,
        "Only {} of 5 downloads succeeded",
        successful.len()
    );
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_cif_format() {
    let pdb_ids = vec!["1UBQ", "8HM2"];
    let results = download_multiple_async(&pdb_ids, FileFormat::Cif, None).await;

    assert_eq!(results.len(), 2);
    for (pdb_id, result) in &results {
        assert!(result.is_ok(), "Failed to download {} as CIF", pdb_id);
    }
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_empty_list() {
    let pdb_ids: Vec<&str> = vec![];
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;
    assert!(results.is_empty());
}

#[tokio::test]
#[ignore = "requires network access"]
async fn test_download_multiple_async_single_item() {
    let pdb_ids = vec!["1UBQ"];
    let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;

    assert_eq!(results.len(), 1);
    assert!(results[0].1.is_ok());
}
