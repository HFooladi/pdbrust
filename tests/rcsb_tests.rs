//! Integration tests for the RCSB module.
//!
//! Note: These tests require network access to the RCSB PDB servers.
//! Tests that require network are marked with #[ignore] and can be run
//! with `cargo test --features rcsb -- --ignored`

#![cfg(feature = "rcsb")]

use pdbrust::rcsb::{
    ExperimentalMethod, FileFormat, PolymerType, SearchQuery, download_pdb_string,
    download_structure, download_to_file, rcsb_search,
};
use tempfile::tempdir;

// ============================================================================
// SearchQuery Builder Tests (no network required)
// ============================================================================

#[test]
fn test_search_query_new() {
    let query = SearchQuery::new();
    assert!(query.is_empty());
}

#[test]
fn test_search_query_with_text() {
    let query = SearchQuery::new().with_text("ubiquitin");
    assert_eq!(query.text, Some("ubiquitin".to_string()));
    assert!(!query.is_empty());
}

#[test]
fn test_search_query_with_organism() {
    let query = SearchQuery::new().with_organism("Homo sapiens");
    assert_eq!(query.organism, Some("Homo sapiens".to_string()));
}

#[test]
fn test_search_query_with_resolution() {
    let query = SearchQuery::new()
        .with_resolution_min(1.0)
        .with_resolution_max(2.5);
    assert_eq!(query.resolution_min, Some(1.0));
    assert_eq!(query.resolution_max, Some(2.5));
}

#[test]
fn test_search_query_with_experimental_method() {
    let query = SearchQuery::new().with_experimental_method(ExperimentalMethod::XRay);
    assert_eq!(query.experimental_method, Some(ExperimentalMethod::XRay));
}

#[test]
fn test_search_query_with_polymer_type() {
    let query = SearchQuery::new().with_polymer_type(PolymerType::Protein);
    assert_eq!(query.polymer_type, Some(PolymerType::Protein));
}

#[test]
fn test_search_query_with_sequence_length() {
    let query = SearchQuery::new()
        .with_sequence_length_min(50)
        .with_sequence_length_max(200);
    assert_eq!(query.sequence_length_min, Some(50));
    assert_eq!(query.sequence_length_max, Some(200));
}

#[test]
fn test_search_query_with_release_date() {
    let query = SearchQuery::new()
        .with_release_date_min("2020-01-01")
        .with_release_date_max("2024-12-31");
    assert_eq!(query.release_date_min, Some("2020-01-01".to_string()));
    assert_eq!(query.release_date_max, Some("2024-12-31".to_string()));
}

#[test]
fn test_search_query_with_ec_number() {
    let query = SearchQuery::new().with_ec_number("2.7.11.1");
    assert_eq!(query.ec_number, Some("2.7.11.1".to_string()));
}

#[test]
fn test_search_query_complex() {
    let query = SearchQuery::new()
        .with_text("kinase")
        .with_organism("Homo sapiens")
        .with_resolution_max(2.0)
        .with_experimental_method(ExperimentalMethod::XRay)
        .with_polymer_type(PolymerType::Protein);

    assert_eq!(query.text, Some("kinase".to_string()));
    assert_eq!(query.organism, Some("Homo sapiens".to_string()));
    assert_eq!(query.resolution_max, Some(2.0));
    assert!(!query.is_empty());
}

#[test]
fn test_search_query_to_json() {
    let query = SearchQuery::new().with_text("ubiquitin");
    let json = query.to_json();

    assert!(json.contains("ubiquitin"));
    assert!(json.contains("full_text"));
    assert!(json.contains("query"));
}

#[test]
fn test_search_query_to_json_multiple_criteria() {
    let query = SearchQuery::new()
        .with_text("kinase")
        .with_resolution_max(2.5);
    let json = query.to_json();

    assert!(json.contains("kinase"));
    assert!(json.contains("resolution_combined"));
    assert!(json.contains("logical_operator"));
}

// ============================================================================
// FileFormat Tests (no network required)
// ============================================================================

#[test]
fn test_file_format_extension() {
    assert_eq!(FileFormat::Pdb.extension(), "pdb");
    assert_eq!(FileFormat::Cif.extension(), "cif");
}

#[test]
fn test_file_format_compressed_extension() {
    assert_eq!(FileFormat::Pdb.compressed_extension(), "pdb.gz");
    assert_eq!(FileFormat::Cif.compressed_extension(), "cif.gz");
}

#[test]
fn test_file_format_display() {
    assert_eq!(FileFormat::Pdb.to_string(), "PDB");
    assert_eq!(FileFormat::Cif.to_string(), "mmCIF");
}

// ============================================================================
// ExperimentalMethod Tests (no network required)
// ============================================================================

#[test]
fn test_experimental_method_api_values() {
    assert_eq!(ExperimentalMethod::XRay.api_value(), "X-RAY DIFFRACTION");
    assert_eq!(ExperimentalMethod::Nmr.api_value(), "SOLUTION NMR");
    assert_eq!(ExperimentalMethod::Em.api_value(), "ELECTRON MICROSCOPY");
    assert_eq!(ExperimentalMethod::Other.api_value(), "OTHER");
}

// ============================================================================
// PolymerType Tests (no network required)
// ============================================================================

#[test]
fn test_polymer_type_api_values() {
    assert_eq!(PolymerType::Protein.api_value(), "Protein");
    assert_eq!(PolymerType::Dna.api_value(), "DNA");
    assert_eq!(PolymerType::Rna.api_value(), "RNA");
    assert_eq!(PolymerType::Hybrid.api_value(), "NA-hybrid");
}

// ============================================================================
// Network Tests (require network access)
// These are marked with #[ignore] and can be run with:
// cargo test --features rcsb -- --ignored
// ============================================================================

#[test]
#[ignore = "requires network access"]
fn test_download_pdb_string_pdb_format() {
    let content =
        download_pdb_string("1UBQ", FileFormat::Pdb).expect("Failed to download 1UBQ.pdb");

    // Check that it looks like a PDB file
    assert!(content.contains("HEADER"));
    assert!(content.contains("ATOM"));
    assert!(content.len() > 1000);
}

#[test]
#[ignore = "requires network access"]
fn test_download_pdb_string_cif_format() {
    let content =
        download_pdb_string("1UBQ", FileFormat::Cif).expect("Failed to download 1UBQ.cif");

    // Check that it looks like an mmCIF file
    assert!(content.contains("_atom_site"));
    assert!(content.len() > 1000);
}

#[test]
#[ignore = "requires network access"]
fn test_download_structure_pdb_format() {
    let structure =
        download_structure("1UBQ", FileFormat::Pdb).expect("Failed to download and parse 1UBQ.pdb");

    // 1UBQ (ubiquitin) has 76 residues and ~600 atoms
    assert!(structure.atoms.len() > 500);
    assert!(!structure.get_chain_ids().is_empty());
}

#[test]
#[ignore = "requires network access"]
fn test_download_structure_cif_format() {
    let structure =
        download_structure("1UBQ", FileFormat::Cif).expect("Failed to download and parse 1UBQ.cif");

    // Should have similar atom count as PDB format
    assert!(structure.atoms.len() > 500);
}

#[test]
#[ignore = "requires network access"]
fn test_download_to_file() {
    let dir = tempdir().expect("Failed to create temp dir");
    let path = dir.path().join("1UBQ.pdb");

    download_to_file("1UBQ", &path, FileFormat::Pdb).expect("Failed to download 1UBQ to file");

    // Verify file exists and has content
    assert!(path.exists());
    let content = std::fs::read_to_string(&path).expect("Failed to read file");
    assert!(content.contains("ATOM"));
}

#[test]
#[ignore = "requires network access"]
fn test_download_nonexistent_pdb() {
    let result = download_pdb_string("XXXX", FileFormat::Pdb);
    assert!(result.is_err());
}

#[test]
#[ignore = "requires network access"]
fn test_rcsb_search_simple() {
    let query = SearchQuery::new().with_text("ubiquitin");
    let result = rcsb_search(&query, 5).expect("Search failed");

    // Should find some structures
    assert!(!result.pdb_ids.is_empty());
    assert!(result.pdb_ids.len() <= 5);
    assert!(result.total_count > 0);
}

#[test]
#[ignore = "requires network access"]
fn test_rcsb_search_with_resolution() {
    let query = SearchQuery::new()
        .with_text("ubiquitin")
        .with_resolution_max(1.5);

    let result = rcsb_search(&query, 10).expect("Search failed");

    // High-resolution ubiquitin structures should exist
    assert!(result.total_count > 0);
}

#[test]
#[ignore = "requires network access"]
fn test_rcsb_search_with_organism() {
    let query = SearchQuery::new()
        .with_text("insulin")
        .with_organism("Homo sapiens");

    let result = rcsb_search(&query, 10).expect("Search failed");

    // Human insulin structures should exist
    assert!(result.total_count > 0);
}

#[test]
#[ignore = "requires network access"]
fn test_rcsb_search_xray_only() {
    let query = SearchQuery::new()
        .with_text("hemoglobin")
        .with_experimental_method(ExperimentalMethod::XRay);

    let result = rcsb_search(&query, 10).expect("Search failed");

    // X-ray hemoglobin structures should exist
    assert!(result.total_count > 0);
}

#[test]
#[ignore = "requires network access"]
fn test_search_and_download_workflow() {
    // Search for small ubiquitin structures
    let query = SearchQuery::new()
        .with_text("ubiquitin")
        .with_resolution_max(2.0);

    let result = rcsb_search(&query, 1).expect("Search failed");
    assert!(!result.pdb_ids.is_empty());

    // Download the first result
    let pdb_id = &result.pdb_ids[0];
    let structure =
        download_structure(pdb_id, FileFormat::Pdb).expect("Failed to download structure");

    // Verify we got a valid structure
    assert!(structure.atoms.len() > 100);
}
