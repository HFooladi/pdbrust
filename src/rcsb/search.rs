//! RCSB PDB Search API v2 client.
//!
//! This module provides functionality for searching the RCSB PDB using their
//! Search API v2.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::fmt;

use super::RCSB_SEARCH_URL;

/// Errors that can occur during RCSB search operations.
#[derive(Debug)]
pub enum SearchError {
    /// HTTP request failed
    RequestFailed(String),
    /// Failed to parse response
    ParseError(String),
    /// API returned an error
    ApiError(String),
    /// Invalid query parameters
    InvalidQuery(String),
}

impl fmt::Display for SearchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SearchError::RequestFailed(msg) => write!(f, "Request failed: {}", msg),
            SearchError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            SearchError::ApiError(msg) => write!(f, "API error: {}", msg),
            SearchError::InvalidQuery(msg) => write!(f, "Invalid query: {}", msg),
        }
    }
}

impl std::error::Error for SearchError {}

impl From<reqwest::Error> for SearchError {
    fn from(err: reqwest::Error) -> Self {
        SearchError::RequestFailed(err.to_string())
    }
}

impl From<serde_json::Error> for SearchError {
    fn from(err: serde_json::Error) -> Self {
        SearchError::ParseError(err.to_string())
    }
}

/// Experimental method filter for searches.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExperimentalMethod {
    /// X-ray crystallography
    XRay,
    /// Nuclear Magnetic Resonance
    Nmr,
    /// Electron Microscopy
    Em,
    /// Other methods
    Other,
}

impl ExperimentalMethod {
    /// Get the RCSB API value for this method.
    pub fn api_value(&self) -> &'static str {
        match self {
            ExperimentalMethod::XRay => "X-RAY DIFFRACTION",
            ExperimentalMethod::Nmr => "SOLUTION NMR",
            ExperimentalMethod::Em => "ELECTRON MICROSCOPY",
            ExperimentalMethod::Other => "OTHER",
        }
    }
}

/// Polymer type filter for searches.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PolymerType {
    /// Protein structures
    Protein,
    /// DNA structures
    Dna,
    /// RNA structures
    Rna,
    /// Hybrid structures
    Hybrid,
}

impl PolymerType {
    /// Get the RCSB API value for this polymer type.
    pub fn api_value(&self) -> &'static str {
        match self {
            PolymerType::Protein => "Protein",
            PolymerType::Dna => "DNA",
            PolymerType::Rna => "RNA",
            PolymerType::Hybrid => "NA-hybrid",
        }
    }
}

/// Builder for RCSB search queries.
///
/// This struct provides a fluent interface for building search queries
/// against the RCSB PDB Search API v2.
///
/// # Examples
///
/// ```
/// use pdbrust::rcsb::SearchQuery;
///
/// let query = SearchQuery::new()
///     .with_text("kinase")
///     .with_organism("Homo sapiens")
///     .with_resolution_max(2.0);
/// ```
#[derive(Debug, Clone, Default)]
pub struct SearchQuery {
    /// Full-text search term
    pub text: Option<String>,
    /// Organism/source filter (scientific name)
    pub organism: Option<String>,
    /// Maximum resolution in Angstroms
    pub resolution_max: Option<f64>,
    /// Minimum resolution in Angstroms
    pub resolution_min: Option<f64>,
    /// Experimental method filter
    pub experimental_method: Option<ExperimentalMethod>,
    /// Polymer type filter
    pub polymer_type: Option<PolymerType>,
    /// Minimum release date (YYYY-MM-DD format)
    pub release_date_min: Option<String>,
    /// Maximum release date (YYYY-MM-DD format)
    pub release_date_max: Option<String>,
    /// Minimum sequence length
    pub sequence_length_min: Option<usize>,
    /// Maximum sequence length
    pub sequence_length_max: Option<usize>,
    /// PDB ID prefix filter
    pub pdb_id_prefix: Option<String>,
    /// Enzyme classification (EC) number
    pub ec_number: Option<String>,
}

impl SearchQuery {
    /// Create a new empty search query.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a full-text search term.
    pub fn with_text(mut self, text: &str) -> Self {
        self.text = Some(text.to_string());
        self
    }

    /// Filter by source organism (scientific name).
    pub fn with_organism(mut self, organism: &str) -> Self {
        self.organism = Some(organism.to_string());
        self
    }

    /// Set maximum resolution (Angstroms).
    pub fn with_resolution_max(mut self, resolution: f64) -> Self {
        self.resolution_max = Some(resolution);
        self
    }

    /// Set minimum resolution (Angstroms).
    pub fn with_resolution_min(mut self, resolution: f64) -> Self {
        self.resolution_min = Some(resolution);
        self
    }

    /// Filter by experimental method.
    pub fn with_experimental_method(mut self, method: ExperimentalMethod) -> Self {
        self.experimental_method = Some(method);
        self
    }

    /// Filter by polymer type.
    pub fn with_polymer_type(mut self, polymer_type: PolymerType) -> Self {
        self.polymer_type = Some(polymer_type);
        self
    }

    /// Set minimum release date (YYYY-MM-DD format).
    pub fn with_release_date_min(mut self, date: &str) -> Self {
        self.release_date_min = Some(date.to_string());
        self
    }

    /// Set maximum release date (YYYY-MM-DD format).
    pub fn with_release_date_max(mut self, date: &str) -> Self {
        self.release_date_max = Some(date.to_string());
        self
    }

    /// Set minimum sequence length.
    pub fn with_sequence_length_min(mut self, length: usize) -> Self {
        self.sequence_length_min = Some(length);
        self
    }

    /// Set maximum sequence length.
    pub fn with_sequence_length_max(mut self, length: usize) -> Self {
        self.sequence_length_max = Some(length);
        self
    }

    /// Filter by PDB ID prefix.
    pub fn with_pdb_id_prefix(mut self, prefix: &str) -> Self {
        self.pdb_id_prefix = Some(prefix.to_string());
        self
    }

    /// Filter by enzyme classification (EC) number.
    pub fn with_ec_number(mut self, ec: &str) -> Self {
        self.ec_number = Some(ec.to_string());
        self
    }

    /// Check if the query has any search criteria.
    pub fn is_empty(&self) -> bool {
        self.text.is_none()
            && self.organism.is_none()
            && self.resolution_max.is_none()
            && self.resolution_min.is_none()
            && self.experimental_method.is_none()
            && self.polymer_type.is_none()
            && self.release_date_min.is_none()
            && self.release_date_max.is_none()
            && self.sequence_length_min.is_none()
            && self.sequence_length_max.is_none()
            && self.pdb_id_prefix.is_none()
            && self.ec_number.is_none()
    }

    /// Convert the query to a JSON string for the RCSB API.
    pub fn to_json(&self) -> String {
        let query_value = self.build_query();
        serde_json::to_string(&query_value).unwrap_or_default()
    }

    /// Build the query JSON value.
    fn build_query(&self) -> Value {
        let mut nodes: Vec<Value> = Vec::new();

        // Full-text search
        if let Some(ref text) = self.text {
            nodes.push(json!({
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": text
                }
            }));
        }

        // Organism filter
        if let Some(ref organism) = self.organism {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                    "operator": "exact_match",
                    "value": organism
                }
            }));
        }

        // Resolution range
        if let Some(max) = self.resolution_max {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": max
                }
            }));
        }

        if let Some(min) = self.resolution_min {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "greater_or_equal",
                    "value": min
                }
            }));
        }

        // Experimental method
        if let Some(ref method) = self.experimental_method {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": method.api_value()
                }
            }));
        }

        // Polymer type
        if let Some(ref polymer) = self.polymer_type {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "value": polymer.api_value()
                }
            }));
        }

        // Release date range
        if let Some(ref date) = self.release_date_min {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_accession_info.initial_release_date",
                    "operator": "greater_or_equal",
                    "value": date
                }
            }));
        }

        if let Some(ref date) = self.release_date_max {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_accession_info.initial_release_date",
                    "operator": "less_or_equal",
                    "value": date
                }
            }));
        }

        // Sequence length range
        if let Some(min) = self.sequence_length_min {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_sample_sequence_length",
                    "operator": "greater_or_equal",
                    "value": min
                }
            }));
        }

        if let Some(max) = self.sequence_length_max {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_sample_sequence_length",
                    "operator": "less_or_equal",
                    "value": max
                }
            }));
        }

        // EC number
        if let Some(ref ec) = self.ec_number {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
                    "operator": "exact_match",
                    "value": ec
                }
            }));
        }

        // PDB ID prefix
        if let Some(ref prefix) = self.pdb_id_prefix {
            nodes.push(json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_container_identifiers.entry_id",
                    "operator": "contains_phrase",
                    "value": prefix
                }
            }));
        }

        // Build the final query
        let query = if nodes.is_empty() {
            // Default to returning all structures if no criteria specified
            json!({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "exists"
                }
            })
        } else if nodes.len() == 1 {
            nodes.into_iter().next().unwrap()
        } else {
            json!({
                "type": "group",
                "logical_operator": "and",
                "nodes": nodes
            })
        };

        json!({
            "query": query,
            "return_type": "entry",
            "request_options": {
                "return_all_hits": true,
                "results_content_type": ["experimental"]
            }
        })
    }
}

/// Result from an RCSB search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResult {
    /// List of PDB IDs matching the query
    pub pdb_ids: Vec<String>,
    /// Total number of results found
    pub total_count: usize,
}

/// Perform a search against the RCSB PDB.
///
/// # Arguments
///
/// * `query` - The search query
/// * `max_results` - Maximum number of results to return (0 for all)
///
/// # Returns
///
/// A `SearchResult` containing the matching PDB IDs.
///
/// # Errors
///
/// Returns a `SearchError` if the request fails or the response cannot be parsed.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{SearchQuery, rcsb_search};
///
/// let query = SearchQuery::new()
///     .with_text("ubiquitin")
///     .with_resolution_max(2.0);
///
/// let result = rcsb_search(&query, 10)?;
/// println!("Found {} structures", result.pdb_ids.len());
/// ```
pub fn rcsb_search(query: &SearchQuery, max_results: usize) -> Result<SearchResult, SearchError> {
    let client = reqwest::blocking::Client::new();

    let json_body = query.to_json();

    let response = client
        .post(RCSB_SEARCH_URL)
        .header("Content-Type", "application/json")
        .body(json_body)
        .send()?;

    if !response.status().is_success() {
        let status = response.status();
        let body = response.text().unwrap_or_default();
        return Err(SearchError::ApiError(format!(
            "HTTP {}: {}",
            status, body
        )));
    }

    let body: Value = response.json()?;

    // Parse the response
    let total_count = body
        .get("total_count")
        .and_then(|v| v.as_u64())
        .unwrap_or(0) as usize;

    let mut pdb_ids: Vec<String> = body
        .get("result_set")
        .and_then(|v| v.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|item| {
                    item.get("identifier")
                        .and_then(|v| v.as_str())
                        .map(|s| s.to_string())
                })
                .collect()
        })
        .unwrap_or_default();

    // Limit results if requested
    if max_results > 0 && pdb_ids.len() > max_results {
        pdb_ids.truncate(max_results);
    }

    Ok(SearchResult {
        pdb_ids,
        total_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_query_default() {
        let query = SearchQuery::new();
        assert!(query.is_empty());
    }

    #[test]
    fn test_search_query_with_text() {
        let query = SearchQuery::new().with_text("kinase");
        assert_eq!(query.text, Some("kinase".to_string()));
        assert!(!query.is_empty());
    }

    #[test]
    fn test_search_query_chain() {
        let query = SearchQuery::new()
            .with_text("ubiquitin")
            .with_organism("Homo sapiens")
            .with_resolution_max(2.0)
            .with_experimental_method(ExperimentalMethod::XRay);

        assert_eq!(query.text, Some("ubiquitin".to_string()));
        assert_eq!(query.organism, Some("Homo sapiens".to_string()));
        assert_eq!(query.resolution_max, Some(2.0));
        assert_eq!(query.experimental_method, Some(ExperimentalMethod::XRay));
    }

    #[test]
    fn test_search_query_to_json() {
        let query = SearchQuery::new().with_text("kinase");
        let json = query.to_json();

        assert!(json.contains("kinase"));
        assert!(json.contains("full_text"));
    }

    #[test]
    fn test_search_query_resolution_range() {
        let query = SearchQuery::new()
            .with_resolution_min(1.5)
            .with_resolution_max(2.5);

        let json = query.to_json();
        assert!(json.contains("resolution_combined"));
    }

    #[test]
    fn test_experimental_method_api_value() {
        assert_eq!(ExperimentalMethod::XRay.api_value(), "X-RAY DIFFRACTION");
        assert_eq!(ExperimentalMethod::Nmr.api_value(), "SOLUTION NMR");
        assert_eq!(ExperimentalMethod::Em.api_value(), "ELECTRON MICROSCOPY");
    }

    #[test]
    fn test_polymer_type_api_value() {
        assert_eq!(PolymerType::Protein.api_value(), "Protein");
        assert_eq!(PolymerType::Dna.api_value(), "DNA");
        assert_eq!(PolymerType::Rna.api_value(), "RNA");
    }

    #[test]
    fn test_search_error_display() {
        let err = SearchError::RequestFailed("timeout".to_string());
        assert_eq!(err.to_string(), "Request failed: timeout");

        let err = SearchError::ApiError("404".to_string());
        assert_eq!(err.to_string(), "API error: 404");
    }
}
