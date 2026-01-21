//! Async RCSB PDB file download utilities.
//!
//! This module provides asynchronous functionality for downloading structure files
//! from the RCSB Protein Data Bank with concurrency control.
//!
//! # Examples
//!
//! ## Download multiple structures concurrently
//!
//! ```ignore
//! use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};
//!
//! #[tokio::main]
//! async fn main() {
//!     let pdb_ids = vec!["1UBQ", "8HM2", "4INS", "1HHB", "2MBP"];
//!
//!     // Download with default options (5 concurrent, 100ms rate limit)
//!     let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;
//!
//!     for (pdb_id, result) in results {
//!         match result {
//!             Ok(structure) => println!("{}: {} atoms", pdb_id, structure.atoms.len()),
//!             Err(e) => eprintln!("{}: {}", pdb_id, e),
//!         }
//!     }
//! }
//! ```
//!
//! ## Download with custom options
//!
//! ```ignore
//! use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};
//!
//! #[tokio::main]
//! async fn main() {
//!     let options = AsyncDownloadOptions {
//!         max_concurrent: 10,
//!         rate_limit_ms: 50,
//!         ..Default::default()
//!     };
//!
//!     let results = download_multiple_async(&pdb_ids, FileFormat::Cif, Some(options)).await;
//! }
//! ```

use std::path::Path;
use std::sync::Arc;
use std::time::Duration;

use futures::future::join_all;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;
use tokio::sync::Semaphore;
use tokio::time::sleep;

use crate::core::PdbStructure;
use crate::parser::{parse_mmcif_string, parse_pdb_string};

use super::{DownloadError, FileFormat, RCSB_DOWNLOAD_URL};

/// Options for controlling async download behavior.
///
/// # Examples
///
/// ```
/// use pdbrust::rcsb::AsyncDownloadOptions;
///
/// // Default options
/// let opts = AsyncDownloadOptions::default();
/// assert_eq!(opts.max_concurrent, 5);
/// assert_eq!(opts.rate_limit_ms, 100);
///
/// // Conservative options for rate-limited APIs
/// let conservative = AsyncDownloadOptions::conservative();
/// assert_eq!(conservative.max_concurrent, 2);
///
/// // Fast options for when you have permission
/// let fast = AsyncDownloadOptions::fast();
/// assert_eq!(fast.max_concurrent, 20);
/// ```
#[derive(Debug, Clone)]
pub struct AsyncDownloadOptions {
    /// Maximum number of concurrent downloads (default: 5)
    pub max_concurrent: usize,
    /// Minimum delay between requests in milliseconds (default: 100)
    pub rate_limit_ms: u64,
    /// Request timeout in seconds (default: 30)
    pub timeout_secs: u64,
    /// Number of retry attempts on failure (default: 2)
    pub retries: usize,
}

impl Default for AsyncDownloadOptions {
    fn default() -> Self {
        Self {
            max_concurrent: 5,
            rate_limit_ms: 100,
            timeout_secs: 30,
            retries: 2,
        }
    }
}

impl AsyncDownloadOptions {
    /// Create new options with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create conservative options suitable for rate-limited APIs.
    ///
    /// Uses 2 concurrent connections and 500ms rate limiting.
    pub fn conservative() -> Self {
        Self {
            max_concurrent: 2,
            rate_limit_ms: 500,
            timeout_secs: 60,
            retries: 3,
        }
    }

    /// Create fast options for when you have appropriate API access.
    ///
    /// Uses 20 concurrent connections and 25ms rate limiting.
    /// Use responsibly and only when you have permission.
    pub fn fast() -> Self {
        Self {
            max_concurrent: 20,
            rate_limit_ms: 25,
            timeout_secs: 30,
            retries: 1,
        }
    }

    /// Set the maximum number of concurrent downloads.
    pub fn with_max_concurrent(mut self, max_concurrent: usize) -> Self {
        self.max_concurrent = max_concurrent;
        self
    }

    /// Set the rate limit delay in milliseconds.
    pub fn with_rate_limit_ms(mut self, rate_limit_ms: u64) -> Self {
        self.rate_limit_ms = rate_limit_ms;
        self
    }

    /// Set the request timeout in seconds.
    pub fn with_timeout_secs(mut self, timeout_secs: u64) -> Self {
        self.timeout_secs = timeout_secs;
        self
    }

    /// Set the number of retry attempts.
    pub fn with_retries(mut self, retries: usize) -> Self {
        self.retries = retries;
        self
    }
}

/// Build the download URL for a PDB structure.
fn build_download_url(pdb_id: &str, format: FileFormat) -> String {
    let pdb_id_upper = pdb_id.to_uppercase();
    match format {
        FileFormat::Pdb => format!("{}/{}.pdb", RCSB_DOWNLOAD_URL, pdb_id_upper),
        FileFormat::Cif => format!("{}/{}.cif", RCSB_DOWNLOAD_URL, pdb_id_upper),
    }
}

/// Download a PDB file as a string asynchronously.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID
/// * `format` - The desired file format
///
/// # Returns
///
/// The file contents as a string.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_pdb_string_async, FileFormat};
///
/// #[tokio::main]
/// async fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let content = download_pdb_string_async("1UBQ", FileFormat::Pdb).await?;
///     println!("Downloaded {} bytes", content.len());
///     Ok(())
/// }
/// ```
pub async fn download_pdb_string_async(
    pdb_id: &str,
    format: FileFormat,
) -> Result<String, DownloadError> {
    let url = build_download_url(pdb_id, format);
    let client = reqwest::Client::new();

    let response = client.get(&url).send().await?;

    if response.status() == reqwest::StatusCode::NOT_FOUND {
        return Err(DownloadError::NotFound(pdb_id.to_string()));
    }

    if !response.status().is_success() {
        return Err(DownloadError::RequestFailed(format!(
            "HTTP {}: {}",
            response.status(),
            response.text().await.unwrap_or_default()
        )));
    }

    Ok(response.text().await?)
}

/// Download and parse a structure from RCSB PDB asynchronously.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID
/// * `format` - The desired file format
///
/// # Returns
///
/// A parsed `PdbStructure`.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_structure_async, FileFormat};
///
/// #[tokio::main]
/// async fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let structure = download_structure_async("1UBQ", FileFormat::Pdb).await?;
///     println!("Downloaded {} atoms", structure.atoms.len());
///     Ok(())
/// }
/// ```
pub async fn download_structure_async(
    pdb_id: &str,
    format: FileFormat,
) -> Result<PdbStructure, DownloadError> {
    let content = download_pdb_string_async(pdb_id, format).await?;

    let structure = match format {
        FileFormat::Pdb => parse_pdb_string(&content)?,
        FileFormat::Cif => parse_mmcif_string(&content)?,
    };

    Ok(structure)
}

/// Download a structure to a file asynchronously.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID
/// * `path` - The path to save the file
/// * `format` - The desired file format
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_to_file_async, FileFormat};
///
/// #[tokio::main]
/// async fn main() -> Result<(), Box<dyn std::error::Error>> {
///     download_to_file_async("1UBQ", "1UBQ.pdb", FileFormat::Pdb).await?;
///     Ok(())
/// }
/// ```
pub async fn download_to_file_async<P: AsRef<Path>>(
    pdb_id: &str,
    path: P,
    format: FileFormat,
) -> Result<(), DownloadError> {
    let content = download_pdb_string_async(pdb_id, format).await?;
    let mut file = File::create(path).await?;
    file.write_all(content.as_bytes()).await?;
    Ok(())
}

/// Download a single structure with retry logic.
async fn download_with_retry(
    client: &reqwest::Client,
    pdb_id: &str,
    format: FileFormat,
    retries: usize,
    timeout: Duration,
) -> Result<PdbStructure, DownloadError> {
    let url = build_download_url(pdb_id, format);
    let mut last_error = DownloadError::RequestFailed("No attempts made".to_string());

    for attempt in 0..=retries {
        if attempt > 0 {
            // Exponential backoff: 1s, 2s, 4s, ...
            let backoff = Duration::from_secs(1 << (attempt - 1));
            sleep(backoff).await;
        }

        let result = async {
            let response = client
                .get(&url)
                .timeout(timeout)
                .send()
                .await?;

            if response.status() == reqwest::StatusCode::NOT_FOUND {
                return Err(DownloadError::NotFound(pdb_id.to_string()));
            }

            if !response.status().is_success() {
                return Err(DownloadError::RequestFailed(format!(
                    "HTTP {}: {}",
                    response.status(),
                    response.text().await.unwrap_or_default()
                )));
            }

            let content = response.text().await?;

            let structure = match format {
                FileFormat::Pdb => parse_pdb_string(&content)?,
                FileFormat::Cif => parse_mmcif_string(&content)?,
            };

            Ok(structure)
        }
        .await;

        match result {
            Ok(structure) => return Ok(structure),
            Err(e) => {
                // Don't retry on 404
                if matches!(e, DownloadError::NotFound(_)) {
                    return Err(e);
                }
                last_error = e;
            }
        }
    }

    Err(last_error)
}

/// Download multiple structures concurrently with controlled parallelism.
///
/// This is the main function for bulk downloading structures from RCSB PDB.
/// It uses a semaphore to limit concurrent connections and implements
/// rate limiting between requests.
///
/// # Arguments
///
/// * `pdb_ids` - A slice of PDB IDs to download
/// * `format` - The desired file format
/// * `options` - Optional download configuration (defaults if None)
///
/// # Returns
///
/// A vector of (pdb_id, result) tuples, maintaining the input order.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};
///
/// #[tokio::main]
/// async fn main() {
///     let pdb_ids = vec!["1UBQ", "8HM2", "4INS"];
///
///     // Simple usage with default options
///     let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;
///
///     let successful: Vec<_> = results.iter()
///         .filter(|(_, r)| r.is_ok())
///         .collect();
///     println!("Downloaded {} of {} structures", successful.len(), pdb_ids.len());
///
///     // With custom options
///     let options = AsyncDownloadOptions::default()
///         .with_max_concurrent(10)
///         .with_rate_limit_ms(50);
///     let results = download_multiple_async(&pdb_ids, FileFormat::Cif, Some(options)).await;
/// }
/// ```
pub async fn download_multiple_async(
    pdb_ids: &[&str],
    format: FileFormat,
    options: Option<AsyncDownloadOptions>,
) -> Vec<(String, Result<PdbStructure, DownloadError>)> {
    let options = options.unwrap_or_default();
    let semaphore = Arc::new(Semaphore::new(options.max_concurrent));
    let client = reqwest::Client::builder()
        .timeout(Duration::from_secs(options.timeout_secs))
        .build()
        .unwrap_or_else(|_| reqwest::Client::new());
    let client = Arc::new(client);
    let rate_limit = Duration::from_millis(options.rate_limit_ms);
    let timeout = Duration::from_secs(options.timeout_secs);
    let retries = options.retries;

    let tasks: Vec<_> = pdb_ids
        .iter()
        .enumerate()
        .map(|(idx, &pdb_id)| {
            let semaphore = Arc::clone(&semaphore);
            let client = Arc::clone(&client);
            let pdb_id = pdb_id.to_string();

            async move {
                // Wait for a permit from the semaphore
                let _permit = semaphore.acquire().await.unwrap();

                // Rate limiting: add delay based on index
                if idx > 0 {
                    sleep(rate_limit).await;
                }

                let result = download_with_retry(&client, &pdb_id, format, retries, timeout).await;

                (pdb_id, result)
            }
        })
        .collect();

    join_all(tasks).await
}

/// Download multiple structures to files concurrently.
///
/// # Arguments
///
/// * `pdb_ids` - A slice of PDB IDs to download
/// * `output_dir` - The directory to save files in
/// * `format` - The desired file format
/// * `options` - Optional download configuration (defaults if None)
///
/// # Returns
///
/// A vector of (pdb_id, result) tuples indicating success or failure.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_multiple_to_files_async, AsyncDownloadOptions, FileFormat};
///
/// #[tokio::main]
/// async fn main() {
///     let pdb_ids = vec!["1UBQ", "8HM2", "4INS"];
///
///     let results = download_multiple_to_files_async(
///         &pdb_ids,
///         "./structures",
///         FileFormat::Pdb,
///         None,
///     ).await;
///
///     for (pdb_id, result) in results {
///         match result {
///             Ok(()) => println!("{}: saved", pdb_id),
///             Err(e) => eprintln!("{}: {}", pdb_id, e),
///         }
///     }
/// }
/// ```
pub async fn download_multiple_to_files_async<P: AsRef<Path> + Sync>(
    pdb_ids: &[&str],
    output_dir: P,
    format: FileFormat,
    options: Option<AsyncDownloadOptions>,
) -> Vec<(String, Result<(), DownloadError>)> {
    let options = options.unwrap_or_default();
    let semaphore = Arc::new(Semaphore::new(options.max_concurrent));
    let client = reqwest::Client::builder()
        .timeout(Duration::from_secs(options.timeout_secs))
        .build()
        .unwrap_or_else(|_| reqwest::Client::new());
    let client = Arc::new(client);
    let rate_limit = Duration::from_millis(options.rate_limit_ms);
    let timeout = Duration::from_secs(options.timeout_secs);
    let retries = options.retries;
    let output_dir = output_dir.as_ref().to_path_buf();

    let tasks: Vec<_> = pdb_ids
        .iter()
        .enumerate()
        .map(|(idx, &pdb_id)| {
            let semaphore = Arc::clone(&semaphore);
            let client = Arc::clone(&client);
            let pdb_id = pdb_id.to_string();
            let output_dir = output_dir.clone();

            async move {
                // Wait for a permit from the semaphore
                let _permit = semaphore.acquire().await.unwrap();

                // Rate limiting
                if idx > 0 {
                    sleep(rate_limit).await;
                }

                let url = build_download_url(&pdb_id, format);
                let filename = format!("{}.{}", pdb_id.to_uppercase(), format.extension());
                let path = output_dir.join(filename);

                let result: Result<(), DownloadError> = async {
                    let mut last_error =
                        DownloadError::RequestFailed("No attempts made".to_string());

                    for attempt in 0..=retries {
                        if attempt > 0 {
                            let backoff = Duration::from_secs(1 << (attempt - 1));
                            sleep(backoff).await;
                        }

                        let download_result = async {
                            let response = client.get(&url).timeout(timeout).send().await?;

                            if response.status() == reqwest::StatusCode::NOT_FOUND {
                                return Err(DownloadError::NotFound(pdb_id.clone()));
                            }

                            if !response.status().is_success() {
                                return Err(DownloadError::RequestFailed(format!(
                                    "HTTP {}: {}",
                                    response.status(),
                                    response.text().await.unwrap_or_default()
                                )));
                            }

                            let content = response.text().await?;
                            let mut file = File::create(&path).await?;
                            file.write_all(content.as_bytes()).await?;
                            Ok(())
                        }
                        .await;

                        match download_result {
                            Ok(()) => return Ok(()),
                            Err(e) => {
                                if matches!(e, DownloadError::NotFound(_)) {
                                    return Err(e);
                                }
                                last_error = e;
                            }
                        }
                    }

                    Err(last_error)
                }
                .await;

                (pdb_id, result)
            }
        })
        .collect();

    join_all(tasks).await
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_async_options_default() {
        let options = AsyncDownloadOptions::default();
        assert_eq!(options.max_concurrent, 5);
        assert_eq!(options.rate_limit_ms, 100);
        assert_eq!(options.timeout_secs, 30);
        assert_eq!(options.retries, 2);
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
    fn test_async_options_builder() {
        let options = AsyncDownloadOptions::new()
            .with_max_concurrent(10)
            .with_rate_limit_ms(50)
            .with_timeout_secs(45)
            .with_retries(5);

        assert_eq!(options.max_concurrent, 10);
        assert_eq!(options.rate_limit_ms, 50);
        assert_eq!(options.timeout_secs, 45);
        assert_eq!(options.retries, 5);
    }

    #[test]
    fn test_build_download_url_pdb() {
        let url = build_download_url("1ubq", FileFormat::Pdb);
        assert_eq!(url, "https://files.rcsb.org/download/1UBQ.pdb");
    }

    #[test]
    fn test_build_download_url_cif() {
        let url = build_download_url("8hm2", FileFormat::Cif);
        assert_eq!(url, "https://files.rcsb.org/download/8HM2.cif");
    }
}
