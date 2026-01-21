//! Python bindings for RCSB PDB search and download

use crate::error::{convert_download_error, convert_search_error};
use crate::structure::PyPdbStructure;
#[cfg(feature = "rcsb-async")]
use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions};
use pdbrust::rcsb::{
    download_structure as rust_download_structure, rcsb_search as rust_rcsb_search,
    ExperimentalMethod, FileFormat, PolymerType, SearchQuery, SearchResult,
};
use pyo3::prelude::*;

/// File format for PDB downloads
#[pyclass(name = "FileFormat")]
#[derive(Clone)]
pub struct PyFileFormat {
    inner: FileFormat,
}

#[pymethods]
impl PyFileFormat {
    /// PDB format
    #[staticmethod]
    fn pdb() -> Self {
        PyFileFormat {
            inner: FileFormat::Pdb,
        }
    }

    /// mmCIF format
    #[staticmethod]
    fn cif() -> Self {
        PyFileFormat {
            inner: FileFormat::Cif,
        }
    }

    /// Get file extension
    fn extension(&self) -> &'static str {
        self.inner.extension()
    }

    /// Get compressed file extension
    fn compressed_extension(&self) -> &'static str {
        self.inner.compressed_extension()
    }

    fn __repr__(&self) -> String {
        match self.inner {
            FileFormat::Pdb => "FileFormat.pdb()".to_string(),
            FileFormat::Cif => "FileFormat.cif()".to_string(),
        }
    }
}

/// Experimental method filter for RCSB search
#[pyclass(name = "ExperimentalMethod")]
#[derive(Clone)]
pub struct PyExperimentalMethod {
    inner: ExperimentalMethod,
}

#[pymethods]
impl PyExperimentalMethod {
    #[staticmethod]
    fn xray() -> Self {
        PyExperimentalMethod {
            inner: ExperimentalMethod::XRay,
        }
    }

    #[staticmethod]
    fn nmr() -> Self {
        PyExperimentalMethod {
            inner: ExperimentalMethod::Nmr,
        }
    }

    #[staticmethod]
    fn em() -> Self {
        PyExperimentalMethod {
            inner: ExperimentalMethod::Em,
        }
    }

    #[staticmethod]
    fn other() -> Self {
        PyExperimentalMethod {
            inner: ExperimentalMethod::Other,
        }
    }

    fn __repr__(&self) -> String {
        format!("ExperimentalMethod.{:?}", self.inner)
    }
}

/// Polymer type filter for RCSB search
#[pyclass(name = "PolymerType")]
#[derive(Clone)]
pub struct PyPolymerType {
    inner: PolymerType,
}

#[pymethods]
impl PyPolymerType {
    #[staticmethod]
    fn protein() -> Self {
        PyPolymerType {
            inner: PolymerType::Protein,
        }
    }

    #[staticmethod]
    fn dna() -> Self {
        PyPolymerType {
            inner: PolymerType::Dna,
        }
    }

    #[staticmethod]
    fn rna() -> Self {
        PyPolymerType {
            inner: PolymerType::Rna,
        }
    }

    fn __repr__(&self) -> String {
        format!("PolymerType.{:?}", self.inner)
    }
}

/// Search query builder for RCSB PDB
#[pyclass(name = "SearchQuery")]
#[derive(Clone)]
pub struct PySearchQuery {
    inner: SearchQuery,
}

#[pymethods]
impl PySearchQuery {
    /// Create a new empty search query
    #[new]
    fn new() -> Self {
        PySearchQuery {
            inner: SearchQuery::new(),
        }
    }

    /// Add text search term
    fn with_text(&self, text: &str) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_text(text),
        }
    }

    /// Filter by organism
    fn with_organism(&self, organism: &str) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_organism(organism),
        }
    }

    /// Set maximum resolution (Angstroms)
    fn with_resolution_max(&self, resolution: f64) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_resolution_max(resolution),
        }
    }

    /// Set minimum resolution (Angstroms)
    fn with_resolution_min(&self, resolution: f64) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_resolution_min(resolution),
        }
    }

    /// Filter by experimental method
    fn with_experimental_method(&self, method: &PyExperimentalMethod) -> PySearchQuery {
        PySearchQuery {
            inner: self
                .inner
                .clone()
                .with_experimental_method(method.inner.clone()),
        }
    }

    /// Filter by polymer type
    fn with_polymer_type(&self, polymer_type: &PyPolymerType) -> PySearchQuery {
        PySearchQuery {
            inner: self
                .inner
                .clone()
                .with_polymer_type(polymer_type.inner.clone()),
        }
    }

    /// Set minimum release date (YYYY-MM-DD format)
    fn with_release_date_min(&self, date: &str) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_release_date_min(date),
        }
    }

    /// Set maximum release date (YYYY-MM-DD format)
    fn with_release_date_max(&self, date: &str) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_release_date_max(date),
        }
    }

    /// Set minimum sequence length
    fn with_sequence_length_min(&self, length: usize) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_sequence_length_min(length),
        }
    }

    /// Set maximum sequence length
    fn with_sequence_length_max(&self, length: usize) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_sequence_length_max(length),
        }
    }

    /// Filter by EC number
    fn with_ec_number(&self, ec: &str) -> PySearchQuery {
        PySearchQuery {
            inner: self.inner.clone().with_ec_number(ec),
        }
    }

    /// Check if query is empty
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn __repr__(&self) -> String {
        format!("SearchQuery(text={:?})", self.inner.text)
    }
}

/// Search result from RCSB PDB
#[pyclass(name = "SearchResult")]
pub struct PySearchResult {
    /// List of PDB IDs matching the query
    #[pyo3(get)]
    pdb_ids: Vec<String>,
    /// Total number of matches (may be more than returned)
    #[pyo3(get)]
    total_count: usize,
}

#[pymethods]
impl PySearchResult {
    fn __repr__(&self) -> String {
        format!(
            "SearchResult(count={}, total={})",
            self.pdb_ids.len(),
            self.total_count
        )
    }

    fn __len__(&self) -> usize {
        self.pdb_ids.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<SearchResultIter>> {
        let iter = SearchResultIter {
            inner: slf.pdb_ids.clone().into_iter(),
        };
        Py::new(slf.py(), iter)
    }
}

#[pyclass]
struct SearchResultIter {
    inner: std::vec::IntoIter<String>,
}

#[pymethods]
impl SearchResultIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<String> {
        slf.inner.next()
    }
}

impl From<SearchResult> for PySearchResult {
    fn from(result: SearchResult) -> Self {
        PySearchResult {
            pdb_ids: result.pdb_ids,
            total_count: result.total_count,
        }
    }
}

/// Search RCSB PDB
///
/// Args:
///     query: SearchQuery object with search criteria
///     max_results: Maximum number of results to return
///
/// Returns:
///     SearchResult with list of PDB IDs
///
/// Raises:
///     RuntimeError: If search request fails
#[pyfunction]
pub fn rcsb_search(query: &PySearchQuery, max_results: usize) -> PyResult<PySearchResult> {
    rust_rcsb_search(&query.inner, max_results)
        .map(PySearchResult::from)
        .map_err(convert_search_error)
}

/// Download a structure from RCSB PDB
///
/// Args:
///     pdb_id: 4-character PDB ID (e.g., "1UBQ")
///     format: FileFormat.pdb() or FileFormat.cif()
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     ValueError: If PDB ID not found
///     RuntimeError: If download fails
#[pyfunction]
pub fn download_structure(pdb_id: &str, format: &PyFileFormat) -> PyResult<PyPdbStructure> {
    rust_download_structure(pdb_id, format.inner)
        .map(PyPdbStructure::from)
        .map_err(convert_download_error)
}

/// Download structure content as string
///
/// Args:
///     pdb_id: 4-character PDB ID
///     format: FileFormat.pdb() or FileFormat.cif()
///
/// Returns:
///     String content of the structure file
#[pyfunction]
pub fn download_pdb_string(pdb_id: &str, format: &PyFileFormat) -> PyResult<String> {
    pdbrust::rcsb::download_pdb_string(pdb_id, format.inner).map_err(convert_download_error)
}

/// Download structure to a file
///
/// Args:
///     pdb_id: 4-character PDB ID
///     path: Output file path
///     format: FileFormat.pdb() or FileFormat.cif()
#[pyfunction]
pub fn download_to_file(pdb_id: &str, path: &str, format: &PyFileFormat) -> PyResult<()> {
    pdbrust::rcsb::download_to_file(pdb_id, path, format.inner).map_err(convert_download_error)
}

/// Options for controlling bulk download behavior
///
/// Example:
///     >>> options = AsyncDownloadOptions(max_concurrent=10, rate_limit_ms=50)
///     >>> results = download_multiple(pdb_ids, FileFormat.pdb(), options)
#[cfg(feature = "rcsb-async")]
#[pyclass(name = "AsyncDownloadOptions")]
#[derive(Clone)]
pub struct PyAsyncDownloadOptions {
    inner: AsyncDownloadOptions,
}

#[cfg(feature = "rcsb-async")]
#[pymethods]
impl PyAsyncDownloadOptions {
    /// Create download options
    ///
    /// Args:
    ///     max_concurrent: Maximum concurrent downloads (default: 5)
    ///     rate_limit_ms: Minimum delay between requests in ms (default: 100)
    ///     timeout_secs: Request timeout in seconds (default: 30)
    ///     retries: Number of retry attempts on failure (default: 2)
    #[new]
    #[pyo3(signature = (max_concurrent=5, rate_limit_ms=100, timeout_secs=30, retries=2))]
    fn new(max_concurrent: usize, rate_limit_ms: u64, timeout_secs: u64, retries: usize) -> Self {
        PyAsyncDownloadOptions {
            inner: AsyncDownloadOptions {
                max_concurrent,
                rate_limit_ms,
                timeout_secs,
                retries,
            },
        }
    }

    /// Create conservative options (2 concurrent, 500ms delay)
    #[staticmethod]
    fn conservative() -> Self {
        PyAsyncDownloadOptions {
            inner: AsyncDownloadOptions::conservative(),
        }
    }

    /// Create fast options (20 concurrent, 25ms delay)
    #[staticmethod]
    fn fast() -> Self {
        PyAsyncDownloadOptions {
            inner: AsyncDownloadOptions::fast(),
        }
    }

    /// Maximum concurrent downloads
    #[getter]
    fn max_concurrent(&self) -> usize {
        self.inner.max_concurrent
    }

    /// Rate limit in milliseconds
    #[getter]
    fn rate_limit_ms(&self) -> u64 {
        self.inner.rate_limit_ms
    }

    /// Timeout in seconds
    #[getter]
    fn timeout_secs(&self) -> u64 {
        self.inner.timeout_secs
    }

    /// Number of retries
    #[getter]
    fn retries(&self) -> usize {
        self.inner.retries
    }

    fn __repr__(&self) -> String {
        format!(
            "AsyncDownloadOptions(max_concurrent={}, rate_limit_ms={}, timeout_secs={}, retries={})",
            self.inner.max_concurrent,
            self.inner.rate_limit_ms,
            self.inner.timeout_secs,
            self.inner.retries
        )
    }
}

/// Download result for a single structure
#[cfg(feature = "rcsb-async")]
#[pyclass(name = "DownloadResult")]
pub struct PyDownloadResult {
    /// PDB ID
    #[pyo3(get)]
    pub pdb_id: String,
    /// Whether download succeeded
    #[pyo3(get)]
    pub success: bool,
    /// Structure (if successful) - stored as inner PdbStructure for cloning
    structure: Option<pdbrust::PdbStructure>,
    /// Error message (if failed)
    #[pyo3(get)]
    pub error: Option<String>,
}

#[cfg(feature = "rcsb-async")]
#[pymethods]
impl PyDownloadResult {
    /// Get the structure (raises if download failed)
    fn get_structure(&self) -> PyResult<PyPdbStructure> {
        self.structure
            .as_ref()
            .map(|s| PyPdbStructure::from(s.clone()))
            .ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "Download failed for {}: {}",
                    self.pdb_id,
                    self.error.as_deref().unwrap_or("Unknown error")
                ))
            })
    }

    fn __repr__(&self) -> String {
        if self.success {
            format!("DownloadResult({}: success)", self.pdb_id)
        } else {
            format!(
                "DownloadResult({}: failed - {})",
                self.pdb_id,
                self.error.as_deref().unwrap_or("Unknown error")
            )
        }
    }
}

/// Download multiple structures concurrently
///
/// This function efficiently downloads multiple PDB structures using
/// concurrent connections with rate limiting. The GIL is released
/// during I/O operations for better performance.
///
/// Args:
///     pdb_ids: List of 4-character PDB IDs to download
///     format: FileFormat.pdb() or FileFormat.cif()
///     options: Optional AsyncDownloadOptions for controlling concurrency
///
/// Returns:
///     List of DownloadResult objects, one for each PDB ID
///
/// Example:
///     >>> results = download_multiple(["1UBQ", "8HM2", "4INS"], FileFormat.pdb())
///     >>> for r in results:
///     ...     if r.success:
///     ...         print(f"{r.pdb_id}: {len(r.get_structure().atoms)} atoms")
///     ...     else:
///     ...         print(f"{r.pdb_id}: {r.error}")
#[cfg(feature = "rcsb-async")]
#[pyfunction]
#[pyo3(signature = (pdb_ids, format, options=None))]
pub fn download_multiple(
    py: Python<'_>,
    pdb_ids: Vec<String>,
    format: &PyFileFormat,
    options: Option<&PyAsyncDownloadOptions>,
) -> PyResult<Vec<PyDownloadResult>> {
    // Convert Python strings to &str slice
    let pdb_id_refs: Vec<&str> = pdb_ids.iter().map(|s| s.as_str()).collect();
    let rust_options = options.map(|o| o.inner.clone());
    let format = format.inner;

    // Release GIL during async I/O
    py.allow_threads(|| {
        // Create a new tokio runtime for this call
        let runtime = tokio::runtime::Runtime::new()
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let results = runtime
            .block_on(async { download_multiple_async(&pdb_id_refs, format, rust_options).await });

        // Convert results to Python types
        let py_results: Vec<PyDownloadResult> = results
            .into_iter()
            .map(|(pdb_id, result)| match result {
                Ok(structure) => PyDownloadResult {
                    pdb_id,
                    success: true,
                    structure: Some(structure),
                    error: None,
                },
                Err(e) => PyDownloadResult {
                    pdb_id,
                    success: false,
                    structure: None,
                    error: Some(e.to_string()),
                },
            })
            .collect();

        Ok(py_results)
    })
}
