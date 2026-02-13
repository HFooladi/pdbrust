//! Python bindings for DockQ v2 interface quality assessment.
//!
//! Provides DockQ scoring for evaluating protein-protein docking quality.

use pdbrust::dockq::{
    ChainMappingStrategy, DockQOptions, DockQQuality, DockQResult, InterfaceResult,
};
use pyo3::prelude::*;

/// Quality classification for a DockQ score.
///
/// DockQ scores are classified into four categories:
///
/// - Incorrect: DockQ < 0.23
/// - Acceptable: 0.23 <= DockQ < 0.49
/// - Medium: 0.49 <= DockQ < 0.80
/// - High: DockQ >= 0.80
///
/// Example:
///     >>> quality = DockQQuality.from_score(0.85)
///     >>> print(quality)  # "High"
///     >>> quality == DockQQuality.HIGH  # True
#[pyclass(name = "DockQQuality")]
#[derive(Clone)]
pub struct PyDockQQuality {
    inner: DockQQuality,
}

#[pymethods]
impl PyDockQQuality {
    /// DockQ < 0.23
    #[classattr]
    fn INCORRECT() -> Self {
        PyDockQQuality {
            inner: DockQQuality::Incorrect,
        }
    }

    /// 0.23 <= DockQ < 0.49
    #[classattr]
    fn ACCEPTABLE() -> Self {
        PyDockQQuality {
            inner: DockQQuality::Acceptable,
        }
    }

    /// 0.49 <= DockQ < 0.80
    #[classattr]
    fn MEDIUM() -> Self {
        PyDockQQuality {
            inner: DockQQuality::Medium,
        }
    }

    /// DockQ >= 0.80
    #[classattr]
    fn HIGH() -> Self {
        PyDockQQuality {
            inner: DockQQuality::High,
        }
    }

    /// Classify a DockQ score into a quality category.
    ///
    /// Args:
    ///     dockq: DockQ score (0.0 to 1.0)
    ///
    /// Returns:
    ///     DockQQuality classification
    #[staticmethod]
    fn from_score(dockq: f64) -> Self {
        PyDockQQuality {
            inner: DockQQuality::from_score(dockq),
        }
    }

    fn __repr__(&self) -> String {
        format!("DockQQuality.{}", self.__str__())
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        match self.inner {
            DockQQuality::Incorrect => 0,
            DockQQuality::Acceptable => 1,
            DockQQuality::Medium => 2,
            DockQQuality::High => 3,
        }
    }
}

impl From<DockQQuality> for PyDockQQuality {
    fn from(inner: DockQQuality) -> Self {
        Self { inner }
    }
}

/// Result for a single chain-pair interface.
///
/// Contains all DockQ score components for one interface:
/// fnat, fnonnat, iRMSD, LRMSD, and the combined DockQ score.
///
/// Attributes:
///     native_receptor_chain (str): Receptor chain ID in native
///     native_ligand_chain (str): Ligand chain ID in native
///     model_receptor_chain (str): Receptor chain ID in model
///     model_ligand_chain (str): Ligand chain ID in model
///     fnat (float): Fraction of native contacts preserved [0, 1]
///     fnonnat (float): Fraction of non-native contacts in model
///     irmsd (float): Interface RMSD in Angstroms
///     lrmsd (float): Ligand RMSD in Angstroms
///     dockq (float): Combined DockQ score [0, 1]
///     quality (DockQQuality): Quality classification
///     f1 (float): F1 score for contacts
///     num_native_contacts (int): Number of contacts in native
///     num_model_contacts (int): Number of contacts in model
///     num_clashes (int): Number of atomic clashes (< 2.0 A)
#[pyclass(name = "InterfaceResult")]
#[derive(Clone)]
pub struct PyInterfaceResult {
    inner: InterfaceResult,
}

#[pymethods]
impl PyInterfaceResult {
    /// Receptor chain ID in the native structure.
    #[getter]
    fn native_receptor_chain(&self) -> &str {
        &self.inner.native_chains.0
    }

    /// Ligand chain ID in the native structure.
    #[getter]
    fn native_ligand_chain(&self) -> &str {
        &self.inner.native_chains.1
    }

    /// Receptor chain ID in the model structure.
    #[getter]
    fn model_receptor_chain(&self) -> &str {
        &self.inner.model_chains.0
    }

    /// Ligand chain ID in the model structure.
    #[getter]
    fn model_ligand_chain(&self) -> &str {
        &self.inner.model_chains.1
    }

    /// Fraction of native contacts preserved [0, 1].
    #[getter]
    fn fnat(&self) -> f64 {
        self.inner.fnat
    }

    /// Fraction of non-native contacts in the model.
    #[getter]
    fn fnonnat(&self) -> f64 {
        self.inner.fnonnat
    }

    /// Interface RMSD in Angstroms.
    #[getter]
    fn irmsd(&self) -> f64 {
        self.inner.irmsd
    }

    /// Ligand RMSD in Angstroms.
    #[getter]
    fn lrmsd(&self) -> f64 {
        self.inner.lrmsd
    }

    /// Combined DockQ score [0, 1].
    #[getter]
    fn dockq(&self) -> f64 {
        self.inner.dockq
    }

    /// Quality classification (Incorrect/Acceptable/Medium/High).
    #[getter]
    fn quality(&self) -> PyDockQQuality {
        PyDockQQuality::from(self.inner.quality)
    }

    /// F1 score for contacts.
    #[getter]
    fn f1(&self) -> f64 {
        self.inner.f1
    }

    /// Number of contacts in the native structure.
    #[getter]
    fn num_native_contacts(&self) -> usize {
        self.inner.num_native_contacts
    }

    /// Number of contacts in the model structure.
    #[getter]
    fn num_model_contacts(&self) -> usize {
        self.inner.num_model_contacts
    }

    /// Number of atomic clashes (< 2.0 A).
    #[getter]
    fn num_clashes(&self) -> usize {
        self.inner.num_clashes
    }

    fn __repr__(&self) -> String {
        format!(
            "InterfaceResult(native={}-{}, model={}-{}, DockQ={:.4}, fnat={:.4}, iRMSD={:.2}, LRMSD={:.2}, {})",
            self.inner.native_chains.0,
            self.inner.native_chains.1,
            self.inner.model_chains.0,
            self.inner.model_chains.1,
            self.inner.dockq,
            self.inner.fnat,
            self.inner.irmsd,
            self.inner.lrmsd,
            self.inner.quality,
        )
    }
}

impl From<InterfaceResult> for PyInterfaceResult {
    fn from(inner: InterfaceResult) -> Self {
        Self { inner }
    }
}

/// Result of DockQ calculation for an entire complex.
///
/// Contains per-interface results and the overall DockQ score
/// (average over all interfaces with equal weighting).
///
/// Attributes:
///     interfaces (list[InterfaceResult]): Per-interface results
///     total_dockq (float): Average DockQ over all interfaces
///     chain_mapping (list[tuple[str, str]]): Model-to-native chain mapping
///     num_interfaces (int): Number of interfaces evaluated
///
/// Example:
///     >>> result = model.dockq_to(native)
///     >>> print(f"DockQ: {result.total_dockq:.4f}")
///     >>> for iface in result:
///     ...     print(f"  {iface.native_receptor_chain}-{iface.native_ligand_chain}: {iface.dockq:.3f}")
#[pyclass(name = "DockQResult")]
#[derive(Clone)]
pub struct PyDockQResult {
    inner: DockQResult,
}

#[pymethods]
impl PyDockQResult {
    /// Per-interface results.
    #[getter]
    fn interfaces(&self) -> Vec<PyInterfaceResult> {
        self.inner
            .interfaces
            .iter()
            .cloned()
            .map(PyInterfaceResult::from)
            .collect()
    }

    /// Average DockQ over all interfaces (equal weighting).
    #[getter]
    fn total_dockq(&self) -> f64 {
        self.inner.total_dockq
    }

    /// Model-to-native chain mapping as list of (model_chain, native_chain) tuples.
    #[getter]
    fn chain_mapping(&self) -> Vec<(String, String)> {
        self.inner.chain_mapping.clone()
    }

    /// Number of interfaces evaluated.
    #[getter]
    fn num_interfaces(&self) -> usize {
        self.inner.num_interfaces
    }

    fn __repr__(&self) -> String {
        format!(
            "DockQResult(total_dockq={:.4}, num_interfaces={})",
            self.inner.total_dockq, self.inner.num_interfaces
        )
    }

    fn __str__(&self) -> String {
        format!(
            "DockQ: {:.4} ({} interface{})",
            self.inner.total_dockq,
            self.inner.num_interfaces,
            if self.inner.num_interfaces != 1 {
                "s"
            } else {
                ""
            }
        )
    }

    fn __len__(&self) -> usize {
        self.inner.interfaces.len()
    }

    fn __getitem__(&self, idx: usize) -> PyResult<PyInterfaceResult> {
        self.inner
            .interfaces
            .get(idx)
            .cloned()
            .map(PyInterfaceResult::from)
            .ok_or_else(|| {
                pyo3::exceptions::PyIndexError::new_err(format!(
                    "index {} out of range for {} interfaces",
                    idx,
                    self.inner.interfaces.len()
                ))
            })
    }
}

impl From<DockQResult> for PyDockQResult {
    fn from(inner: DockQResult) -> Self {
        Self { inner }
    }
}

/// Strategy for establishing chain correspondence between model and native.
///
/// Use ``ChainMappingStrategy.auto()`` for automatic sequence-based mapping,
/// or ``ChainMappingStrategy.explicit(mapping)`` for manual chain assignment.
///
/// Example:
///     >>> # Automatic mapping (default)
///     >>> strategy = ChainMappingStrategy.auto()
///     >>>
///     >>> # Explicit mapping
///     >>> strategy = ChainMappingStrategy.explicit([("A", "A"), ("B", "B")])
#[pyclass(name = "ChainMappingStrategy")]
#[derive(Clone)]
pub struct PyChainMappingStrategy {
    pub(crate) inner: ChainMappingStrategy,
}

#[pymethods]
impl PyChainMappingStrategy {
    /// Create an automatic chain mapping strategy.
    ///
    /// Uses sequence alignment to find the best chain correspondence.
    ///
    /// Returns:
    ///     ChainMappingStrategy: Automatic mapping strategy
    #[staticmethod]
    fn auto() -> Self {
        PyChainMappingStrategy {
            inner: ChainMappingStrategy::Auto,
        }
    }

    /// Create an explicit chain mapping strategy.
    ///
    /// Args:
    ///     mapping: List of (model_chain, native_chain) tuples
    ///
    /// Returns:
    ///     ChainMappingStrategy: Explicit mapping strategy
    #[staticmethod]
    fn explicit(mapping: Vec<(String, String)>) -> Self {
        PyChainMappingStrategy {
            inner: ChainMappingStrategy::Explicit(mapping),
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            ChainMappingStrategy::Auto => "ChainMappingStrategy.auto()".to_string(),
            ChainMappingStrategy::Explicit(m) => {
                format!("ChainMappingStrategy.explicit({:?})", m)
            }
        }
    }
}

/// Options for DockQ calculation.
///
/// Configure contact thresholds and chain mapping strategy.
///
/// Args:
///     contact_threshold: Distance threshold for native contact detection (default: 5.0 A)
///     interface_threshold: Distance threshold for interface residue identification (default: 10.0 A)
///     chain_mapping: Chain mapping strategy (default: automatic)
///
/// Example:
///     >>> options = DockQOptions()  # Default options
///     >>> options = DockQOptions(contact_threshold=4.0)
///     >>> options = DockQOptions(chain_mapping=ChainMappingStrategy.explicit([("A", "A"), ("B", "B")]))
#[pyclass(name = "DockQOptions")]
#[derive(Clone)]
pub struct PyDockQOptions {
    pub(crate) inner: DockQOptions,
}

#[pymethods]
impl PyDockQOptions {
    /// Create DockQ options.
    ///
    /// Args:
    ///     contact_threshold: Distance threshold for fnat contact detection (default: 5.0 A).
    ///                        Atom pairs within this distance are considered contacts.
    ///     interface_threshold: Distance threshold for identifying interface residues
    ///                          for iRMSD (default: 10.0 A).
    ///     chain_mapping: Strategy for chain correspondence (default: automatic).
    #[new]
    #[pyo3(signature = (contact_threshold=5.0, interface_threshold=10.0, chain_mapping=None))]
    fn new(
        contact_threshold: f64,
        interface_threshold: f64,
        chain_mapping: Option<PyChainMappingStrategy>,
    ) -> Self {
        let chain_mapping_strategy = chain_mapping
            .map(|cm| cm.inner)
            .unwrap_or(ChainMappingStrategy::Auto);
        PyDockQOptions {
            inner: DockQOptions {
                contact_threshold,
                interface_threshold,
                chain_mapping: chain_mapping_strategy,
                ..Default::default()
            },
        }
    }

    /// Distance threshold for fnat contact detection (Angstroms).
    #[getter]
    fn contact_threshold(&self) -> f64 {
        self.inner.contact_threshold
    }

    /// Distance threshold for interface residue identification (Angstroms).
    #[getter]
    fn interface_threshold(&self) -> f64 {
        self.inner.interface_threshold
    }

    fn __repr__(&self) -> String {
        let mapping = match &self.inner.chain_mapping {
            ChainMappingStrategy::Auto => "Auto".to_string(),
            ChainMappingStrategy::Explicit(m) => format!("Explicit({:?})", m),
        };
        format!(
            "DockQOptions(contact_threshold={:.1}, interface_threshold={:.1}, chain_mapping={})",
            self.inner.contact_threshold, self.inner.interface_threshold, mapping
        )
    }
}
