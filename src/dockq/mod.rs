//! DockQ v2 interface quality assessment for protein-protein complexes.
//!
//! This module implements the DockQ scoring method for evaluating the quality
//! of modeled protein-protein interfaces. DockQ is the standard metric used
//! in CAPRI and CASP-multimer evaluations.
//!
//! # Feature Flag
//!
//! Enable this module in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.6", features = ["dockq"] }
//! ```
//!
//! # Quick Start
//!
//! ```rust,ignore
//! use pdbrust::{parse_pdb_file, PdbStructure};
//!
//! let model = parse_pdb_file("model.pdb")?;
//! let native = parse_pdb_file("native.pdb")?;
//!
//! // Compute DockQ with automatic chain mapping
//! let result = model.dockq_to(&native)?;
//! println!("DockQ: {:.4}", result.total_dockq);
//!
//! for iface in &result.interfaces {
//!     println!("Interface {}-{}: DockQ={:.3} fnat={:.3} iRMSD={:.2} LRMSD={:.2}",
//!         iface.native_chains.0, iface.native_chains.1,
//!         iface.dockq, iface.fnat, iface.irmsd, iface.lrmsd);
//! }
//! ```
//!
//! # DockQ Score Components
//!
//! DockQ combines three metrics:
//!
//! - **fnat**: Fraction of native contacts preserved in the model
//! - **iRMSD**: Interface RMSD (backbone atoms at the interface)
//! - **LRMSD**: Ligand RMSD (backbone RMSD of the smaller chain after
//!   receptor alignment)
//!
//! The combined score: `DockQ = (fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD/8.5)^2)) / 3`
//!
//! # Quality Classification
//!
//! | Quality | DockQ Range |
//! |---------|-------------|
//! | Incorrect | < 0.23 |
//! | Acceptable | 0.23 - 0.49 |
//! | Medium | 0.49 - 0.80 |
//! | High | >= 0.80 |

mod chain_mapping;
mod contacts;
mod metrics;
pub mod sequence_align;

use std::collections::HashMap;

use crate::core::PdbStructure;
use crate::error::PdbError;
use crate::geometry::AtomSelection;

pub use chain_mapping::find_chain_mapping;
pub use contacts::ResidueKey;
pub use sequence_align::{AlignmentParams, SequenceAlignment, align_sequences, sequence_identity};

/// Options for DockQ calculation.
#[derive(Debug, Clone)]
pub struct DockQOptions {
    /// Distance threshold for fnat contact detection (default: 5.0 A).
    pub contact_threshold: f64,
    /// Distance threshold for identifying interface residues (default: 10.0 A).
    pub interface_threshold: f64,
    /// Atom selection for fnat contact detection (default: AllAtoms for heavy atoms).
    pub fnat_atom_selection: AtomSelection,
    /// Atom selection for iRMSD/LRMSD (default: Backbone).
    pub rmsd_atom_selection: AtomSelection,
    /// Strategy for establishing chain correspondence.
    pub chain_mapping: ChainMappingStrategy,
}

impl Default for DockQOptions {
    fn default() -> Self {
        Self {
            contact_threshold: 5.0,
            interface_threshold: 10.0,
            fnat_atom_selection: AtomSelection::AllAtoms,
            rmsd_atom_selection: AtomSelection::Backbone,
            chain_mapping: ChainMappingStrategy::Auto,
        }
    }
}

/// Strategy for establishing chain correspondence between model and native.
#[derive(Debug, Clone)]
pub enum ChainMappingStrategy {
    /// Automatic: use sequence alignment to determine chain mapping.
    Auto,
    /// Explicit mapping: vec of (model_chain, native_chain) pairs.
    Explicit(Vec<(String, String)>),
}

/// Result for the entire complex comparison.
#[derive(Debug, Clone)]
pub struct DockQResult {
    /// Per-interface results.
    pub interfaces: Vec<InterfaceResult>,
    /// Average DockQ over all native interfaces (equal weighting).
    pub total_dockq: f64,
    /// Final model-to-native chain mapping used.
    pub chain_mapping: Vec<(String, String)>,
    /// Number of interfaces evaluated.
    pub num_interfaces: usize,
}

/// Result for a single chain-pair interface.
#[derive(Debug, Clone)]
pub struct InterfaceResult {
    /// Chain pair in the native structure (receptor, ligand).
    pub native_chains: (String, String),
    /// Corresponding chain pair in the model.
    pub model_chains: (String, String),
    /// Fraction of native contacts preserved [0, 1].
    pub fnat: f64,
    /// Fraction of non-native contacts in the model.
    pub fnonnat: f64,
    /// Interface RMSD in Angstroms.
    pub irmsd: f64,
    /// Ligand RMSD in Angstroms.
    pub lrmsd: f64,
    /// DockQ score [0, 1].
    pub dockq: f64,
    /// Quality classification.
    pub quality: DockQQuality,
    /// F1 score for contacts.
    pub f1: f64,
    /// Number of contacts in native.
    pub num_native_contacts: usize,
    /// Number of contacts in model.
    pub num_model_contacts: usize,
    /// Number of atomic clashes (< 2.0 A).
    pub num_clashes: usize,
}

/// DockQ quality classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DockQQuality {
    /// DockQ < 0.23
    Incorrect,
    /// 0.23 <= DockQ < 0.49
    Acceptable,
    /// 0.49 <= DockQ < 0.80
    Medium,
    /// DockQ >= 0.80
    High,
}

impl DockQQuality {
    /// Classify a DockQ score.
    pub fn from_score(dockq: f64) -> Self {
        if dockq >= 0.80 {
            DockQQuality::High
        } else if dockq >= 0.49 {
            DockQQuality::Medium
        } else if dockq >= 0.23 {
            DockQQuality::Acceptable
        } else {
            DockQQuality::Incorrect
        }
    }
}

impl std::fmt::Display for DockQQuality {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DockQQuality::Incorrect => write!(f, "Incorrect"),
            DockQQuality::Acceptable => write!(f, "Acceptable"),
            DockQQuality::Medium => write!(f, "Medium"),
            DockQQuality::High => write!(f, "High"),
        }
    }
}

/// Compute DockQ for a specific interface (two chain pairs).
pub fn calculate_interface_dockq(
    model: &PdbStructure,
    native: &PdbStructure,
    model_chains: (&str, &str),
    native_chains: (&str, &str),
    options: &DockQOptions,
) -> Result<InterfaceResult, PdbError> {
    // Build residue mapping for both chain pairs
    let mut residue_mapping: HashMap<ResidueKey, ResidueKey> = HashMap::new();
    let map_a =
        chain_mapping::build_residue_mapping(model, native, model_chains.0, native_chains.0);
    let map_b =
        chain_mapping::build_residue_mapping(model, native, model_chains.1, native_chains.1);
    residue_mapping.extend(map_a);
    residue_mapping.extend(map_b);

    // Compute fnat
    let native_contacts = contacts::find_interface_contacts(
        native,
        native_chains.0,
        native_chains.1,
        options.contact_threshold,
    );
    let model_contacts = contacts::find_interface_contacts(
        model,
        model_chains.0,
        model_chains.1,
        options.contact_threshold,
    );

    if native_contacts.is_empty() {
        return Err(PdbError::NoInterfaceContacts(format!(
            "No contacts between chains {} and {} in native (threshold: {} A)",
            native_chains.0, native_chains.1, options.contact_threshold
        )));
    }

    let contact_result =
        contacts::compare_contacts(&native_contacts, &model_contacts, &residue_mapping);

    // Count clashes
    let num_clashes = contacts::count_clashes(model, model_chains.0, model_chains.1);

    // Compute iRMSD
    let irmsd = metrics::compute_irmsd(
        model,
        native,
        native_chains.0,
        native_chains.1,
        model_chains.0,
        model_chains.1,
        &residue_mapping,
        options.interface_threshold,
    )
    .unwrap_or(f64::INFINITY);

    // Compute LRMSD
    let lrmsd = metrics::compute_lrmsd(
        model,
        native,
        native_chains.0,
        native_chains.1,
        model_chains.0,
        model_chains.1,
    )
    .unwrap_or(f64::INFINITY);

    let dockq = metrics::compute_dockq(contact_result.fnat, irmsd, lrmsd);
    let quality = DockQQuality::from_score(dockq);

    Ok(InterfaceResult {
        native_chains: (native_chains.0.to_string(), native_chains.1.to_string()),
        model_chains: (model_chains.0.to_string(), model_chains.1.to_string()),
        fnat: contact_result.fnat,
        fnonnat: contact_result.fnonnat,
        irmsd,
        lrmsd,
        dockq,
        quality,
        f1: contact_result.f1,
        num_native_contacts: contact_result.num_native_contacts,
        num_model_contacts: contact_result.num_model_contacts,
        num_clashes,
    })
}

/// Compute DockQ for all interfaces in a complex.
fn compute_all_interfaces(
    model: &PdbStructure,
    native: &PdbStructure,
    mapping: &[(String, String)],
    options: &DockQOptions,
) -> Result<DockQResult, PdbError> {
    // Find all native chain pairs that have interface contacts
    let native_chain_ids: Vec<&str> = mapping.iter().map(|(_, n)| n.as_str()).collect();

    let mut interfaces = Vec::new();

    for i in 0..native_chain_ids.len() {
        for j in (i + 1)..native_chain_ids.len() {
            let nc_a = native_chain_ids[i];
            let nc_b = native_chain_ids[j];

            // Check if there are contacts at this interface
            let native_contacts =
                contacts::find_interface_contacts(native, nc_a, nc_b, options.contact_threshold);

            if native_contacts.is_empty() {
                continue;
            }

            // Find corresponding model chains
            let mc_a = mapping
                .iter()
                .find(|(_, n)| n == nc_a)
                .map(|(m, _)| m.as_str());
            let mc_b = mapping
                .iter()
                .find(|(_, n)| n == nc_b)
                .map(|(m, _)| m.as_str());

            if let (Some(mc_a), Some(mc_b)) = (mc_a, mc_b) {
                match calculate_interface_dockq(model, native, (mc_a, mc_b), (nc_a, nc_b), options)
                {
                    Ok(iface_result) => interfaces.push(iface_result),
                    Err(_) => continue, // Skip interfaces that can't be computed
                }
            }
        }
    }

    if interfaces.is_empty() {
        return Err(PdbError::NoInterfaceContacts(
            "No interfaces with contacts found between mapped chains".to_string(),
        ));
    }

    // Average DockQ over all native interfaces (equal weighting, matching DockQ v2)
    let total_dockq = interfaces.iter().map(|i| i.dockq).sum::<f64>() / interfaces.len() as f64;

    let num_interfaces = interfaces.len();

    Ok(DockQResult {
        interfaces,
        total_dockq,
        chain_mapping: mapping.to_vec(),
        num_interfaces,
    })
}

/// Extension methods for PdbStructure providing DockQ functionality.
impl PdbStructure {
    /// Compute DockQ between this (model) and a reference (native) complex.
    ///
    /// Automatically detects chain mapping using sequence alignment and
    /// evaluates all interfaces with contacts.
    ///
    /// # Arguments
    /// * `native` - The reference (native/crystal) structure
    ///
    /// # Returns
    /// [`DockQResult`] containing per-interface scores and overall DockQ
    ///
    /// # Errors
    /// - `NoChainMapping` if no chain correspondence can be found
    /// - `NoInterfaceContacts` if no interfaces with contacts are found
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let result = model.dockq_to(&native)?;
    /// println!("DockQ: {:.4} ({})", result.total_dockq,
    ///     result.interfaces[0].quality);
    /// ```
    pub fn dockq_to(&self, native: &PdbStructure) -> Result<DockQResult, PdbError> {
        self.dockq_to_with_options(native, DockQOptions::default())
    }

    /// Compute DockQ with custom options.
    ///
    /// # Arguments
    /// * `native` - The reference (native/crystal) structure
    /// * `options` - Custom options for thresholds, atom selection, and chain mapping
    ///
    /// # Returns
    /// [`DockQResult`] containing per-interface scores and overall DockQ
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use pdbrust::dockq::{DockQOptions, ChainMappingStrategy};
    ///
    /// let options = DockQOptions {
    ///     chain_mapping: ChainMappingStrategy::Explicit(vec![
    ///         ("A".to_string(), "A".to_string()),
    ///         ("B".to_string(), "B".to_string()),
    ///     ]),
    ///     ..Default::default()
    /// };
    ///
    /// let result = model.dockq_to_with_options(&native, options)?;
    /// ```
    pub fn dockq_to_with_options(
        &self,
        native: &PdbStructure,
        options: DockQOptions,
    ) -> Result<DockQResult, PdbError> {
        let mapping = match &options.chain_mapping {
            ChainMappingStrategy::Auto => find_chain_mapping(self, native)?,
            ChainMappingStrategy::Explicit(m) => m.clone(),
        };

        compute_all_interfaces(self, native, &mapping, &options)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    #[allow(clippy::too_many_arguments)]
    fn create_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
        element: &str,
    ) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: element.to_string(),
        }
    }

    fn create_two_chain_complex() -> PdbStructure {
        let mut structure = PdbStructure::new();
        // Chain A: 4 residues along x-axis
        // Chain B: 4 residues along x-axis, offset in y by 4.0 A (within contact distance)
        let mut serial = 1;
        for i in 0..4 {
            let x = i as f64 * 3.8;
            let resname = match i {
                0 => "ALA",
                1 => "GLY",
                2 => "VAL",
                _ => "LEU",
            };
            // Chain A backbone
            structure.atoms.push(create_atom(
                serial,
                "N",
                resname,
                "A",
                i + 1,
                x - 0.5,
                0.0,
                0.0,
                "N",
            ));
            serial += 1;
            structure.atoms.push(create_atom(
                serial,
                "CA",
                resname,
                "A",
                i + 1,
                x,
                0.0,
                0.0,
                "C",
            ));
            serial += 1;
            structure.atoms.push(create_atom(
                serial,
                "C",
                resname,
                "A",
                i + 1,
                x + 0.5,
                0.0,
                0.0,
                "C",
            ));
            serial += 1;
            structure.atoms.push(create_atom(
                serial,
                "O",
                resname,
                "A",
                i + 1,
                x + 0.5,
                0.5,
                0.0,
                "O",
            ));
            serial += 1;
        }
        for i in 0..4 {
            let x = i as f64 * 3.8;
            let resname = match i {
                0 => "ILE",
                1 => "PHE",
                2 => "TRP",
                _ => "TYR",
            };
            // Chain B backbone, 4.0 A above chain A in y
            structure.atoms.push(create_atom(
                serial,
                "N",
                resname,
                "B",
                i + 1,
                x - 0.5,
                4.0,
                0.0,
                "N",
            ));
            serial += 1;
            structure.atoms.push(create_atom(
                serial,
                "CA",
                resname,
                "B",
                i + 1,
                x,
                4.0,
                0.0,
                "C",
            ));
            serial += 1;
            structure.atoms.push(create_atom(
                serial,
                "C",
                resname,
                "B",
                i + 1,
                x + 0.5,
                4.0,
                0.0,
                "C",
            ));
            serial += 1;
            structure.atoms.push(create_atom(
                serial,
                "O",
                resname,
                "B",
                i + 1,
                x + 0.5,
                4.5,
                0.0,
                "O",
            ));
            serial += 1;
        }
        structure
    }

    #[test]
    fn test_dockq_self_comparison() {
        let structure = create_two_chain_complex();
        let result = structure.dockq_to(&structure).unwrap();

        assert!(
            (result.total_dockq - 1.0).abs() < 1e-6,
            "Self-comparison DockQ should be 1.0, got {}",
            result.total_dockq
        );

        for iface in &result.interfaces {
            assert!(
                (iface.fnat - 1.0).abs() < 1e-6,
                "Self-comparison fnat should be 1.0"
            );
            assert!(iface.irmsd < 1e-6, "Self-comparison iRMSD should be ~0");
            assert!(iface.lrmsd < 1e-6, "Self-comparison LRMSD should be ~0");
            assert_eq!(iface.quality, DockQQuality::High);
        }
    }

    #[test]
    fn test_dockq_quality_classification() {
        assert_eq!(DockQQuality::from_score(0.0), DockQQuality::Incorrect);
        assert_eq!(DockQQuality::from_score(0.22), DockQQuality::Incorrect);
        assert_eq!(DockQQuality::from_score(0.23), DockQQuality::Acceptable);
        assert_eq!(DockQQuality::from_score(0.48), DockQQuality::Acceptable);
        assert_eq!(DockQQuality::from_score(0.49), DockQQuality::Medium);
        assert_eq!(DockQQuality::from_score(0.79), DockQQuality::Medium);
        assert_eq!(DockQQuality::from_score(0.80), DockQQuality::High);
        assert_eq!(DockQQuality::from_score(1.0), DockQQuality::High);
    }

    #[test]
    fn test_dockq_explicit_mapping() {
        let structure = create_two_chain_complex();
        let options = DockQOptions {
            chain_mapping: ChainMappingStrategy::Explicit(vec![
                ("A".to_string(), "A".to_string()),
                ("B".to_string(), "B".to_string()),
            ]),
            ..Default::default()
        };

        let result = structure
            .dockq_to_with_options(&structure, options)
            .unwrap();
        assert!((result.total_dockq - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dockq_translated_model() {
        let native = create_two_chain_complex();
        let mut model = native.clone();

        // Translate all atoms - contacts should be preserved (fnat=1)
        // but iRMSD and LRMSD should be > 0
        for atom in &mut model.atoms {
            atom.x += 2.0;
            atom.y += 1.0;
        }

        let options = DockQOptions {
            chain_mapping: ChainMappingStrategy::Explicit(vec![
                ("A".to_string(), "A".to_string()),
                ("B".to_string(), "B".to_string()),
            ]),
            ..Default::default()
        };

        let result = model.dockq_to_with_options(&native, options).unwrap();

        // Contacts should be preserved since both chains moved together
        for iface in &result.interfaces {
            assert!(
                (iface.fnat - 1.0).abs() < 1e-6,
                "fnat should be 1.0 for rigid body translation, got {}",
                iface.fnat
            );
        }
    }
}
