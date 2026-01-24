//! Backbone dihedral angle analysis and Ramachandran classification.
//!
//! This module provides high-level APIs for computing backbone dihedral angles
//! (φ, ψ, ω) and performing Ramachandran analysis for protein structure validation.
//!
//! # Overview
//!
//! Backbone dihedral angles are fundamental descriptors of protein conformation:
//! - **φ (phi)**: C(i-1)-N(i)-Cα(i)-C(i), typically -180° to +180°
//! - **ψ (psi)**: N(i)-Cα(i)-C(i)-N(i+1), typically -180° to +180°
//! - **ω (omega)**: Cα(i-1)-C(i-1)-N(i)-Cα(i), usually ~180° (trans) or ~0° (cis)
//!
//! # Ramachandran Regions
//!
//! The Ramachandran plot divides φ/ψ space into favored, allowed, and outlier regions:
//! - **Core/Favored**: Most common conformations (~98% in high-quality structures)
//! - **Allowed**: Less common but sterically acceptable (~2%)
//! - **Outlier**: Sterically strained, often indicates errors (<0.2%)
//!
//! Special residue types have different allowed regions:
//! - **Glycine**: More flexibility due to no side chain
//! - **Proline**: Restricted φ due to cyclic side chain
//! - **Pre-proline**: Residues before proline have modified allowed regions
//!
//! # Example
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Get all backbone dihedrals
//! let dihedrals = structure.phi_psi_angles();
//! for d in &dihedrals {
//!     if let (Some(phi), Some(psi)) = (d.phi, d.psi) {
//!         println!("{}{}: φ={:.1}°, ψ={:.1}°, region={:?}",
//!             d.chain_id, d.residue_seq, phi, psi, d.ramachandran_region);
//!     }
//! }
//!
//! // Find Ramachandran outliers
//! let outliers = structure.ramachandran_outliers();
//! println!("Found {} Ramachandran outliers", outliers.len());
//!
//! // Detect cis peptide bonds
//! let cis_bonds = structure.cis_peptide_bonds();
//! for (res1, res2) in &cis_bonds {
//!     println!("Cis peptide: {}{}-{}{}",
//!         res1.chain_id, res1.residue_seq,
//!         res2.chain_id, res2.residue_seq);
//! }
//!
//! // Get Ramachandran statistics
//! let stats = structure.ramachandran_statistics();
//! println!("Favored: {:.1}%, Allowed: {:.1}%, Outliers: {:.1}%",
//!     stats.favored_fraction * 100.0,
//!     stats.allowed_fraction * 100.0,
//!     stats.outlier_fraction * 100.0);
//! ```

use crate::PdbStructure;

#[cfg(feature = "dssp")]
use crate::dssp::{calculate_all_dihedrals, extract_backbone_atoms};

/// Ramachandran region classification.
///
/// Classifies residues based on their position in the Ramachandran plot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RamachandranRegion {
    /// Core/favored region - most common conformations
    Core,
    /// Allowed region - less common but acceptable
    Allowed,
    /// Generously allowed region - marginal but possible
    Generous,
    /// Outlier region - sterically strained
    Outlier,
    /// Glycine-specific region (more flexibility)
    Glycine,
    /// Proline-specific region (restricted φ)
    Proline,
    /// Pre-proline residue (modified allowed regions)
    PrePro,
    /// Unknown - could not determine (missing angles)
    Unknown,
}

impl RamachandranRegion {
    /// Returns true if this is a favorable region (Core, Allowed, or residue-specific).
    pub fn is_favorable(&self) -> bool {
        matches!(
            self,
            RamachandranRegion::Core
                | RamachandranRegion::Allowed
                | RamachandranRegion::Glycine
                | RamachandranRegion::Proline
                | RamachandranRegion::PrePro
        )
    }

    /// Returns true if this is an outlier.
    pub fn is_outlier(&self) -> bool {
        matches!(self, RamachandranRegion::Outlier)
    }
}

/// Per-residue dihedral angles with Ramachandran classification.
#[derive(Debug, Clone)]
pub struct ResidueDihedrals {
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Residue name (3-letter code)
    pub residue_name: String,
    /// Insertion code (if any)
    pub ins_code: Option<char>,
    /// φ (phi) angle in degrees, or None if not calculable
    pub phi: Option<f64>,
    /// ψ (psi) angle in degrees, or None if not calculable
    pub psi: Option<f64>,
    /// ω (omega) angle in degrees, or None if not calculable
    pub omega: Option<f64>,
    /// Ramachandran region classification
    pub ramachandran_region: RamachandranRegion,
}

impl ResidueDihedrals {
    /// Returns true if this residue has both φ and ψ angles.
    pub fn has_phi_psi(&self) -> bool {
        self.phi.is_some() && self.psi.is_some()
    }

    /// Returns true if this is a cis peptide bond (ω ≈ 0°).
    ///
    /// A peptide bond is considered cis if |ω| < 30°.
    pub fn is_cis_peptide(&self) -> bool {
        if let Some(omega) = self.omega {
            omega.abs() < 30.0
        } else {
            false
        }
    }

    /// Returns true if this is a trans peptide bond (ω ≈ 180°).
    pub fn is_trans_peptide(&self) -> bool {
        if let Some(omega) = self.omega {
            (omega.abs() - 180.0).abs() < 30.0
        } else {
            false
        }
    }
}

/// Reference to a residue for cis peptide bond reporting.
#[derive(Debug, Clone)]
pub struct ResidueRef {
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Residue name (3-letter code)
    pub residue_name: String,
    /// Insertion code (if any)
    pub ins_code: Option<char>,
}

/// Ramachandran plot statistics for structure validation.
#[derive(Debug, Clone)]
pub struct RamachandranStats {
    /// Total number of residues analyzed
    pub total_residues: usize,
    /// Number of residues in favored/core region
    pub favored_count: usize,
    /// Number of residues in allowed region
    pub allowed_count: usize,
    /// Number of residues in outlier region
    pub outlier_count: usize,
    /// Number of glycine residues
    pub glycine_count: usize,
    /// Number of proline residues
    pub proline_count: usize,
    /// Number of pre-proline residues
    pub prepro_count: usize,
    /// Fraction of residues in favored region
    pub favored_fraction: f64,
    /// Fraction of residues in allowed region
    pub allowed_fraction: f64,
    /// Fraction of residues in outlier region
    pub outlier_fraction: f64,
    /// Number of cis peptide bonds
    pub cis_peptide_count: usize,
    /// Number of cis non-proline peptide bonds (unusual)
    pub cis_nonpro_count: usize,
}

impl Default for RamachandranStats {
    fn default() -> Self {
        Self {
            total_residues: 0,
            favored_count: 0,
            allowed_count: 0,
            outlier_count: 0,
            glycine_count: 0,
            proline_count: 0,
            prepro_count: 0,
            favored_fraction: 0.0,
            allowed_fraction: 0.0,
            outlier_fraction: 0.0,
            cis_peptide_count: 0,
            cis_nonpro_count: 0,
        }
    }
}

// Ramachandran region boundaries (in degrees)
// These are simplified boundaries based on common validation criteria

/// Check if φ/ψ is in the α-helix favored region
fn is_alpha_helix_core(phi: f64, psi: f64) -> bool {
    // α-helix: φ ≈ -60°, ψ ≈ -45°
    let phi_diff = (phi - (-62.0)).abs();
    let psi_diff = (psi - (-41.0)).abs();
    phi_diff < 25.0 && psi_diff < 25.0
}

/// Check if φ/ψ is in the β-sheet favored region
fn is_beta_sheet_core(phi: f64, psi: f64) -> bool {
    // β-sheet: φ ≈ -120°, ψ ≈ +130°
    let phi_diff = (phi - (-119.0)).abs();
    let psi_in_range = psi > 100.0 && psi < 180.0;
    phi_diff < 35.0 && psi_in_range
}

/// Check if φ/ψ is in the left-handed helix region (for glycine)
fn is_left_helix(phi: f64, psi: f64) -> bool {
    // Left-handed helix: φ ≈ +60°, ψ ≈ +45°
    let phi_diff = (phi - 60.0).abs();
    let psi_diff = (psi - 45.0).abs();
    phi_diff < 30.0 && psi_diff < 30.0
}

/// Check if φ/ψ is in the PPII (polyproline II) region
fn is_ppii_region(phi: f64, psi: f64) -> bool {
    // PPII: φ ≈ -75°, ψ ≈ +145°
    let phi_diff = (phi - (-75.0)).abs();
    let psi_diff = (psi - 145.0).abs();
    phi_diff < 29.0 && psi_diff < 29.0
}

/// Check if φ is in the proline-allowed range
fn is_proline_phi(phi: f64) -> bool {
    // Proline: φ restricted to around -60° due to cyclic side chain
    phi > -80.0 && phi < -40.0
}

/// Classify a residue's Ramachandran region
fn classify_ramachandran(
    phi: f64,
    psi: f64,
    residue_name: &str,
    next_residue_name: Option<&str>,
) -> RamachandranRegion {
    // Handle special residue types
    if residue_name == "GLY" {
        // Glycine has more flexibility
        if is_alpha_helix_core(phi, psi) || is_beta_sheet_core(phi, psi) || is_left_helix(phi, psi)
        {
            return RamachandranRegion::Glycine;
        }
        // Glycine has extended allowed regions
        if phi > -180.0 && phi < 180.0 && psi > -180.0 && psi < 180.0 {
            // Very permissive for glycine, but flag obvious outliers
            if (phi.abs() < 20.0 && psi.abs() > 150.0) || (phi > 100.0 && psi < -100.0) {
                return RamachandranRegion::Outlier;
            }
            return RamachandranRegion::Glycine;
        }
    }

    if residue_name == "PRO" {
        // Proline has restricted φ
        if is_proline_phi(phi) {
            if is_alpha_helix_core(phi, psi) || is_ppii_region(phi, psi) {
                return RamachandranRegion::Proline;
            }
            // Extended allowed for proline
            if psi > -60.0 && psi < 180.0 {
                return RamachandranRegion::Proline;
            }
        }
        // Proline with unusual φ is an outlier
        return RamachandranRegion::Outlier;
    }

    // Check if this is a pre-proline residue
    if let Some(next) = next_residue_name {
        if next == "PRO" {
            // Pre-proline residues have shifted allowed regions
            if is_beta_sheet_core(phi, psi) || is_ppii_region(phi, psi) {
                return RamachandranRegion::PrePro;
            }
            // Slightly different alpha region for pre-pro
            let phi_diff = (phi - (-60.0)).abs();
            let psi_diff = (psi - (-30.0)).abs();
            if phi_diff < 30.0 && psi_diff < 30.0 {
                return RamachandranRegion::PrePro;
            }
        }
    }

    // General residues
    if is_alpha_helix_core(phi, psi) || is_beta_sheet_core(phi, psi) {
        return RamachandranRegion::Core;
    }

    // Check allowed regions (broader boundaries)
    // Extended alpha region
    let phi_diff_alpha = (phi - (-62.0)).abs();
    let psi_diff_alpha = (psi - (-41.0)).abs();
    if phi_diff_alpha < 40.0 && psi_diff_alpha < 40.0 {
        return RamachandranRegion::Allowed;
    }

    // Extended beta region
    if phi < -60.0 && phi > -180.0 && psi > 80.0 && psi < 180.0 {
        return RamachandranRegion::Allowed;
    }

    // PPII region
    if is_ppii_region(phi, psi) {
        return RamachandranRegion::Allowed;
    }

    // Check generously allowed (even broader)
    // Generously allowed alpha
    if phi > -150.0 && phi < 0.0 && psi > -100.0 && psi < 50.0 {
        return RamachandranRegion::Generous;
    }

    // Generously allowed beta
    if phi < -60.0 && psi > 60.0 {
        return RamachandranRegion::Generous;
    }

    // Everything else is an outlier
    RamachandranRegion::Outlier
}

#[cfg(feature = "dssp")]
impl PdbStructure {
    /// Computes backbone dihedral angles (φ, ψ, ω) for all residues.
    ///
    /// Returns a vector of `ResidueDihedrals` containing the backbone angles
    /// and Ramachandran classification for each residue.
    ///
    /// Terminal residues will have `None` for angles that cannot be calculated
    /// (φ for N-terminus, ψ for C-terminus).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let dihedrals = structure.phi_psi_angles();
    /// for d in &dihedrals {
    ///     if let (Some(phi), Some(psi)) = (d.phi, d.psi) {
    ///         println!("{}: φ={:.1}°, ψ={:.1}°", d.residue_name, phi, psi);
    ///     }
    /// }
    /// ```
    pub fn phi_psi_angles(&self) -> Vec<ResidueDihedrals> {
        let backbone_atoms = extract_backbone_atoms(&self.atoms);
        if backbone_atoms.is_empty() {
            return Vec::new();
        }

        let dihedrals = calculate_all_dihedrals(&backbone_atoms);

        backbone_atoms
            .iter()
            .zip(dihedrals.iter())
            .enumerate()
            .map(|(i, (backbone, dihedral))| {
                // Get next residue name for pre-proline detection
                let next_res_name = if i + 1 < backbone_atoms.len()
                    && backbone_atoms[i + 1].chain_id == backbone.chain_id
                {
                    Some(backbone_atoms[i + 1].residue_name.as_str())
                } else {
                    None
                };

                // Convert phi and psi to standard biochemistry convention (IUPAC)
                // The internal DSSP calculation uses an opposite sign convention
                let phi_std = dihedral.phi.map(|p| -p);
                let psi_std = dihedral.psi.map(|p| -p);

                let region = match (phi_std, psi_std) {
                    (Some(phi), Some(psi)) => {
                        classify_ramachandran(phi, psi, &backbone.residue_name, next_res_name)
                    }
                    _ => RamachandranRegion::Unknown,
                };

                ResidueDihedrals {
                    chain_id: backbone.chain_id.clone(),
                    residue_seq: backbone.residue_seq,
                    residue_name: backbone.residue_name.clone(),
                    ins_code: backbone.ins_code,
                    phi: phi_std,
                    psi: psi_std,
                    omega: dihedral.omega,
                    ramachandran_region: region,
                }
            })
            .collect()
    }

    /// Returns residues with Ramachandran outliers.
    ///
    /// These are residues with φ/ψ angles in sterically strained regions,
    /// which may indicate modeling errors or genuine strain.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let outliers = structure.ramachandran_outliers();
    /// for res in &outliers {
    ///     println!("Outlier: {}{} {} (φ={:.1}°, ψ={:.1}°)",
    ///         res.chain_id, res.residue_seq, res.residue_name,
    ///         res.phi.unwrap(), res.psi.unwrap());
    /// }
    /// ```
    pub fn ramachandran_outliers(&self) -> Vec<ResidueDihedrals> {
        self.phi_psi_angles()
            .into_iter()
            .filter(|d| d.ramachandran_region == RamachandranRegion::Outlier)
            .collect()
    }

    /// Detects cis peptide bonds in the structure.
    ///
    /// Cis peptide bonds have ω ≈ 0° instead of the usual ω ≈ 180° (trans).
    /// Cis-proline is relatively common (~5% of prolines), but cis non-proline
    /// is rare and may indicate an error.
    ///
    /// Returns pairs of (residue_i-1, residue_i) where the peptide bond
    /// between them is cis.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let cis_bonds = structure.cis_peptide_bonds();
    /// for (res1, res2) in &cis_bonds {
    ///     println!("Cis peptide: {}-{}", res1.residue_name, res2.residue_name);
    ///     if res2.residue_name != "PRO" {
    ///         println!("  Warning: cis non-proline!");
    ///     }
    /// }
    /// ```
    pub fn cis_peptide_bonds(&self) -> Vec<(ResidueRef, ResidueRef)> {
        let dihedrals = self.phi_psi_angles();
        let mut cis_bonds = Vec::new();

        for i in 1..dihedrals.len() {
            if dihedrals[i].is_cis_peptide() && dihedrals[i].chain_id == dihedrals[i - 1].chain_id {
                let res1 = ResidueRef {
                    chain_id: dihedrals[i - 1].chain_id.clone(),
                    residue_seq: dihedrals[i - 1].residue_seq,
                    residue_name: dihedrals[i - 1].residue_name.clone(),
                    ins_code: dihedrals[i - 1].ins_code,
                };
                let res2 = ResidueRef {
                    chain_id: dihedrals[i].chain_id.clone(),
                    residue_seq: dihedrals[i].residue_seq,
                    residue_name: dihedrals[i].residue_name.clone(),
                    ins_code: dihedrals[i].ins_code,
                };
                cis_bonds.push((res1, res2));
            }
        }

        cis_bonds
    }

    /// Computes Ramachandran plot statistics for structure validation.
    ///
    /// Returns comprehensive statistics about the distribution of residues
    /// across Ramachandran regions.
    ///
    /// # Quality Criteria
    ///
    /// For well-refined structures:
    /// - Favored: >98%
    /// - Allowed: <2%
    /// - Outliers: <0.2%
    ///
    /// # Example
    ///
    /// ```ignore
    /// let stats = structure.ramachandran_statistics();
    /// println!("Ramachandran Statistics:");
    /// println!("  Favored: {:.1}% ({}/{})",
    ///     stats.favored_fraction * 100.0, stats.favored_count, stats.total_residues);
    /// println!("  Allowed: {:.1}%", stats.allowed_fraction * 100.0);
    /// println!("  Outliers: {:.1}% ({})",
    ///     stats.outlier_fraction * 100.0, stats.outlier_count);
    /// println!("  Cis peptides: {} ({} non-Pro)",
    ///     stats.cis_peptide_count, stats.cis_nonpro_count);
    /// ```
    pub fn ramachandran_statistics(&self) -> RamachandranStats {
        let dihedrals = self.phi_psi_angles();
        let mut stats = RamachandranStats::default();

        // Count residues with both φ and ψ (terminal residues excluded)
        let valid_residues: Vec<_> = dihedrals.iter().filter(|d| d.has_phi_psi()).collect();
        stats.total_residues = valid_residues.len();

        if stats.total_residues == 0 {
            return stats;
        }

        for d in &valid_residues {
            match d.ramachandran_region {
                RamachandranRegion::Core => stats.favored_count += 1,
                RamachandranRegion::Allowed | RamachandranRegion::Generous => {
                    stats.allowed_count += 1
                }
                RamachandranRegion::Outlier => stats.outlier_count += 1,
                RamachandranRegion::Glycine => {
                    stats.glycine_count += 1;
                    stats.favored_count += 1; // Glycine in allowed region counts as favored
                }
                RamachandranRegion::Proline => {
                    stats.proline_count += 1;
                    stats.favored_count += 1;
                }
                RamachandranRegion::PrePro => {
                    stats.prepro_count += 1;
                    stats.favored_count += 1;
                }
                RamachandranRegion::Unknown => {}
            }
        }

        // Count cis peptides
        for i in 1..dihedrals.len() {
            if dihedrals[i].is_cis_peptide() && dihedrals[i].chain_id == dihedrals[i - 1].chain_id {
                stats.cis_peptide_count += 1;
                if dihedrals[i].residue_name != "PRO" {
                    stats.cis_nonpro_count += 1;
                }
            }
        }

        // Calculate fractions
        let total = stats.total_residues as f64;
        stats.favored_fraction = stats.favored_count as f64 / total;
        stats.allowed_fraction = stats.allowed_count as f64 / total;
        stats.outlier_fraction = stats.outlier_count as f64 / total;

        stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ramachandran_region_methods() {
        assert!(RamachandranRegion::Core.is_favorable());
        assert!(RamachandranRegion::Allowed.is_favorable());
        assert!(RamachandranRegion::Glycine.is_favorable());
        assert!(RamachandranRegion::Proline.is_favorable());
        assert!(!RamachandranRegion::Outlier.is_favorable());
        assert!(!RamachandranRegion::Generous.is_favorable());

        assert!(RamachandranRegion::Outlier.is_outlier());
        assert!(!RamachandranRegion::Core.is_outlier());
    }

    #[test]
    fn test_classify_alpha_helix() {
        // Typical α-helix: φ ≈ -60°, ψ ≈ -45°
        let region = classify_ramachandran(-62.0, -41.0, "ALA", None);
        assert_eq!(region, RamachandranRegion::Core);
    }

    #[test]
    fn test_classify_beta_sheet() {
        // Typical β-sheet: φ ≈ -120°, ψ ≈ +130°
        let region = classify_ramachandran(-120.0, 130.0, "ALA", None);
        assert_eq!(region, RamachandranRegion::Core);
    }

    #[test]
    fn test_classify_glycine() {
        // Glycine in left-handed helix region
        let region = classify_ramachandran(60.0, 45.0, "GLY", None);
        assert_eq!(region, RamachandranRegion::Glycine);

        // Glycine in alpha region
        let region = classify_ramachandran(-62.0, -41.0, "GLY", None);
        assert_eq!(region, RamachandranRegion::Glycine);
    }

    #[test]
    fn test_classify_proline() {
        // Proline in allowed region
        let region = classify_ramachandran(-60.0, 145.0, "PRO", None);
        assert_eq!(region, RamachandranRegion::Proline);
    }

    #[test]
    fn test_classify_outlier() {
        // Clear outlier for general residue
        let region = classify_ramachandran(60.0, 60.0, "ALA", None);
        assert_eq!(region, RamachandranRegion::Outlier);
    }

    #[test]
    fn test_residue_dihedrals_cis_trans() {
        let mut d = ResidueDihedrals {
            chain_id: "A".to_string(),
            residue_seq: 1,
            residue_name: "ALA".to_string(),
            ins_code: None,
            phi: Some(-60.0),
            psi: Some(-45.0),
            omega: Some(180.0),
            ramachandran_region: RamachandranRegion::Core,
        };

        assert!(d.is_trans_peptide());
        assert!(!d.is_cis_peptide());

        d.omega = Some(0.0);
        assert!(d.is_cis_peptide());
        assert!(!d.is_trans_peptide());
    }

    #[test]
    fn test_ramachandran_stats_default() {
        let stats = RamachandranStats::default();
        assert_eq!(stats.total_residues, 0);
        assert_eq!(stats.favored_count, 0);
        assert_eq!(stats.outlier_count, 0);
    }
}
