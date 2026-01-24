//! Hydrogen bond network analysis for protein structures.
//!
//! This module provides high-level APIs for analyzing mainchain (backbone)
//! hydrogen bonds in protein structures, building on the DSSP algorithm.
//!
//! # Overview
//!
//! Hydrogen bonds are fundamental to protein secondary structure:
//! - **α-helices**: i → i+4 H-bond pattern (N-H...O=C)
//! - **3₁₀-helices**: i → i+3 H-bond pattern
//! - **π-helices**: i → i+5 H-bond pattern
//! - **β-sheets**: H-bonds between strands (parallel or antiparallel)
//! - **Turns**: Short-range H-bonds
//!
//! # H-Bond Energy
//!
//! H-bonds are detected using the Kabsch-Sander electrostatic model:
//!
//! ```text
//! E = 27.888 × (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol
//! ```
//!
//! An H-bond is assigned when E < -0.5 kcal/mol. More negative energies
//! indicate stronger bonds.
//!
//! # Example
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Get all mainchain H-bonds
//! let hbonds = structure.mainchain_hbonds();
//! for hb in &hbonds {
//!     println!("{}{} {} → {}{} {} ({:.2} kcal/mol, {:?})",
//!         hb.donor_chain, hb.donor_resid, hb.donor_resname,
//!         hb.acceptor_chain, hb.acceptor_resid, hb.acceptor_resname,
//!         hb.energy, hb.hbond_type);
//! }
//!
//! // Get H-bonds for a specific residue
//! let res_hbonds = structure.hbonds_for_residue("A", 42);
//! println!("Residue 42 donates: {:?}", res_hbonds.donated.len());
//! println!("Residue 42 accepts: {:?}", res_hbonds.accepted.len());
//!
//! // Get H-bond statistics
//! let stats = structure.hbond_statistics();
//! println!("Total H-bonds: {}", stats.total_hbonds);
//! println!("Intra-helical: {}", stats.intra_helical);
//! println!("Beta-sheet: {}", stats.beta_sheet);
//! ```

use crate::PdbStructure;

#[cfg(feature = "dssp")]
use crate::dssp::{
    BackboneAtoms, compute_all_virtual_hydrogens, detect_hydrogen_bonds, extract_backbone_atoms,
};

/// Classification of hydrogen bond types based on sequence separation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HBondType {
    /// Intra-helical H-bond (i → i+3, i+4, or i+5)
    IntraHelical,
    /// H-bond between β-strands
    BetaSheet,
    /// Short-range turn H-bond (i → i+2)
    Turn,
    /// Long-range H-bond (sequence separation > 5)
    LongRange,
    /// Inter-chain H-bond
    InterChain,
}

impl HBondType {
    /// Returns true if this is a secondary structure-forming H-bond.
    pub fn is_secondary_structure(&self) -> bool {
        matches!(
            self,
            HBondType::IntraHelical | HBondType::BetaSheet | HBondType::Turn
        )
    }
}

/// A mainchain hydrogen bond with detailed information.
#[derive(Debug, Clone)]
pub struct MainchainHBond {
    /// Chain ID of the donor residue (N-H)
    pub donor_chain: String,
    /// Residue sequence number of the donor
    pub donor_resid: i32,
    /// Residue name of the donor
    pub donor_resname: String,
    /// Insertion code of the donor (if any)
    pub donor_ins_code: Option<char>,
    /// Chain ID of the acceptor residue (C=O)
    pub acceptor_chain: String,
    /// Residue sequence number of the acceptor
    pub acceptor_resid: i32,
    /// Residue name of the acceptor
    pub acceptor_resname: String,
    /// Insertion code of the acceptor (if any)
    pub acceptor_ins_code: Option<char>,
    /// H-bond energy in kcal/mol (more negative = stronger)
    pub energy: f64,
    /// N...O distance in Angstroms
    pub n_o_distance: f64,
    /// Sequence separation (acceptor_resid - donor_resid)
    pub sequence_separation: i32,
    /// Classification of H-bond type
    pub hbond_type: HBondType,
}

impl MainchainHBond {
    /// Returns true if this is a strong H-bond (E < -1.0 kcal/mol).
    pub fn is_strong(&self) -> bool {
        self.energy < -1.0
    }

    /// Returns true if this is a helical H-bond (i → i+3, i+4, or i+5).
    pub fn is_helical(&self) -> bool {
        matches!(self.hbond_type, HBondType::IntraHelical)
    }

    /// Returns true if this is a β-sheet H-bond.
    pub fn is_beta_sheet(&self) -> bool {
        matches!(self.hbond_type, HBondType::BetaSheet)
    }
}

/// H-bonds donated and accepted by a specific residue.
#[derive(Debug, Clone, Default)]
pub struct ResidueHBonds {
    /// H-bonds where this residue is the donor (N-H)
    pub donated: Vec<MainchainHBond>,
    /// H-bonds where this residue is the acceptor (C=O)
    pub accepted: Vec<MainchainHBond>,
}

impl ResidueHBonds {
    /// Returns the total number of H-bonds involving this residue.
    pub fn total(&self) -> usize {
        self.donated.len() + self.accepted.len()
    }

    /// Returns true if this residue participates in any H-bonds.
    pub fn has_hbonds(&self) -> bool {
        !self.donated.is_empty() || !self.accepted.is_empty()
    }
}

/// Statistics about the hydrogen bond network.
#[derive(Debug, Clone, Default)]
pub struct HBondStats {
    /// Total number of H-bonds
    pub total_hbonds: usize,
    /// Number of intra-helical H-bonds (i → i+3, i+4, i+5)
    pub intra_helical: usize,
    /// Number of β-sheet H-bonds
    pub beta_sheet: usize,
    /// Number of turn H-bonds
    pub turn: usize,
    /// Number of long-range H-bonds
    pub long_range: usize,
    /// Number of inter-chain H-bonds
    pub inter_chain: usize,
    /// Mean H-bond energy (kcal/mol)
    pub mean_energy: f64,
    /// Minimum (most negative) H-bond energy
    pub min_energy: f64,
    /// Maximum (least negative) H-bond energy
    pub max_energy: f64,
    /// Number of residues that donate H-bonds
    pub donor_residues: usize,
    /// Number of residues that accept H-bonds
    pub acceptor_residues: usize,
}

/// Calculate N-O distance between donor and acceptor
#[cfg(feature = "dssp")]
fn calculate_n_o_distance(donor: &BackboneAtoms, acceptor: &BackboneAtoms) -> Option<f64> {
    let n = donor.n?;
    let o = acceptor.o?;
    let dx = n[0] - o[0];
    let dy = n[1] - o[1];
    let dz = n[2] - o[2];
    Some((dx * dx + dy * dy + dz * dz).sqrt())
}

/// Classify H-bond type based on sequence separation and chain
fn classify_hbond(
    donor_chain: &str,
    donor_resid: i32,
    acceptor_chain: &str,
    acceptor_resid: i32,
) -> HBondType {
    // Inter-chain H-bonds
    if donor_chain != acceptor_chain {
        return HBondType::InterChain;
    }

    let sep = (acceptor_resid - donor_resid).abs();

    match sep {
        2 => HBondType::Turn,                 // i → i+2 turn
        3..=5 => HBondType::IntraHelical,     // Helical patterns
        _ if sep > 5 => HBondType::BetaSheet, // Long-range = likely β-sheet
        _ => HBondType::LongRange,
    }
}

#[cfg(feature = "dssp")]
impl PdbStructure {
    /// Returns all mainchain (backbone) hydrogen bonds.
    ///
    /// H-bonds are detected using the DSSP algorithm based on electrostatic
    /// energy calculations. Returns bonds sorted by donor residue.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let hbonds = structure.mainchain_hbonds();
    /// for hb in &hbonds {
    ///     println!("{}->{}: {:.2} kcal/mol",
    ///         hb.donor_resid, hb.acceptor_resid, hb.energy);
    /// }
    /// ```
    pub fn mainchain_hbonds(&self) -> Vec<MainchainHBond> {
        let mut backbone_atoms = extract_backbone_atoms(&self.atoms);
        if backbone_atoms.is_empty() {
            return Vec::new();
        }

        // Compute virtual hydrogens where needed
        compute_all_virtual_hydrogens(&mut backbone_atoms);

        // Detect H-bonds
        let raw_hbonds = detect_hydrogen_bonds(&backbone_atoms);

        // Convert to user-friendly format
        raw_hbonds
            .iter()
            .filter_map(|hb| {
                let donor = &backbone_atoms[hb.donor_index];
                let acceptor = &backbone_atoms[hb.acceptor_index];

                let n_o_distance = calculate_n_o_distance(donor, acceptor)?;

                let hbond_type = classify_hbond(
                    &donor.chain_id,
                    donor.residue_seq,
                    &acceptor.chain_id,
                    acceptor.residue_seq,
                );

                let sequence_separation = if donor.chain_id == acceptor.chain_id {
                    acceptor.residue_seq - donor.residue_seq
                } else {
                    0
                };

                Some(MainchainHBond {
                    donor_chain: donor.chain_id.clone(),
                    donor_resid: donor.residue_seq,
                    donor_resname: donor.residue_name.clone(),
                    donor_ins_code: donor.ins_code,
                    acceptor_chain: acceptor.chain_id.clone(),
                    acceptor_resid: acceptor.residue_seq,
                    acceptor_resname: acceptor.residue_name.clone(),
                    acceptor_ins_code: acceptor.ins_code,
                    energy: hb.energy,
                    n_o_distance,
                    sequence_separation,
                    hbond_type,
                })
            })
            .collect()
    }

    /// Returns H-bonds involving a specific residue.
    ///
    /// # Arguments
    ///
    /// * `chain` - Chain identifier (e.g., "A")
    /// * `resid` - Residue sequence number
    ///
    /// # Returns
    ///
    /// A `ResidueHBonds` struct containing H-bonds where the specified
    /// residue is either the donor (N-H) or acceptor (C=O).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let hbonds = structure.hbonds_for_residue("A", 42);
    /// println!("Residue A42 donates {} H-bonds", hbonds.donated.len());
    /// println!("Residue A42 accepts {} H-bonds", hbonds.accepted.len());
    /// ```
    pub fn hbonds_for_residue(&self, chain: &str, resid: i32) -> ResidueHBonds {
        let all_hbonds = self.mainchain_hbonds();
        let mut result = ResidueHBonds::default();

        for hb in all_hbonds {
            if hb.donor_chain == chain && hb.donor_resid == resid {
                result.donated.push(hb);
            } else if hb.acceptor_chain == chain && hb.acceptor_resid == resid {
                result.accepted.push(hb);
            }
        }

        result
    }

    /// Returns statistics about the hydrogen bond network.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let stats = structure.hbond_statistics();
    /// println!("Total H-bonds: {}", stats.total_hbonds);
    /// println!("Helical: {}", stats.intra_helical);
    /// println!("β-sheet: {}", stats.beta_sheet);
    /// println!("Mean energy: {:.2} kcal/mol", stats.mean_energy);
    /// ```
    pub fn hbond_statistics(&self) -> HBondStats {
        let hbonds = self.mainchain_hbonds();
        let mut stats = HBondStats::default();

        if hbonds.is_empty() {
            return stats;
        }

        stats.total_hbonds = hbonds.len();

        let mut sum_energy = 0.0;
        stats.min_energy = f64::MAX;
        stats.max_energy = f64::MIN;

        let mut donor_set = std::collections::HashSet::new();
        let mut acceptor_set = std::collections::HashSet::new();

        for hb in &hbonds {
            // Count by type
            match hb.hbond_type {
                HBondType::IntraHelical => stats.intra_helical += 1,
                HBondType::BetaSheet => stats.beta_sheet += 1,
                HBondType::Turn => stats.turn += 1,
                HBondType::LongRange => stats.long_range += 1,
                HBondType::InterChain => stats.inter_chain += 1,
            }

            // Energy statistics
            sum_energy += hb.energy;
            if hb.energy < stats.min_energy {
                stats.min_energy = hb.energy;
            }
            if hb.energy > stats.max_energy {
                stats.max_energy = hb.energy;
            }

            // Track unique donors and acceptors
            donor_set.insert((hb.donor_chain.clone(), hb.donor_resid));
            acceptor_set.insert((hb.acceptor_chain.clone(), hb.acceptor_resid));
        }

        stats.mean_energy = sum_energy / stats.total_hbonds as f64;
        stats.donor_residues = donor_set.len();
        stats.acceptor_residues = acceptor_set.len();

        stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hbond_type_methods() {
        assert!(HBondType::IntraHelical.is_secondary_structure());
        assert!(HBondType::BetaSheet.is_secondary_structure());
        assert!(HBondType::Turn.is_secondary_structure());
        assert!(!HBondType::LongRange.is_secondary_structure());
        assert!(!HBondType::InterChain.is_secondary_structure());
    }

    #[test]
    fn test_classify_hbond_helical() {
        // i → i+4 (α-helix)
        let hb_type = classify_hbond("A", 1, "A", 5);
        assert_eq!(hb_type, HBondType::IntraHelical);

        // i → i+3 (3₁₀-helix)
        let hb_type = classify_hbond("A", 1, "A", 4);
        assert_eq!(hb_type, HBondType::IntraHelical);
    }

    #[test]
    fn test_classify_hbond_beta_sheet() {
        // Long-range same chain
        let hb_type = classify_hbond("A", 1, "A", 20);
        assert_eq!(hb_type, HBondType::BetaSheet);
    }

    #[test]
    fn test_classify_hbond_turn() {
        let hb_type = classify_hbond("A", 1, "A", 3);
        assert_eq!(hb_type, HBondType::Turn);
    }

    #[test]
    fn test_classify_hbond_inter_chain() {
        let hb_type = classify_hbond("A", 1, "B", 5);
        assert_eq!(hb_type, HBondType::InterChain);
    }

    #[test]
    fn test_mainchain_hbond_methods() {
        let hb = MainchainHBond {
            donor_chain: "A".to_string(),
            donor_resid: 1,
            donor_resname: "ALA".to_string(),
            donor_ins_code: None,
            acceptor_chain: "A".to_string(),
            acceptor_resid: 5,
            acceptor_resname: "LEU".to_string(),
            acceptor_ins_code: None,
            energy: -1.5,
            n_o_distance: 2.9,
            sequence_separation: 4,
            hbond_type: HBondType::IntraHelical,
        };

        assert!(hb.is_strong());
        assert!(hb.is_helical());
        assert!(!hb.is_beta_sheet());
    }

    #[test]
    fn test_residue_hbonds() {
        let mut res_hbonds = ResidueHBonds::default();
        assert!(!res_hbonds.has_hbonds());
        assert_eq!(res_hbonds.total(), 0);

        res_hbonds.donated.push(MainchainHBond {
            donor_chain: "A".to_string(),
            donor_resid: 1,
            donor_resname: "ALA".to_string(),
            donor_ins_code: None,
            acceptor_chain: "A".to_string(),
            acceptor_resid: 5,
            acceptor_resname: "LEU".to_string(),
            acceptor_ins_code: None,
            energy: -1.5,
            n_o_distance: 2.9,
            sequence_separation: 4,
            hbond_type: HBondType::IntraHelical,
        });

        assert!(res_hbonds.has_hbonds());
        assert_eq!(res_hbonds.total(), 1);
    }

    #[test]
    fn test_hbond_stats_default() {
        let stats = HBondStats::default();
        assert_eq!(stats.total_hbonds, 0);
        assert_eq!(stats.intra_helical, 0);
        assert_eq!(stats.beta_sheet, 0);
    }
}
