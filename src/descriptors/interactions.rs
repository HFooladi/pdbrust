//! Protein-ligand interaction analysis.
//!
//! This module provides functionality for analyzing interactions between
//! proteins and small molecule ligands, useful for drug discovery and
//! binding site characterization.
//!
//! # Interaction Types
//!
//! The following interaction types are detected:
//!
//! - **Hydrogen bonds**: Donor-acceptor pairs within 3.5 Å
//! - **Salt bridges**: Charged groups within 4.0 Å
//! - **Hydrophobic contacts**: Non-polar atoms within 4.0 Å
//!
//! # Example
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein_ligand.pdb")?;
//!
//! // Get binding site residues
//! let binding_site = structure.binding_site("ATP", 5.0);
//! for res in &binding_site.contact_residues {
//!     println!("{}{} {} (min distance: {:.2} Å)",
//!         res.chain_id, res.residue_seq, res.residue_name, res.min_distance);
//! }
//!
//! // Analyze specific ligand interactions
//! if let Some(interactions) = structure.ligand_interactions("ATP") {
//!     println!("Contact residues: {}", interactions.contact_residues.len());
//!     println!("H-bonds: {}", interactions.hydrogen_bonds.len());
//!     println!("Salt bridges: {}", interactions.salt_bridges.len());
//!     println!("Hydrophobic contacts: {}", interactions.hydrophobic_contacts.len());
//! }
//! ```

use crate::PdbStructure;
use crate::records::Atom;
use std::collections::{HashMap, HashSet};

/// Type alias for residue contact tracking: (min_distance, num_contacts, residue_name)
type ResidueContactInfo = (f64, usize, String);

/// Residue in contact with a ligand.
#[derive(Debug, Clone)]
pub struct ContactResidue {
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Residue name (3-letter code)
    pub residue_name: String,
    /// Insertion code (if any)
    pub ins_code: Option<char>,
    /// Minimum distance to any ligand atom (Angstroms)
    pub min_distance: f64,
    /// Number of atom-atom contacts
    pub num_contacts: usize,
}

/// A potential hydrogen bond between protein and ligand.
#[derive(Debug, Clone)]
pub struct ProteinLigandHBond {
    /// Chain ID of protein residue
    pub protein_chain: String,
    /// Residue number of protein residue
    pub protein_resid: i32,
    /// Residue name of protein residue
    pub protein_resname: String,
    /// Atom name in protein
    pub protein_atom: String,
    /// Ligand residue name
    pub ligand_name: String,
    /// Atom name in ligand
    pub ligand_atom: String,
    /// Distance between donor and acceptor (Angstroms)
    pub distance: f64,
    /// True if protein is the donor (N-H...O=ligand)
    pub is_protein_donor: bool,
}

/// A salt bridge between protein and ligand.
#[derive(Debug, Clone)]
pub struct SaltBridge {
    /// Chain ID of protein residue
    pub protein_chain: String,
    /// Residue number of protein residue
    pub protein_resid: i32,
    /// Residue name of protein residue
    pub protein_resname: String,
    /// Charged atom in protein (e.g., NZ for LYS, OD1 for ASP)
    pub protein_atom: String,
    /// Ligand residue name
    pub ligand_name: String,
    /// Charged atom in ligand
    pub ligand_atom: String,
    /// Distance between charged groups (Angstroms)
    pub distance: f64,
    /// True if protein residue is positively charged
    pub protein_positive: bool,
}

/// A hydrophobic contact between protein and ligand.
#[derive(Debug, Clone)]
pub struct HydrophobicContact {
    /// Chain ID of protein residue
    pub protein_chain: String,
    /// Residue number of protein residue
    pub protein_resid: i32,
    /// Residue name of protein residue
    pub protein_resname: String,
    /// Hydrophobic atom in protein
    pub protein_atom: String,
    /// Ligand residue name
    pub ligand_name: String,
    /// Hydrophobic atom in ligand
    pub ligand_atom: String,
    /// Distance between atoms (Angstroms)
    pub distance: f64,
}

/// Complete interaction profile for a ligand.
#[derive(Debug, Clone)]
pub struct LigandInteractionProfile {
    /// Ligand residue name (3-letter code)
    pub ligand_name: String,
    /// Chain where ligand is found
    pub ligand_chain: String,
    /// Residue number of ligand
    pub ligand_resid: i32,
    /// Residues in contact with the ligand
    pub contact_residues: Vec<ContactResidue>,
    /// Hydrogen bonds with protein
    pub hydrogen_bonds: Vec<ProteinLigandHBond>,
    /// Salt bridges with protein
    pub salt_bridges: Vec<SaltBridge>,
    /// Hydrophobic contacts with protein
    pub hydrophobic_contacts: Vec<HydrophobicContact>,
}

impl LigandInteractionProfile {
    /// Returns the total number of interactions detected.
    pub fn total_interactions(&self) -> usize {
        self.hydrogen_bonds.len() + self.salt_bridges.len() + self.hydrophobic_contacts.len()
    }

    /// Returns true if any interactions were detected.
    pub fn has_interactions(&self) -> bool {
        !self.hydrogen_bonds.is_empty()
            || !self.salt_bridges.is_empty()
            || !self.hydrophobic_contacts.is_empty()
    }
}

/// Binding site definition (residues around a ligand).
#[derive(Debug, Clone)]
pub struct BindingSite {
    /// Ligand residue name
    pub ligand_name: String,
    /// Chain where ligand is found
    pub ligand_chain: String,
    /// Residue number of ligand
    pub ligand_resid: i32,
    /// Distance cutoff used for detection (Angstroms)
    pub distance_cutoff: f64,
    /// Residues within the distance cutoff
    pub contact_residues: Vec<ContactResidue>,
}

impl BindingSite {
    /// Returns the number of contact residues.
    pub fn num_residues(&self) -> usize {
        self.contact_residues.len()
    }

    /// Returns residues sorted by minimum distance to ligand.
    pub fn residues_by_distance(&self) -> Vec<&ContactResidue> {
        let mut sorted: Vec<_> = self.contact_residues.iter().collect();
        sorted.sort_by(|a, b| a.min_distance.partial_cmp(&b.min_distance).unwrap());
        sorted
    }
}

// Distance thresholds for interaction detection
const HBOND_DISTANCE_CUTOFF: f64 = 3.5;
const SALT_BRIDGE_CUTOFF: f64 = 4.0;
const HYDROPHOBIC_CUTOFF: f64 = 4.0;

// Common H-bond donor atoms
const HBOND_DONORS: &[&str] = &[
    "N", "NE", "NH1", "NH2", "NZ", "ND1", "ND2", "NE1", "NE2", "OG", "OG1", "OH",
];

// Common H-bond acceptor atoms
const HBOND_ACCEPTORS: &[&str] = &[
    "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2",
];

// Positively charged atoms (basic residues)
const POSITIVE_ATOMS: &[(&str, &str)] = &[
    ("LYS", "NZ"),
    ("ARG", "NH1"),
    ("ARG", "NH2"),
    ("ARG", "NE"),
    ("HIS", "ND1"),
    ("HIS", "NE2"),
];

// Negatively charged atoms (acidic residues)
const NEGATIVE_ATOMS: &[(&str, &str)] = &[
    ("ASP", "OD1"),
    ("ASP", "OD2"),
    ("GLU", "OE1"),
    ("GLU", "OE2"),
];

// Hydrophobic residues
const HYDROPHOBIC_RESIDUES: &[&str] = &["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"];

/// Calculate distance between two atoms
fn atom_distance(a1: &Atom, a2: &Atom) -> f64 {
    let dx = a1.x - a2.x;
    let dy = a1.y - a2.y;
    let dz = a1.z - a2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Check if an atom is a potential H-bond donor
fn is_hbond_donor(atom: &Atom) -> bool {
    HBOND_DONORS.contains(&atom.name.trim())
}

/// Check if an atom is a potential H-bond acceptor
fn is_hbond_acceptor(atom: &Atom) -> bool {
    HBOND_ACCEPTORS.contains(&atom.name.trim())
}

/// Check if an atom is positively charged
fn is_positive_atom(atom: &Atom) -> bool {
    POSITIVE_ATOMS
        .iter()
        .any(|(res, atm)| atom.residue_name == *res && atom.name.trim() == *atm)
}

/// Check if an atom is negatively charged
fn is_negative_atom(atom: &Atom) -> bool {
    NEGATIVE_ATOMS
        .iter()
        .any(|(res, atm)| atom.residue_name == *res && atom.name.trim() == *atm)
}

/// Check if an atom is hydrophobic (carbon in hydrophobic residue)
fn is_hydrophobic_atom(atom: &Atom) -> bool {
    // Carbon atoms in hydrophobic residues
    if atom.element.trim() == "C" && HYDROPHOBIC_RESIDUES.contains(&atom.residue_name.as_str()) {
        return true;
    }
    // Side chain carbons in other residues
    if atom.element.trim() == "C" && !["CA", "C"].contains(&atom.name.trim()) {
        return true;
    }
    false
}

/// Check if a ligand atom might be an H-bond donor (has N or O with H)
fn ligand_might_donate(atom: &Atom) -> bool {
    let elem = atom.element.trim();
    elem == "N" || elem == "O"
}

/// Check if a ligand atom might be an H-bond acceptor
fn ligand_might_accept(atom: &Atom) -> bool {
    let elem = atom.element.trim();
    elem == "N" || elem == "O" || elem == "S"
}

/// Check if a ligand atom might be charged
fn ligand_might_be_charged(atom: &Atom) -> bool {
    let elem = atom.element.trim();
    let name = atom.name.trim();
    // Typical charged groups in ligands
    (elem == "N" && (name.contains("+") || name.starts_with("N")))
        || (elem == "O" && (name.contains("-") || name.starts_with("O")))
}

/// Check if a ligand atom is hydrophobic
fn ligand_is_hydrophobic(atom: &Atom) -> bool {
    atom.element.trim() == "C"
}

impl PdbStructure {
    /// Returns the binding site around a ligand.
    ///
    /// # Arguments
    ///
    /// * `ligand_name` - 3-letter code of the ligand (e.g., "ATP", "HEM")
    /// * `distance_cutoff` - Maximum distance in Angstroms (typically 4.0-6.0)
    ///
    /// # Returns
    ///
    /// A `BindingSite` with all protein residues within the cutoff distance.
    /// Returns `None` if the ligand is not found.
    ///
    /// # Example
    ///
    /// ```ignore
    /// if let Some(site) = structure.binding_site("ATP", 5.0) {
    ///     println!("Found {} residues in binding site", site.num_residues());
    ///     for res in site.residues_by_distance() {
    ///         println!("  {} {}: {:.2} Å", res.chain_id, res.residue_seq, res.min_distance);
    ///     }
    /// }
    /// ```
    pub fn binding_site(&self, ligand_name: &str, distance_cutoff: f64) -> Option<BindingSite> {
        // Find ligand atoms (HETATM records with matching name)
        let ligand_atoms: Vec<&Atom> = self
            .atoms
            .iter()
            .filter(|a| a.residue_name == ligand_name)
            .collect();

        if ligand_atoms.is_empty() {
            return None;
        }

        // Get ligand location info from first atom
        let ligand_chain = ligand_atoms[0].chain_id.clone();
        let ligand_resid = ligand_atoms[0].residue_seq;

        // Find protein atoms (standard amino acids)
        let protein_atoms: Vec<&Atom> = self
            .atoms
            .iter()
            .filter(|a| is_standard_amino_acid(&a.residue_name))
            .collect();

        // Track contacts per residue
        let mut residue_contacts: HashMap<(String, i32, Option<char>), ResidueContactInfo> =
            HashMap::new();

        for prot_atom in &protein_atoms {
            for lig_atom in &ligand_atoms {
                let dist = atom_distance(prot_atom, lig_atom);
                if dist <= distance_cutoff {
                    let key = (
                        prot_atom.chain_id.clone(),
                        prot_atom.residue_seq,
                        prot_atom.ins_code,
                    );
                    let entry = residue_contacts.entry(key).or_insert((
                        f64::MAX,
                        0,
                        prot_atom.residue_name.clone(),
                    ));
                    if dist < entry.0 {
                        entry.0 = dist;
                    }
                    entry.1 += 1;
                }
            }
        }

        let contact_residues: Vec<ContactResidue> = residue_contacts
            .into_iter()
            .map(
                |((chain_id, residue_seq, ins_code), (min_dist, contacts, res_name))| {
                    ContactResidue {
                        chain_id,
                        residue_seq,
                        residue_name: res_name,
                        ins_code,
                        min_distance: min_dist,
                        num_contacts: contacts,
                    }
                },
            )
            .collect();

        Some(BindingSite {
            ligand_name: ligand_name.to_string(),
            ligand_chain,
            ligand_resid,
            distance_cutoff,
            contact_residues,
        })
    }

    /// Analyzes interactions between protein and a specific ligand.
    ///
    /// Detects hydrogen bonds, salt bridges, and hydrophobic contacts.
    ///
    /// # Arguments
    ///
    /// * `ligand_name` - 3-letter code of the ligand
    ///
    /// # Returns
    ///
    /// A `LigandInteractionProfile` with all detected interactions.
    /// Returns `None` if the ligand is not found.
    ///
    /// # Example
    ///
    /// ```ignore
    /// if let Some(profile) = structure.ligand_interactions("ATP") {
    ///     for hb in &profile.hydrogen_bonds {
    ///         println!("H-bond: {}:{} - {}:{} ({:.2} Å)",
    ///             hb.protein_resname, hb.protein_atom,
    ///             hb.ligand_name, hb.ligand_atom,
    ///             hb.distance);
    ///     }
    /// }
    /// ```
    pub fn ligand_interactions(&self, ligand_name: &str) -> Option<LigandInteractionProfile> {
        // Find ligand atoms
        let ligand_atoms: Vec<&Atom> = self
            .atoms
            .iter()
            .filter(|a| a.residue_name == ligand_name)
            .collect();

        if ligand_atoms.is_empty() {
            return None;
        }

        let ligand_chain = ligand_atoms[0].chain_id.clone();
        let ligand_resid = ligand_atoms[0].residue_seq;

        // Get binding site first
        let binding_site = self.binding_site(ligand_name, 6.0)?;

        // Find protein atoms
        let protein_atoms: Vec<&Atom> = self
            .atoms
            .iter()
            .filter(|a| is_standard_amino_acid(&a.residue_name))
            .collect();

        let mut hydrogen_bonds = Vec::new();
        let mut salt_bridges = Vec::new();
        let mut hydrophobic_contacts = Vec::new();

        // Detect interactions
        for prot_atom in &protein_atoms {
            for lig_atom in &ligand_atoms {
                let dist = atom_distance(prot_atom, lig_atom);

                // H-bond detection
                if dist <= HBOND_DISTANCE_CUTOFF {
                    // Protein donor, ligand acceptor
                    if is_hbond_donor(prot_atom) && ligand_might_accept(lig_atom) {
                        hydrogen_bonds.push(ProteinLigandHBond {
                            protein_chain: prot_atom.chain_id.clone(),
                            protein_resid: prot_atom.residue_seq,
                            protein_resname: prot_atom.residue_name.clone(),
                            protein_atom: prot_atom.name.clone(),
                            ligand_name: ligand_name.to_string(),
                            ligand_atom: lig_atom.name.clone(),
                            distance: dist,
                            is_protein_donor: true,
                        });
                    }
                    // Protein acceptor, ligand donor
                    if is_hbond_acceptor(prot_atom) && ligand_might_donate(lig_atom) {
                        hydrogen_bonds.push(ProteinLigandHBond {
                            protein_chain: prot_atom.chain_id.clone(),
                            protein_resid: prot_atom.residue_seq,
                            protein_resname: prot_atom.residue_name.clone(),
                            protein_atom: prot_atom.name.clone(),
                            ligand_name: ligand_name.to_string(),
                            ligand_atom: lig_atom.name.clone(),
                            distance: dist,
                            is_protein_donor: false,
                        });
                    }
                }

                // Salt bridge detection
                if dist <= SALT_BRIDGE_CUTOFF {
                    let prot_positive = is_positive_atom(prot_atom);
                    let prot_negative = is_negative_atom(prot_atom);
                    let lig_charged = ligand_might_be_charged(lig_atom);

                    if (prot_positive || prot_negative) && lig_charged {
                        salt_bridges.push(SaltBridge {
                            protein_chain: prot_atom.chain_id.clone(),
                            protein_resid: prot_atom.residue_seq,
                            protein_resname: prot_atom.residue_name.clone(),
                            protein_atom: prot_atom.name.clone(),
                            ligand_name: ligand_name.to_string(),
                            ligand_atom: lig_atom.name.clone(),
                            distance: dist,
                            protein_positive: prot_positive,
                        });
                    }
                }

                // Hydrophobic contact detection
                if dist <= HYDROPHOBIC_CUTOFF
                    && is_hydrophobic_atom(prot_atom)
                    && ligand_is_hydrophobic(lig_atom)
                {
                    hydrophobic_contacts.push(HydrophobicContact {
                        protein_chain: prot_atom.chain_id.clone(),
                        protein_resid: prot_atom.residue_seq,
                        protein_resname: prot_atom.residue_name.clone(),
                        protein_atom: prot_atom.name.clone(),
                        ligand_name: ligand_name.to_string(),
                        ligand_atom: lig_atom.name.clone(),
                        distance: dist,
                    });
                }
            }
        }

        Some(LigandInteractionProfile {
            ligand_name: ligand_name.to_string(),
            ligand_chain,
            ligand_resid,
            contact_residues: binding_site.contact_residues,
            hydrogen_bonds,
            salt_bridges,
            hydrophobic_contacts,
        })
    }

    /// Returns interaction profiles for all ligands in the structure.
    ///
    /// A ligand is identified as a non-standard residue (HETATM) that is
    /// not water (HOH/WAT) and has contacts with protein.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let all_interactions = structure.all_ligand_interactions();
    /// for profile in &all_interactions {
    ///     println!("Ligand {}: {} contacts, {} H-bonds",
    ///         profile.ligand_name,
    ///         profile.contact_residues.len(),
    ///         profile.hydrogen_bonds.len());
    /// }
    /// ```
    pub fn all_ligand_interactions(&self) -> Vec<LigandInteractionProfile> {
        // Find unique ligand names (excluding water)
        let ligand_names: HashSet<String> = self
            .atoms
            .iter()
            .filter(|a| !is_standard_amino_acid(&a.residue_name))
            .filter(|a| !["HOH", "WAT", "H2O", "DOD"].contains(&a.residue_name.as_str()))
            .map(|a| a.residue_name.clone())
            .collect();

        ligand_names
            .into_iter()
            .filter_map(|name| self.ligand_interactions(&name))
            .filter(|profile| !profile.contact_residues.is_empty())
            .collect()
    }
}

/// Check if a residue name is a standard amino acid
fn is_standard_amino_acid(name: &str) -> bool {
    const AMINO_ACIDS: &[&str] = &[
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
    ];
    AMINO_ACIDS.contains(&name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    fn make_atom(
        chain: &str,
        resid: i32,
        res_name: &str,
        name: &str,
        element: &str,
        x: f64,
        y: f64,
        z: f64,
    ) -> Atom {
        Atom {
            serial: resid,
            name: name.to_string(),
            alt_loc: None,
            residue_name: res_name.to_string(),
            chain_id: chain.to_string(),
            residue_seq: resid,
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

    #[test]
    fn test_is_standard_amino_acid() {
        assert!(is_standard_amino_acid("ALA"));
        assert!(is_standard_amino_acid("GLY"));
        assert!(is_standard_amino_acid("TRP"));
        assert!(!is_standard_amino_acid("ATP"));
        assert!(!is_standard_amino_acid("HOH"));
    }

    #[test]
    fn test_atom_distance() {
        let a1 = make_atom("A", 1, "ALA", "CA", "C", 0.0, 0.0, 0.0);
        let a2 = make_atom("A", 2, "GLY", "CA", "C", 3.0, 4.0, 0.0);
        assert!((atom_distance(&a1, &a2) - 5.0).abs() < 0.001);
    }

    #[test]
    fn test_is_hbond_donor() {
        let donor = make_atom("A", 1, "ALA", "N", "N", 0.0, 0.0, 0.0);
        assert!(is_hbond_donor(&donor));

        let non_donor = make_atom("A", 1, "ALA", "CA", "C", 0.0, 0.0, 0.0);
        assert!(!is_hbond_donor(&non_donor));
    }

    #[test]
    fn test_is_positive_atom() {
        let lys_nz = make_atom("A", 1, "LYS", "NZ", "N", 0.0, 0.0, 0.0);
        assert!(is_positive_atom(&lys_nz));

        let ala_n = make_atom("A", 1, "ALA", "N", "N", 0.0, 0.0, 0.0);
        assert!(!is_positive_atom(&ala_n));
    }

    #[test]
    fn test_is_negative_atom() {
        let asp_od1 = make_atom("A", 1, "ASP", "OD1", "O", 0.0, 0.0, 0.0);
        assert!(is_negative_atom(&asp_od1));

        let ala_o = make_atom("A", 1, "ALA", "O", "O", 0.0, 0.0, 0.0);
        assert!(!is_negative_atom(&ala_o));
    }

    #[test]
    fn test_binding_site_not_found() {
        let structure = PdbStructure::new();
        assert!(structure.binding_site("ATP", 5.0).is_none());
    }

    #[test]
    fn test_binding_site_detection() {
        let mut structure = PdbStructure::new();

        // Add protein atoms
        structure
            .atoms
            .push(make_atom("A", 1, "ALA", "CA", "C", 0.0, 0.0, 0.0));
        structure
            .atoms
            .push(make_atom("A", 2, "GLY", "CA", "C", 3.0, 0.0, 0.0));
        structure
            .atoms
            .push(make_atom("A", 3, "VAL", "CA", "C", 20.0, 0.0, 0.0)); // Far away

        // Add ligand atom close to residue 1
        structure
            .atoms
            .push(make_atom("A", 100, "ATP", "C1", "C", 2.0, 0.0, 0.0));

        let site = structure.binding_site("ATP", 5.0);
        assert!(site.is_some());

        let site = site.unwrap();
        assert_eq!(site.ligand_name, "ATP");
        assert!(!site.contact_residues.is_empty());

        // Should find residues 1 and 2, but not 3
        let res_ids: Vec<i32> = site
            .contact_residues
            .iter()
            .map(|r| r.residue_seq)
            .collect();
        assert!(res_ids.contains(&1));
        assert!(res_ids.contains(&2));
        assert!(!res_ids.contains(&3));
    }

    #[test]
    fn test_ligand_interaction_profile() {
        let profile = LigandInteractionProfile {
            ligand_name: "ATP".to_string(),
            ligand_chain: "A".to_string(),
            ligand_resid: 100,
            contact_residues: vec![],
            hydrogen_bonds: vec![ProteinLigandHBond {
                protein_chain: "A".to_string(),
                protein_resid: 1,
                protein_resname: "SER".to_string(),
                protein_atom: "OG".to_string(),
                ligand_name: "ATP".to_string(),
                ligand_atom: "O1".to_string(),
                distance: 2.8,
                is_protein_donor: true,
            }],
            salt_bridges: vec![],
            hydrophobic_contacts: vec![],
        };

        assert_eq!(profile.total_interactions(), 1);
        assert!(profile.has_interactions());
    }

    #[test]
    fn test_binding_site_by_distance() {
        let site = BindingSite {
            ligand_name: "ATP".to_string(),
            ligand_chain: "A".to_string(),
            ligand_resid: 100,
            distance_cutoff: 5.0,
            contact_residues: vec![
                ContactResidue {
                    chain_id: "A".to_string(),
                    residue_seq: 1,
                    residue_name: "ALA".to_string(),
                    ins_code: None,
                    min_distance: 3.5,
                    num_contacts: 2,
                },
                ContactResidue {
                    chain_id: "A".to_string(),
                    residue_seq: 2,
                    residue_name: "GLY".to_string(),
                    ins_code: None,
                    min_distance: 2.5,
                    num_contacts: 1,
                },
            ],
        };

        let sorted = site.residues_by_distance();
        assert_eq!(sorted[0].residue_seq, 2); // Closer
        assert_eq!(sorted[1].residue_seq, 1);
    }
}
