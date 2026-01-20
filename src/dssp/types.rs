//! Type definitions for DSSP secondary structure assignment.
//!
//! This module contains the core types used for representing secondary structure
//! assignments computed by the DSSP algorithm.

use std::fmt;

/// Secondary structure classification based on DSSP 4 (Kabsch & Sander).
///
/// The DSSP algorithm classifies residues into 8 states based on hydrogen
/// bonding patterns and backbone geometry:
///
/// - **H**: α-helix (i → i+4 H-bond pattern)
/// - **G**: 3₁₀-helix (i → i+3 H-bond pattern)
/// - **I**: π-helix (i → i+5 H-bond pattern)
/// - **P**: κ-helix / PPII helix (dihedral-based, no H-bonds)
/// - **E**: Extended strand in β-sheet
/// - **B**: Isolated β-bridge
/// - **T**: Hydrogen-bonded turn
/// - **S**: Bend (high backbone curvature)
/// - **C**: Coil (none of the above)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SecondaryStructure {
    /// α-helix: i → i+4 hydrogen bond pattern
    AlphaHelix,
    /// 3₁₀-helix: i → i+3 hydrogen bond pattern
    Helix310,
    /// π-helix: i → i+5 hydrogen bond pattern
    PiHelix,
    /// κ-helix / PPII helix: polyproline II helix (dihedral-based)
    KappaHelix,
    /// Extended strand: part of β-sheet
    ExtendedStrand,
    /// Isolated β-bridge: single H-bond pair
    BetaBridge,
    /// Turn: H-bonded turn, not part of helix
    Turn,
    /// Bend: high backbone curvature
    Bend,
    /// Coil: none of the above
    #[default]
    Coil,
}

impl SecondaryStructure {
    /// Returns the single-character DSSP code for this secondary structure.
    ///
    /// The codes follow standard DSSP notation:
    /// - H = α-helix
    /// - G = 3₁₀-helix
    /// - I = π-helix
    /// - P = κ-helix / PPII
    /// - E = Extended strand
    /// - B = Isolated β-bridge
    /// - T = Turn
    /// - S = Bend
    /// - C = Coil
    pub fn code(&self) -> char {
        match self {
            SecondaryStructure::AlphaHelix => 'H',
            SecondaryStructure::Helix310 => 'G',
            SecondaryStructure::PiHelix => 'I',
            SecondaryStructure::KappaHelix => 'P',
            SecondaryStructure::ExtendedStrand => 'E',
            SecondaryStructure::BetaBridge => 'B',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Bend => 'S',
            SecondaryStructure::Coil => 'C',
        }
    }

    /// Creates a SecondaryStructure from a DSSP code character.
    ///
    /// Returns `None` if the character is not a valid DSSP code.
    pub fn from_code(code: char) -> Option<Self> {
        match code {
            'H' => Some(SecondaryStructure::AlphaHelix),
            'G' => Some(SecondaryStructure::Helix310),
            'I' => Some(SecondaryStructure::PiHelix),
            'P' => Some(SecondaryStructure::KappaHelix),
            'E' => Some(SecondaryStructure::ExtendedStrand),
            'B' => Some(SecondaryStructure::BetaBridge),
            'T' => Some(SecondaryStructure::Turn),
            'S' => Some(SecondaryStructure::Bend),
            'C' => Some(SecondaryStructure::Coil),
            _ => None,
        }
    }

    /// Returns true if this is a helical secondary structure (H, G, I, or P).
    pub fn is_helix(&self) -> bool {
        matches!(
            self,
            SecondaryStructure::AlphaHelix
                | SecondaryStructure::Helix310
                | SecondaryStructure::PiHelix
                | SecondaryStructure::KappaHelix
        )
    }

    /// Returns true if this is a β-structure (E or B).
    pub fn is_sheet(&self) -> bool {
        matches!(
            self,
            SecondaryStructure::ExtendedStrand | SecondaryStructure::BetaBridge
        )
    }

    /// Returns true if this is coil/loop (C, T, or S).
    pub fn is_coil(&self) -> bool {
        matches!(
            self,
            SecondaryStructure::Coil | SecondaryStructure::Turn | SecondaryStructure::Bend
        )
    }
}

impl fmt::Display for SecondaryStructure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.code())
    }
}

/// Secondary structure assignment for a single residue.
#[derive(Debug, Clone)]
pub struct ResidueSSAssignment {
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Residue name (3-letter code)
    pub residue_name: String,
    /// Insertion code (if any)
    pub ins_code: Option<char>,
    /// Assigned secondary structure
    pub ss: SecondaryStructure,
}

impl ResidueSSAssignment {
    /// Creates a new residue secondary structure assignment.
    pub fn new(
        chain_id: String,
        residue_seq: i32,
        residue_name: String,
        ins_code: Option<char>,
        ss: SecondaryStructure,
    ) -> Self {
        Self {
            chain_id,
            residue_seq,
            residue_name,
            ins_code,
            ss,
        }
    }
}

/// Complete secondary structure assignment for a protein structure.
///
/// Contains per-residue assignments and summary statistics.
#[derive(Debug, Clone)]
pub struct SecondaryStructureAssignment {
    /// Per-residue secondary structure assignments
    pub residue_assignments: Vec<ResidueSSAssignment>,
    /// Number of residues assigned to helix (H, G, I, P)
    pub helix_count: usize,
    /// Number of residues assigned to sheet (E, B)
    pub sheet_count: usize,
    /// Number of residues assigned to coil (C, T, S)
    pub coil_count: usize,
    /// Fraction of residues in helical conformation
    pub helix_fraction: f64,
    /// Fraction of residues in sheet conformation
    pub sheet_fraction: f64,
    /// Fraction of residues in coil conformation
    pub coil_fraction: f64,
    /// Warnings generated during assignment
    pub warnings: Vec<String>,
}

impl SecondaryStructureAssignment {
    /// Creates a new empty secondary structure assignment.
    pub fn new() -> Self {
        Self {
            residue_assignments: Vec::new(),
            helix_count: 0,
            sheet_count: 0,
            coil_count: 0,
            helix_fraction: 0.0,
            sheet_fraction: 0.0,
            coil_fraction: 0.0,
            warnings: Vec::new(),
        }
    }

    /// Returns the secondary structure as a string of single-character codes.
    ///
    /// The string contains one character per residue in sequence order.
    pub fn as_codes(&self) -> String {
        self.residue_assignments
            .iter()
            .map(|r| r.ss.code())
            .collect()
    }

    /// Returns the number of residues with assignments.
    pub fn len(&self) -> usize {
        self.residue_assignments.len()
    }

    /// Returns true if there are no residue assignments.
    pub fn is_empty(&self) -> bool {
        self.residue_assignments.is_empty()
    }

    /// Computes summary statistics from the residue assignments.
    pub fn compute_statistics(&mut self) {
        let total = self.residue_assignments.len();
        if total == 0 {
            return;
        }

        self.helix_count = self
            .residue_assignments
            .iter()
            .filter(|r| r.ss.is_helix())
            .count();
        self.sheet_count = self
            .residue_assignments
            .iter()
            .filter(|r| r.ss.is_sheet())
            .count();
        self.coil_count = self
            .residue_assignments
            .iter()
            .filter(|r| r.ss.is_coil())
            .count();

        let total_f64 = total as f64;
        self.helix_fraction = self.helix_count as f64 / total_f64;
        self.sheet_fraction = self.sheet_count as f64 / total_f64;
        self.coil_fraction = self.coil_count as f64 / total_f64;
    }

    /// Adds a warning message to the assignment.
    pub fn add_warning(&mut self, warning: String) {
        self.warnings.push(warning);
    }

    /// Returns the secondary structure composition as (helix_fraction, sheet_fraction, coil_fraction).
    pub fn composition(&self) -> (f64, f64, f64) {
        (self.helix_fraction, self.sheet_fraction, self.coil_fraction)
    }
}

impl Default for SecondaryStructureAssignment {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for SecondaryStructureAssignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_codes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secondary_structure_codes() {
        assert_eq!(SecondaryStructure::AlphaHelix.code(), 'H');
        assert_eq!(SecondaryStructure::Helix310.code(), 'G');
        assert_eq!(SecondaryStructure::PiHelix.code(), 'I');
        assert_eq!(SecondaryStructure::KappaHelix.code(), 'P');
        assert_eq!(SecondaryStructure::ExtendedStrand.code(), 'E');
        assert_eq!(SecondaryStructure::BetaBridge.code(), 'B');
        assert_eq!(SecondaryStructure::Turn.code(), 'T');
        assert_eq!(SecondaryStructure::Bend.code(), 'S');
        assert_eq!(SecondaryStructure::Coil.code(), 'C');
    }

    #[test]
    fn test_from_code() {
        assert_eq!(
            SecondaryStructure::from_code('H'),
            Some(SecondaryStructure::AlphaHelix)
        );
        assert_eq!(
            SecondaryStructure::from_code('E'),
            Some(SecondaryStructure::ExtendedStrand)
        );
        assert_eq!(
            SecondaryStructure::from_code('P'),
            Some(SecondaryStructure::KappaHelix)
        );
        assert_eq!(SecondaryStructure::from_code('X'), None);
    }

    #[test]
    fn test_is_helix() {
        assert!(SecondaryStructure::AlphaHelix.is_helix());
        assert!(SecondaryStructure::Helix310.is_helix());
        assert!(SecondaryStructure::PiHelix.is_helix());
        assert!(SecondaryStructure::KappaHelix.is_helix());
        assert!(!SecondaryStructure::ExtendedStrand.is_helix());
        assert!(!SecondaryStructure::Coil.is_helix());
    }

    #[test]
    fn test_is_sheet() {
        assert!(SecondaryStructure::ExtendedStrand.is_sheet());
        assert!(SecondaryStructure::BetaBridge.is_sheet());
        assert!(!SecondaryStructure::AlphaHelix.is_sheet());
        assert!(!SecondaryStructure::Coil.is_sheet());
    }

    #[test]
    fn test_is_coil() {
        assert!(SecondaryStructure::Coil.is_coil());
        assert!(SecondaryStructure::Turn.is_coil());
        assert!(SecondaryStructure::Bend.is_coil());
        assert!(!SecondaryStructure::AlphaHelix.is_coil());
        assert!(!SecondaryStructure::ExtendedStrand.is_coil());
    }

    #[test]
    fn test_secondary_structure_assignment() {
        let mut assignment = SecondaryStructureAssignment::new();
        assert!(assignment.is_empty());

        assignment
            .residue_assignments
            .push(ResidueSSAssignment::new(
                "A".to_string(),
                1,
                "ALA".to_string(),
                None,
                SecondaryStructure::AlphaHelix,
            ));
        assignment
            .residue_assignments
            .push(ResidueSSAssignment::new(
                "A".to_string(),
                2,
                "GLY".to_string(),
                None,
                SecondaryStructure::ExtendedStrand,
            ));
        assignment
            .residue_assignments
            .push(ResidueSSAssignment::new(
                "A".to_string(),
                3,
                "VAL".to_string(),
                None,
                SecondaryStructure::Coil,
            ));

        assignment.compute_statistics();

        assert_eq!(assignment.len(), 3);
        assert_eq!(assignment.helix_count, 1);
        assert_eq!(assignment.sheet_count, 1);
        assert_eq!(assignment.coil_count, 1);
        assert!((assignment.helix_fraction - 0.333).abs() < 0.01);
        assert_eq!(assignment.as_codes(), "HEC");
    }

    #[test]
    fn test_default() {
        assert_eq!(SecondaryStructure::default(), SecondaryStructure::Coil);
    }

    #[test]
    fn test_display() {
        let ss = SecondaryStructure::AlphaHelix;
        assert_eq!(format!("{}", ss), "H");
    }
}
