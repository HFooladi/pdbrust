//! Canonical molecular classification constants and helpers.
//!
//! This module centralizes the definitions for identifying molecular entity types
//! (water, protein, nucleic acid, ligand, ion) so that all other modules use
//! consistent criteria. Import from here instead of defining local lists.

use crate::records::Atom;

/// Residue names recognized as water across all PDB conventions.
pub const WATER_RESIDUES: &[&str] = &["HOH", "WAT", "H2O", "DOD", "TIP", "TIP3"];

/// The 20 standard amino acid residue names (3-letter codes).
pub const STANDARD_AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Standard nucleotide residue names (RNA + DNA).
pub const STANDARD_NUCLEOTIDES: &[&str] = &[
    "A", "C", "G", "U", // RNA
    "DA", "DC", "DG", "DT", // DNA
];

/// Common ion residue names.
pub const COMMON_IONS: &[&str] = &[
    "NA", "CL", "MG", "ZN", "CA", "FE", "MN", "CO", "CU", "K", "NI", "CD", "HG",
];

/// Check if a residue name is water.
#[inline]
pub fn is_water(residue_name: &str) -> bool {
    let name = residue_name.trim();
    WATER_RESIDUES.contains(&name)
}

/// Check if a residue name is a standard amino acid (the canonical 20).
#[inline]
pub fn is_standard_amino_acid(residue_name: &str) -> bool {
    STANDARD_AMINO_ACIDS.contains(&residue_name.trim())
}

/// Check if a residue name is a standard nucleotide.
#[inline]
pub fn is_standard_nucleotide(residue_name: &str) -> bool {
    STANDARD_NUCLEOTIDES.contains(&residue_name.trim())
}

/// Check if a residue name is a standard biological residue (amino acid or nucleotide).
#[inline]
pub fn is_standard_residue(residue_name: &str) -> bool {
    is_standard_amino_acid(residue_name) || is_standard_nucleotide(residue_name)
}

/// Check if an atom is a protein atom (standard amino acid, not HETATM).
///
/// This uses the `is_hetatm` flag as the primary discriminator, which correctly
/// excludes modified residues like MSE/SEP that are marked as HETATM in PDB files.
#[inline]
pub fn is_protein_atom(atom: &Atom) -> bool {
    !atom.is_hetatm && is_standard_amino_acid(&atom.residue_name)
}

/// Check if an atom is a ligand atom (HETATM, not water, not ion).
#[inline]
pub fn is_ligand_atom(atom: &Atom) -> bool {
    atom.is_hetatm && !is_water(&atom.residue_name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_water() {
        assert!(is_water("HOH"));
        assert!(is_water("WAT"));
        assert!(is_water("H2O"));
        assert!(is_water("DOD"));
        assert!(is_water("TIP"));
        assert!(is_water("TIP3"));
        assert!(!is_water("ALA"));
        assert!(!is_water("LIG"));
    }

    #[test]
    fn test_is_standard_amino_acid() {
        assert!(is_standard_amino_acid("ALA"));
        assert!(is_standard_amino_acid("GLY"));
        assert!(is_standard_amino_acid(" ALA "));
        assert!(!is_standard_amino_acid("HOH"));
        assert!(!is_standard_amino_acid("ATP"));
        // SEC and PYL are non-standard (not in the canonical 20)
        assert!(!is_standard_amino_acid("SEC"));
        assert!(!is_standard_amino_acid("PYL"));
    }

    #[test]
    fn test_is_standard_nucleotide() {
        assert!(is_standard_nucleotide("A"));
        assert!(is_standard_nucleotide("DA"));
        assert!(!is_standard_nucleotide("ALA"));
    }

    #[test]
    fn test_is_standard_residue() {
        assert!(is_standard_residue("ALA"));
        assert!(is_standard_residue("A"));
        assert!(!is_standard_residue("HOH"));
    }

    #[test]
    fn test_is_protein_atom() {
        let protein = Atom {
            serial: 1,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            ins_code: None,
            is_hetatm: false,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        };
        assert!(is_protein_atom(&protein));

        let hetatm = Atom {
            is_hetatm: true,
            residue_name: "MSE".to_string(),
            ..protein.clone()
        };
        assert!(!is_protein_atom(&hetatm));
    }

    #[test]
    fn test_is_ligand_atom() {
        let ligand = Atom {
            serial: 1,
            name: "C1".to_string(),
            alt_loc: None,
            residue_name: "LIG".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 100,
            ins_code: None,
            is_hetatm: true,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        };
        assert!(is_ligand_atom(&ligand));

        let water = Atom {
            residue_name: "HOH".to_string(),
            ..ligand.clone()
        };
        assert!(!is_ligand_atom(&water));
    }
}
