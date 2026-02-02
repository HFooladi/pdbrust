//! Van der Waals and covalent radii tables.
//!
//! This module provides atomic radii data commonly used in structural biology
//! and drug discovery for clash detection and volume calculations.
//!
//! # Van der Waals Radii
//!
//! The van der Waals radii are based on the Bondi radii (1964) with updates
//! from Alvarez (2013) for elements commonly found in drug-like molecules.
//! These are the standard radii used in programs like PyMOL, Chimera, and RDKit.
//!
//! # Covalent Radii
//!
//! Covalent radii are used for metal coordination distance checks and are
//! based on the compilation by Cordero et al. (2008).
//!
//! # References
//!
//! - Bondi, A. (1964). J. Phys. Chem. 68, 441-451.
//! - Alvarez, S. (2013). Dalton Trans. 42, 8617-8636.
//! - Cordero, B. et al. (2008). Dalton Trans. 2832-2838.

use std::collections::HashMap;
use std::sync::LazyLock;

/// Van der Waals radii in Angstroms (Bondi/Alvarez radii).
///
/// Contains radii for elements commonly found in biological molecules
/// and drug-like compounds.
static VDW_RADII: LazyLock<HashMap<&'static str, f64>> = LazyLock::new(|| {
    let mut m = HashMap::new();

    // First row (H)
    m.insert("H", 1.20);

    // Second row (C, N, O, F)
    m.insert("C", 1.70);
    m.insert("N", 1.55);
    m.insert("O", 1.52);
    m.insert("F", 1.47);

    // Third row (Si, P, S, Cl)
    m.insert("SI", 2.10);
    m.insert("P", 1.80);
    m.insert("S", 1.80);
    m.insert("CL", 1.75);

    // Fourth row (As, Se, Br)
    m.insert("AS", 1.85);
    m.insert("SE", 1.90);
    m.insert("BR", 1.85);

    // Fifth row (I)
    m.insert("I", 1.98);

    // Metals commonly found in biological systems
    m.insert("MG", 1.73);
    m.insert("CA", 2.31);
    m.insert("FE", 2.05);
    m.insert("ZN", 1.39);
    m.insert("CU", 1.40);
    m.insert("MN", 2.05);
    m.insert("CO", 2.00);
    m.insert("NI", 1.63);
    m.insert("NA", 2.27);
    m.insert("K", 2.75);

    // Noble gases (rarely in structures but included for completeness)
    m.insert("HE", 1.40);
    m.insert("NE", 1.54);
    m.insert("AR", 1.88);
    m.insert("KR", 2.02);
    m.insert("XE", 2.16);

    // Additional elements found in drugs/ligands
    m.insert("B", 1.92);
    m.insert("AL", 1.84);
    m.insert("GA", 1.87);
    m.insert("GE", 2.11);
    m.insert("SN", 2.17);
    m.insert("PB", 2.02);
    m.insert("BI", 2.07);
    m.insert("SB", 2.06);
    m.insert("TE", 2.06);

    m
});

/// Covalent radii in Angstroms (Cordero et al. 2008).
///
/// Used for metal coordination and bond distance checks.
static COVALENT_RADII: LazyLock<HashMap<&'static str, f64>> = LazyLock::new(|| {
    let mut m = HashMap::new();

    // First row
    m.insert("H", 0.31);

    // Second row
    m.insert("C", 0.76);
    m.insert("N", 0.71);
    m.insert("O", 0.66);
    m.insert("F", 0.57);

    // Third row
    m.insert("SI", 1.11);
    m.insert("P", 1.07);
    m.insert("S", 1.05);
    m.insert("CL", 1.02);

    // Fourth row
    m.insert("AS", 1.19);
    m.insert("SE", 1.20);
    m.insert("BR", 1.20);

    // Fifth row
    m.insert("I", 1.39);

    // Metals
    m.insert("MG", 1.41);
    m.insert("CA", 1.76);
    m.insert("FE", 1.52);
    m.insert("ZN", 1.22);
    m.insert("CU", 1.32);
    m.insert("MN", 1.61);
    m.insert("CO", 1.50);
    m.insert("NI", 1.24);
    m.insert("NA", 1.66);
    m.insert("K", 2.03);

    // Additional elements
    m.insert("B", 0.84);
    m.insert("AL", 1.21);

    m
});

/// Default van der Waals radius for unknown elements (Angstroms).
/// This is a conservative value close to carbon.
pub const DEFAULT_VDW_RADIUS: f64 = 1.70;

/// Default covalent radius for unknown elements (Angstroms).
pub const DEFAULT_COVALENT_RADIUS: f64 = 1.00;

/// Get the van der Waals radius for an element.
///
/// Returns the Bondi/Alvarez radius for known elements, or a conservative
/// default value (1.70 Å) for unknown elements.
///
/// # Arguments
///
/// * `element` - Element symbol (case-insensitive)
///
/// # Returns
///
/// The van der Waals radius in Angstroms.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::ligand_quality::vdw_radius;
///
/// assert_eq!(vdw_radius("C"), 1.70);
/// assert_eq!(vdw_radius("N"), 1.55);
/// assert_eq!(vdw_radius("O"), 1.52);
/// assert_eq!(vdw_radius("h"), 1.20); // case-insensitive
/// ```
pub fn vdw_radius(element: &str) -> f64 {
    let elem_upper = element.trim().to_uppercase();
    *VDW_RADII
        .get(elem_upper.as_str())
        .unwrap_or(&DEFAULT_VDW_RADIUS)
}

/// Get the covalent radius for an element.
///
/// Returns the Cordero radius for known elements, or a conservative
/// default value (1.00 Å) for unknown elements.
///
/// # Arguments
///
/// * `element` - Element symbol (case-insensitive)
///
/// # Returns
///
/// The covalent radius in Angstroms.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::ligand_quality::covalent_radius;
///
/// assert_eq!(covalent_radius("C"), 0.76);
/// assert_eq!(covalent_radius("FE"), 1.52);
/// ```
pub fn covalent_radius(element: &str) -> f64 {
    let elem_upper = element.trim().to_uppercase();
    *COVALENT_RADII
        .get(elem_upper.as_str())
        .unwrap_or(&DEFAULT_COVALENT_RADIUS)
}

/// Check if an element is a metal.
///
/// # Arguments
///
/// * `element` - Element symbol (case-insensitive)
///
/// # Returns
///
/// `true` if the element is a common biological metal.
pub fn is_metal(element: &str) -> bool {
    let elem_upper = element.trim().to_uppercase();
    matches!(
        elem_upper.as_str(),
        "MG" | "CA"
            | "FE"
            | "ZN"
            | "CU"
            | "MN"
            | "CO"
            | "NI"
            | "NA"
            | "K"
            | "AL"
            | "GA"
            | "SN"
            | "PB"
    )
}

/// Calculate the expected minimum distance between two atoms.
///
/// Uses the van der Waals radii with a scaling factor to determine
/// the minimum allowed distance (typically 0.75 × sum of vdW radii
/// per PoseBusters).
///
/// # Arguments
///
/// * `element1` - First element symbol
/// * `element2` - Second element symbol
/// * `scale` - Scaling factor (typically 0.75 for clash detection)
///
/// # Returns
///
/// The minimum expected distance in Angstroms.
pub fn min_contact_distance(element1: &str, element2: &str, scale: f64) -> f64 {
    (vdw_radius(element1) + vdw_radius(element2)) * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vdw_radius_common_elements() {
        assert!((vdw_radius("H") - 1.20).abs() < 1e-10);
        assert!((vdw_radius("C") - 1.70).abs() < 1e-10);
        assert!((vdw_radius("N") - 1.55).abs() < 1e-10);
        assert!((vdw_radius("O") - 1.52).abs() < 1e-10);
        assert!((vdw_radius("S") - 1.80).abs() < 1e-10);
    }

    #[test]
    fn test_vdw_radius_case_insensitive() {
        assert_eq!(vdw_radius("C"), vdw_radius("c"));
        assert_eq!(vdw_radius("CL"), vdw_radius("cl"));
        assert_eq!(vdw_radius("CL"), vdw_radius("Cl"));
    }

    #[test]
    fn test_vdw_radius_unknown_element() {
        assert_eq!(vdw_radius("XX"), DEFAULT_VDW_RADIUS);
        assert_eq!(vdw_radius(""), DEFAULT_VDW_RADIUS);
    }

    #[test]
    fn test_covalent_radius_common_elements() {
        assert!((covalent_radius("C") - 0.76).abs() < 1e-10);
        assert!((covalent_radius("N") - 0.71).abs() < 1e-10);
        assert!((covalent_radius("O") - 0.66).abs() < 1e-10);
        assert!((covalent_radius("FE") - 1.52).abs() < 1e-10);
    }

    #[test]
    fn test_covalent_radius_unknown_element() {
        assert_eq!(covalent_radius("XX"), DEFAULT_COVALENT_RADIUS);
    }

    #[test]
    fn test_is_metal() {
        assert!(is_metal("FE"));
        assert!(is_metal("ZN"));
        assert!(is_metal("MG"));
        assert!(is_metal("fe")); // case-insensitive
        assert!(!is_metal("C"));
        assert!(!is_metal("N"));
        assert!(!is_metal("O"));
    }

    #[test]
    fn test_min_contact_distance() {
        // C-C contact: (1.70 + 1.70) * 0.75 = 2.55
        let dist = min_contact_distance("C", "C", 0.75);
        assert!((dist - 2.55).abs() < 1e-10);

        // C-N contact: (1.70 + 1.55) * 0.75 = 2.4375
        let dist = min_contact_distance("C", "N", 0.75);
        assert!((dist - 2.4375).abs() < 1e-10);
    }

    #[test]
    fn test_metals_have_radii() {
        // All metals should have defined vdW radii
        let metals = ["MG", "CA", "FE", "ZN", "CU", "MN", "CO", "NI", "NA", "K"];
        for metal in metals {
            let r = vdw_radius(metal);
            assert!(
                r > 1.0 && r != DEFAULT_VDW_RADIUS,
                "Metal {} should have a defined radius, got {}",
                metal,
                r
            );
        }
    }

    #[test]
    fn test_halogens_have_radii() {
        let halogens = ["F", "CL", "BR", "I"];
        for hal in halogens {
            let r = vdw_radius(hal);
            assert!(
                r > 1.0 && r != DEFAULT_VDW_RADIUS,
                "Halogen {} should have a defined radius, got {}",
                hal,
                r
            );
        }
    }
}
