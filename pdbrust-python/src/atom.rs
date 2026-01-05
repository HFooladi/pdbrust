//! Python bindings for the Atom struct

use pdbrust::Atom;
use pyo3::prelude::*;

/// Represents an atom record from a PDB file.
///
/// Contains all standard PDB ATOM record fields including position,
/// identification, and thermal factor information.
#[pyclass(name = "Atom")]
#[derive(Clone)]
pub struct PyAtom {
    pub(crate) inner: Atom,
}

#[pymethods]
impl PyAtom {
    /// Atom serial number
    #[getter]
    fn serial(&self) -> i32 {
        self.inner.serial
    }

    /// Atom name (e.g., "CA", "N", "O")
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Residue name (e.g., "ALA", "GLY")
    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    /// Chain identifier (e.g., "A", "B")
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    /// Residue sequence number
    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_seq
    }

    /// X coordinate in Angstroms
    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    /// Y coordinate in Angstroms
    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    /// Z coordinate in Angstroms
    #[getter]
    fn z(&self) -> f64 {
        self.inner.z
    }

    /// Occupancy factor
    #[getter]
    fn occupancy(&self) -> f64 {
        self.inner.occupancy
    }

    /// Temperature factor (B-factor)
    #[getter]
    fn temp_factor(&self) -> f64 {
        self.inner.temp_factor
    }

    /// Element symbol
    #[getter]
    fn element(&self) -> &str {
        &self.inner.element
    }

    /// Alternate location indicator
    #[getter]
    fn alt_loc(&self) -> Option<char> {
        self.inner.alt_loc
    }

    /// Insertion code
    #[getter]
    fn ins_code(&self) -> Option<char> {
        self.inner.ins_code
    }

    /// Get coordinates as (x, y, z) tuple
    ///
    /// Returns:
    ///     Tuple of (x, y, z) coordinates in Angstroms
    fn get_coordinates(&self) -> (f64, f64, f64) {
        self.inner.get_coordinates()
    }

    /// Calculate distance to another atom
    ///
    /// Args:
    ///     other: Another Atom object
    ///
    /// Returns:
    ///     Distance in Angstroms
    fn distance_to(&self, other: &PyAtom) -> f64 {
        self.inner.calculate_distance_to(&other.inner)
    }

    /// Calculate squared distance to another atom (faster, no sqrt)
    ///
    /// Args:
    ///     other: Another Atom object
    ///
    /// Returns:
    ///     Squared distance in Angstroms^2
    fn distance_squared_to(&self, other: &PyAtom) -> f64 {
        self.inner.calculate_distance_squared_to(&other.inner)
    }

    /// Calculate angle between three atoms (self -> atom2 -> atom3)
    ///
    /// Args:
    ///     atom2: Middle atom
    ///     atom3: End atom
    ///
    /// Returns:
    ///     Angle in degrees
    fn angle_to(&self, atom2: &PyAtom, atom3: &PyAtom) -> f64 {
        self.inner
            .calculate_angle_between(&atom2.inner, &atom3.inner)
    }

    /// Check if this is a backbone atom (N, CA, C, O)
    ///
    /// Returns:
    ///     True if backbone atom
    fn is_backbone(&self) -> bool {
        self.inner.is_backbone()
    }

    /// Check if this is a hydrogen atom
    ///
    /// Returns:
    ///     True if hydrogen
    fn is_hydrogen(&self) -> bool {
        self.inner.is_hydrogen()
    }

    /// Check if this is a heavy atom (not hydrogen)
    ///
    /// Returns:
    ///     True if heavy atom
    fn is_heavy_atom(&self) -> bool {
        self.inner.is_heavy_atom()
    }

    /// Get residue identifier as (chain_id, residue_seq, ins_code)
    fn get_residue_id(&self) -> (String, i32, Option<char>) {
        let (chain, seq, ins) = self.inner.get_residue_id();
        (chain.to_string(), seq, ins)
    }

    fn __repr__(&self) -> String {
        format!(
            "Atom(serial={}, name='{}', residue='{}{}', chain='{}', coords=({:.3}, {:.3}, {:.3}))",
            self.inner.serial,
            self.inner.name.trim(),
            self.inner.residue_name,
            self.inner.residue_seq,
            self.inner.chain_id,
            self.inner.x,
            self.inner.y,
            self.inner.z
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{} {} {} {}{}",
            self.inner.name.trim(),
            self.inner.element,
            self.inner.chain_id,
            self.inner.residue_name,
            self.inner.residue_seq
        )
    }

    fn __eq__(&self, other: &PyAtom) -> bool {
        self.inner.serial == other.inner.serial
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.serial.hash(&mut hasher);
        hasher.finish()
    }
}

impl From<Atom> for PyAtom {
    fn from(atom: Atom) -> Self {
        PyAtom { inner: atom }
    }
}

impl From<&Atom> for PyAtom {
    fn from(atom: &Atom) -> Self {
        PyAtom {
            inner: atom.clone(),
        }
    }
}
