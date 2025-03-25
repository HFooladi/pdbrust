//! MODEL record structure and implementations

use super::{Atom, Remark};

/// Represents a MODEL record from a PDB file.
///
/// Contains atoms and remarks specific to one model in a multi-model structure.
#[derive(Debug, Clone)]
pub struct Model {
    /// Serial number of the model.
    pub serial: i32,
    /// List of atoms in the model.
    pub atoms: Vec<Atom>,
    /// List of remarks specific to this model.
    pub remarks: Vec<Remark>,
}
