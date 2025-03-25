//! REMARK record structure and implementations

/// Represents a REMARK record from a PDB file.
///
/// Contains additional information about the structure.
#[derive(Debug, Clone)]
pub struct Remark {
    /// Remark number identifying the type of information.
    pub number: i32,
    /// Content of the remark.
    pub content: String,
}
