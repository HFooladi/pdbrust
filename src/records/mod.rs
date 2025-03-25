//! Data structures for different PDB record types

mod atom;
mod model;
mod seqres;
mod conect;
mod ssbond;
mod remark;

// Re-export record types for easy access
pub use atom::Atom;
pub use model::Model;
pub use seqres::SeqRes;
pub use conect::Conect;
pub use ssbond::SSBond;
pub use remark::Remark;