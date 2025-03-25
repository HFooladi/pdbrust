//! Data structures for different PDB record types

mod atom;
mod conect;
mod model;
mod remark;
mod seqres;
mod ssbond;

// Re-export record types for easy access
pub use atom::Atom;
pub use conect::Conect;
pub use model::Model;
pub use remark::Remark;
pub use seqres::SeqRes;
pub use ssbond::SSBond;
