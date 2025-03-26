//! Data structures for different PDB record types

mod atom;
mod chain;
mod conect;
mod model;
mod remark;
mod residue;
mod seqres;
mod ssbond;

// Re-export record types for easy access
pub use atom::Atom;
pub use chain::Chain;
pub use conect::Conect;
pub use model::Model;
pub use remark::Remark;
pub use residue::Residue;
pub use seqres::SeqRes;
pub use ssbond::SSBond;
