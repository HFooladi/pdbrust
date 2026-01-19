//! PyMOL/VMD-style selection language for atom filtering.
//!
//! This module provides a powerful query language for selecting atoms
//! from PDB structures, inspired by PyMOL, VMD, and MDAnalysis.
//!
//! # Syntax
//!
//! ## Basic Selectors
//!
//! | Selector | Description | Example |
//! |----------|-------------|---------|
//! | `chain X` | Select atoms from chain X | `chain A` |
//! | `name X` | Select atoms named X | `name CA` |
//! | `resname X` | Select residue type X | `resname ALA` |
//! | `resid N` | Select residue number N | `resid 50` |
//! | `resid N:M` | Select residues N through M | `resid 1:100` |
//! | `element X` | Select element type X | `element C` |
//!
//! ## Keywords
//!
//! | Keyword | Description |
//! |---------|-------------|
//! | `backbone` | N, CA, C, O atoms |
//! | `protein` | Standard amino acids |
//! | `nucleic` | Standard nucleotides |
//! | `water` | Water molecules (HOH, WAT) |
//! | `hetero` | HETATM records (non-protein, non-nucleic) |
//! | `hydrogen` | Hydrogen atoms |
//! | `all` or `*` | All atoms |
//!
//! ## Numeric Comparisons
//!
//! | Selector | Description |
//! |----------|-------------|
//! | `bfactor < 30.0` | B-factor less than 30 |
//! | `bfactor >= 20.0` | B-factor at least 20 |
//! | `occupancy > 0.5` | Occupancy greater than 0.5 |
//!
//! ## Boolean Operations
//!
//! - `and` - Both conditions must match
//! - `or` - Either condition matches
//! - `not` - Negation
//! - `()` - Grouping with parentheses
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::{parse_pdb_file, PdbStructure};
//!
//! let structure = parse_pdb_file("protein.pdb")?;
//!
//! // Select CA atoms from chain A
//! let selected = structure.select("chain A and name CA")?;
//!
//! // Select protein atoms with low B-factor
//! let low_bfactor = structure.select("protein and bfactor < 30.0")?;
//!
//! // Complex selection with grouping
//! let complex = structure.select("(chain A or chain B) and backbone and not hydrogen")?;
//!
//! // Select residue range
//! let residues = structure.select("resid 1:50 and name CA")?;
//! ```
//!
//! # Alternative Keyword Forms
//!
//! The selection language supports alternative forms for keywords:
//!
//! - `chain` / `chainid`
//! - `name` / `atomname`
//! - `resname` / `resn`
//! - `resid` / `resi` / `resnum`
//! - `element` / `elem`
//! - `bfactor` / `b` / `tempfactor`
//! - `occupancy` / `occ`
//! - `backbone` / `bb`
//! - `water` / `waters` / `solvent`
//! - `hetero` / `hetatm` / `het`
//! - `hydrogen` / `hydrogens` / `h`

mod ast;
mod error;
mod evaluator;
mod lexer;
mod parser;

pub use ast::{ComparisonOp, SelectionExpr};
pub use error::SelectionError;

use crate::core::PdbStructure;
use lexer::Lexer;
use parser::Parser;

/// Parse a selection string into an AST.
///
/// This function parses the selection string and returns the abstract syntax tree.
/// It does not execute the selection against a structure.
///
/// # Arguments
///
/// * `selection` - A selection string (e.g., "chain A and name CA")
///
/// # Returns
///
/// The parsed AST, or a `SelectionError` if the selection string is invalid.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::filter::selection::parse_selection;
///
/// let ast = parse_selection("chain A and name CA")?;
/// ```
pub fn parse_selection(selection: &str) -> Result<SelectionExpr, SelectionError> {
    let trimmed = selection.trim();
    if trimmed.is_empty() {
        return Err(SelectionError::EmptySelection);
    }

    let tokens = Lexer::new(trimmed).tokenize()?;
    Parser::new(tokens).parse()
}

/// Execute a selection expression against a structure.
///
/// This function filters atoms from a structure based on the given AST.
///
/// # Arguments
///
/// * `structure` - The structure to filter
/// * `expr` - The selection expression AST
///
/// # Returns
///
/// A new `PdbStructure` containing only the selected atoms.
pub fn execute_selection(structure: &PdbStructure, expr: &SelectionExpr) -> PdbStructure {
    let filtered_atoms: Vec<_> = structure
        .atoms
        .iter()
        .filter(|atom| evaluator::evaluate(expr, atom))
        .cloned()
        .collect();

    PdbStructure {
        atoms: filtered_atoms,
        header: structure.header.clone(),
        title: structure.title.clone(),
        seqres: structure.seqres.clone(),
        connects: Vec::new(), // Connectivity may be invalid after filtering
        ssbonds: structure.ssbonds.clone(),
        remarks: structure.remarks.clone(),
        models: Vec::new(),
        current_model: structure.current_model,
    }
}

/// Validate a selection string without executing it.
///
/// This function parses the selection string and returns an error if it's invalid.
/// Useful for validation in user interfaces or configuration parsing.
///
/// # Arguments
///
/// * `selection` - A selection string to validate
///
/// # Returns
///
/// `Ok(())` if the selection is valid, or a `SelectionError` if invalid.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::filter::selection::validate_selection;
///
/// assert!(validate_selection("chain A and name CA").is_ok());
/// assert!(validate_selection("invalid syntax (").is_err());
/// ```
pub fn validate_selection(selection: &str) -> Result<(), SelectionError> {
    let _ = parse_selection(selection)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_selection_valid() {
        assert!(parse_selection("chain A").is_ok());
        assert!(parse_selection("chain A and name CA").is_ok());
        assert!(parse_selection("(chain A or chain B) and backbone").is_ok());
    }

    #[test]
    fn test_parse_selection_empty() {
        assert!(matches!(
            parse_selection(""),
            Err(SelectionError::EmptySelection)
        ));
        assert!(matches!(
            parse_selection("   "),
            Err(SelectionError::EmptySelection)
        ));
    }

    #[test]
    fn test_validate_selection() {
        assert!(validate_selection("chain A").is_ok());
        assert!(validate_selection("(chain A").is_err());
        assert!(validate_selection("unknown_keyword").is_err());
    }
}
