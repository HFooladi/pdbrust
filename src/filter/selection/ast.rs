//! Abstract Syntax Tree types for the selection language.

/// Comparison operators for numeric fields.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ComparisonOp {
    /// Less than (<)
    Lt,
    /// Greater than (>)
    Gt,
    /// Less than or equal (<=)
    Le,
    /// Greater than or equal (>=)
    Ge,
    /// Equal (=)
    Eq,
}

/// A selection expression in the AST.
#[derive(Debug, Clone, PartialEq)]
pub enum SelectionExpr {
    // Property-based selections
    /// Select atoms by chain ID (e.g., "chain A")
    Chain(String),
    /// Select atoms by atom name (e.g., "name CA")
    Name(String),
    /// Select atoms by residue name (e.g., "resname ALA")
    Resname(String),
    /// Select atoms by residue number (e.g., "resid 50")
    Resid(i32),
    /// Select atoms by residue number range (e.g., "resid 1:100")
    ResidRange {
        /// Start of range (inclusive)
        start: i32,
        /// End of range (inclusive)
        end: i32,
    },
    /// Select atoms by element type (e.g., "element C")
    Element(String),

    // Numeric comparisons
    /// Select atoms by B-factor comparison (e.g., "bfactor < 30.0")
    Bfactor(ComparisonOp, f64),
    /// Select atoms by occupancy comparison (e.g., "occupancy >= 0.5")
    Occupancy(ComparisonOp, f64),

    // Pre-defined selections (keywords)
    /// Select backbone atoms (N, CA, C, O)
    Backbone,
    /// Select standard amino acid atoms
    Protein,
    /// Select standard nucleotide atoms
    Nucleic,
    /// Select water molecules (HOH, WAT, TIP3)
    Water,
    /// Select heteroatoms (non-protein, non-nucleic)
    Hetero,
    /// Select hydrogen atoms
    Hydrogen,

    // Boolean combinations
    /// Logical AND of two expressions
    And(Box<SelectionExpr>, Box<SelectionExpr>),
    /// Logical OR of two expressions
    Or(Box<SelectionExpr>, Box<SelectionExpr>),
    /// Logical NOT of an expression
    Not(Box<SelectionExpr>),

    /// Select all atoms
    All,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ast_construction() {
        let expr = SelectionExpr::And(
            Box::new(SelectionExpr::Chain("A".to_string())),
            Box::new(SelectionExpr::Name("CA".to_string())),
        );
        assert!(matches!(expr, SelectionExpr::And(_, _)));
    }

    #[test]
    fn test_resid_range() {
        let expr = SelectionExpr::ResidRange { start: 1, end: 100 };
        if let SelectionExpr::ResidRange { start, end } = expr {
            assert_eq!(start, 1);
            assert_eq!(end, 100);
        } else {
            panic!("Expected ResidRange");
        }
    }

    #[test]
    fn test_comparison_op() {
        let expr = SelectionExpr::Bfactor(ComparisonOp::Lt, 30.0);
        if let SelectionExpr::Bfactor(op, val) = expr {
            assert_eq!(op, ComparisonOp::Lt);
            assert!((val - 30.0).abs() < f64::EPSILON);
        } else {
            panic!("Expected Bfactor");
        }
    }
}
