//! Evaluator for selection expressions against atoms.

use crate::filter::{is_standard_amino_acid, is_standard_nucleotide};
use crate::records::Atom;

use super::ast::{ComparisonOp, SelectionExpr};

/// Evaluate if an atom matches the selection expression.
pub fn evaluate(expr: &SelectionExpr, atom: &Atom) -> bool {
    match expr {
        SelectionExpr::Chain(id) => {
            // Case-insensitive chain comparison
            atom.chain_id.eq_ignore_ascii_case(id)
        }

        SelectionExpr::Name(name) => {
            // Case-insensitive, trim whitespace
            atom.name.trim().eq_ignore_ascii_case(name)
        }

        SelectionExpr::Resname(name) => {
            // Case-insensitive, trim whitespace
            atom.residue_name.trim().eq_ignore_ascii_case(name)
        }

        SelectionExpr::Resid(num) => atom.residue_seq == *num,

        SelectionExpr::ResidRange { start, end } => {
            atom.residue_seq >= *start && atom.residue_seq <= *end
        }

        SelectionExpr::Element(elem) => {
            // Case-insensitive, trim whitespace
            atom.element.trim().eq_ignore_ascii_case(elem)
        }

        SelectionExpr::Bfactor(op, val) => compare_f64(atom.temp_factor, *op, *val),

        SelectionExpr::Occupancy(op, val) => compare_f64(atom.occupancy, *op, *val),

        SelectionExpr::Backbone => atom.is_backbone(),

        SelectionExpr::Protein => is_standard_amino_acid(&atom.residue_name),

        SelectionExpr::Nucleic => is_standard_nucleotide(&atom.residue_name),

        SelectionExpr::Water => {
            let resname = atom.residue_name.trim().to_uppercase();
            resname == "HOH" || resname == "WAT" || resname == "TIP3" || resname == "TIP"
        }

        SelectionExpr::Hetero => {
            // HETATM records are non-standard residues
            !is_standard_amino_acid(&atom.residue_name)
                && !is_standard_nucleotide(&atom.residue_name)
        }

        SelectionExpr::Hydrogen => atom.is_hydrogen(),

        SelectionExpr::And(left, right) => evaluate(left, atom) && evaluate(right, atom),

        SelectionExpr::Or(left, right) => evaluate(left, atom) || evaluate(right, atom),

        SelectionExpr::Not(inner) => !evaluate(inner, atom),

        SelectionExpr::All => true,
    }
}

#[inline]
fn compare_f64(value: f64, op: ComparisonOp, target: f64) -> bool {
    match op {
        ComparisonOp::Lt => value < target,
        ComparisonOp::Gt => value > target,
        ComparisonOp::Le => value <= target,
        ComparisonOp::Ge => value >= target,
        ComparisonOp::Eq => (value - target).abs() < f64::EPSILON,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_atom(chain: &str, name: &str, resname: &str, resid: i32) -> Atom {
        Atom {
            serial: 1,
            name: name.to_string(),
            alt_loc: None,
            residue_name: resname.to_string(),
            chain_id: chain.to_string(),
            residue_seq: resid,
            ins_code: None,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    #[test]
    fn test_eval_chain() {
        let atom = make_atom("A", "CA", "ALA", 1);
        assert!(evaluate(&SelectionExpr::Chain("A".to_string()), &atom));
        assert!(evaluate(&SelectionExpr::Chain("a".to_string()), &atom)); // case insensitive
        assert!(!evaluate(&SelectionExpr::Chain("B".to_string()), &atom));
    }

    #[test]
    fn test_eval_name() {
        let atom = make_atom("A", "CA", "ALA", 1);
        assert!(evaluate(&SelectionExpr::Name("CA".to_string()), &atom));
        assert!(evaluate(&SelectionExpr::Name("ca".to_string()), &atom)); // case insensitive
        assert!(!evaluate(&SelectionExpr::Name("N".to_string()), &atom));
    }

    #[test]
    fn test_eval_resname() {
        let atom = make_atom("A", "CA", "ALA", 1);
        assert!(evaluate(&SelectionExpr::Resname("ALA".to_string()), &atom));
        assert!(evaluate(&SelectionExpr::Resname("ala".to_string()), &atom)); // case insensitive
        assert!(!evaluate(&SelectionExpr::Resname("GLY".to_string()), &atom));
    }

    #[test]
    fn test_eval_resid() {
        let atom = make_atom("A", "CA", "ALA", 50);
        assert!(evaluate(&SelectionExpr::Resid(50), &atom));
        assert!(!evaluate(&SelectionExpr::Resid(51), &atom));
    }

    #[test]
    fn test_eval_resid_range() {
        let atom = make_atom("A", "CA", "ALA", 50);
        assert!(evaluate(
            &SelectionExpr::ResidRange { start: 1, end: 100 },
            &atom
        ));
        assert!(evaluate(
            &SelectionExpr::ResidRange { start: 50, end: 50 },
            &atom
        ));
        assert!(!evaluate(
            &SelectionExpr::ResidRange { start: 1, end: 49 },
            &atom
        ));
        assert!(!evaluate(
            &SelectionExpr::ResidRange { start: 51, end: 100 },
            &atom
        ));
    }

    #[test]
    fn test_eval_element() {
        let atom = make_atom("A", "CA", "ALA", 1);
        assert!(evaluate(&SelectionExpr::Element("C".to_string()), &atom));
        assert!(evaluate(&SelectionExpr::Element("c".to_string()), &atom)); // case insensitive
        assert!(!evaluate(&SelectionExpr::Element("N".to_string()), &atom));
    }

    #[test]
    fn test_eval_bfactor() {
        let mut atom = make_atom("A", "CA", "ALA", 1);
        atom.temp_factor = 25.0;

        assert!(evaluate(
            &SelectionExpr::Bfactor(ComparisonOp::Lt, 30.0),
            &atom
        ));
        assert!(evaluate(
            &SelectionExpr::Bfactor(ComparisonOp::Le, 25.0),
            &atom
        ));
        assert!(evaluate(
            &SelectionExpr::Bfactor(ComparisonOp::Eq, 25.0),
            &atom
        ));
        assert!(evaluate(
            &SelectionExpr::Bfactor(ComparisonOp::Ge, 25.0),
            &atom
        ));
        assert!(evaluate(
            &SelectionExpr::Bfactor(ComparisonOp::Gt, 20.0),
            &atom
        ));
        assert!(!evaluate(
            &SelectionExpr::Bfactor(ComparisonOp::Gt, 30.0),
            &atom
        ));
    }

    #[test]
    fn test_eval_occupancy() {
        let mut atom = make_atom("A", "CA", "ALA", 1);
        atom.occupancy = 0.75;

        assert!(evaluate(
            &SelectionExpr::Occupancy(ComparisonOp::Ge, 0.5),
            &atom
        ));
        assert!(!evaluate(
            &SelectionExpr::Occupancy(ComparisonOp::Lt, 0.5),
            &atom
        ));
    }

    #[test]
    fn test_eval_backbone() {
        let ca = make_atom("A", "CA", "ALA", 1);
        let n = make_atom("A", "N", "ALA", 1);
        let c = make_atom("A", "C", "ALA", 1);
        let o = make_atom("A", "O", "ALA", 1);
        let cb = make_atom("A", "CB", "ALA", 1);

        assert!(evaluate(&SelectionExpr::Backbone, &ca));
        assert!(evaluate(&SelectionExpr::Backbone, &n));
        assert!(evaluate(&SelectionExpr::Backbone, &c));
        assert!(evaluate(&SelectionExpr::Backbone, &o));
        assert!(!evaluate(&SelectionExpr::Backbone, &cb));
    }

    #[test]
    fn test_eval_protein() {
        let protein_atom = make_atom("A", "CA", "ALA", 1);
        let water_atom = make_atom("A", "O", "HOH", 1);

        assert!(evaluate(&SelectionExpr::Protein, &protein_atom));
        assert!(!evaluate(&SelectionExpr::Protein, &water_atom));
    }

    #[test]
    fn test_eval_water() {
        let water_hoh = make_atom("A", "O", "HOH", 1);
        let water_wat = make_atom("A", "O", "WAT", 1);
        let protein_atom = make_atom("A", "CA", "ALA", 1);

        assert!(evaluate(&SelectionExpr::Water, &water_hoh));
        assert!(evaluate(&SelectionExpr::Water, &water_wat));
        assert!(!evaluate(&SelectionExpr::Water, &protein_atom));
    }

    #[test]
    fn test_eval_hydrogen() {
        let mut h_atom = make_atom("A", "H", "ALA", 1);
        h_atom.element = "H".to_string();
        let ca = make_atom("A", "CA", "ALA", 1);

        assert!(evaluate(&SelectionExpr::Hydrogen, &h_atom));
        assert!(!evaluate(&SelectionExpr::Hydrogen, &ca));
    }

    #[test]
    fn test_eval_and() {
        let atom = make_atom("A", "CA", "ALA", 1);
        let expr = SelectionExpr::And(
            Box::new(SelectionExpr::Chain("A".to_string())),
            Box::new(SelectionExpr::Name("CA".to_string())),
        );
        assert!(evaluate(&expr, &atom));

        let expr_false = SelectionExpr::And(
            Box::new(SelectionExpr::Chain("B".to_string())),
            Box::new(SelectionExpr::Name("CA".to_string())),
        );
        assert!(!evaluate(&expr_false, &atom));
    }

    #[test]
    fn test_eval_or() {
        let atom = make_atom("A", "CA", "ALA", 1);
        let expr = SelectionExpr::Or(
            Box::new(SelectionExpr::Chain("A".to_string())),
            Box::new(SelectionExpr::Chain("B".to_string())),
        );
        assert!(evaluate(&expr, &atom));

        let expr_false = SelectionExpr::Or(
            Box::new(SelectionExpr::Chain("B".to_string())),
            Box::new(SelectionExpr::Chain("C".to_string())),
        );
        assert!(!evaluate(&expr_false, &atom));
    }

    #[test]
    fn test_eval_not() {
        let atom = make_atom("A", "CA", "ALA", 1);
        let expr = SelectionExpr::Not(Box::new(SelectionExpr::Chain("B".to_string())));
        assert!(evaluate(&expr, &atom));

        let expr_false = SelectionExpr::Not(Box::new(SelectionExpr::Chain("A".to_string())));
        assert!(!evaluate(&expr_false, &atom));
    }

    #[test]
    fn test_eval_all() {
        let atom = make_atom("A", "CA", "ALA", 1);
        assert!(evaluate(&SelectionExpr::All, &atom));
    }

    #[test]
    fn test_eval_complex() {
        let atom = make_atom("A", "CA", "ALA", 1);
        // (chain A or chain B) and backbone
        let expr = SelectionExpr::And(
            Box::new(SelectionExpr::Or(
                Box::new(SelectionExpr::Chain("A".to_string())),
                Box::new(SelectionExpr::Chain("B".to_string())),
            )),
            Box::new(SelectionExpr::Backbone),
        );
        assert!(evaluate(&expr, &atom));
    }
}
