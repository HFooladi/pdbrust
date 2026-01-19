//! Recursive descent parser for the selection language.

use super::ast::{ComparisonOp, SelectionExpr};
use super::error::SelectionError;
use super::lexer::{SpannedToken, Token};

/// Recursive descent parser for selection expressions.
pub struct Parser {
    tokens: Vec<SpannedToken>,
    position: usize,
}

impl Parser {
    /// Create a new parser from a token stream.
    pub fn new(tokens: Vec<SpannedToken>) -> Self {
        Self {
            tokens,
            position: 0,
        }
    }

    /// Parse the token stream into an AST.
    pub fn parse(&mut self) -> Result<SelectionExpr, SelectionError> {
        if self.tokens.is_empty() {
            return Err(SelectionError::EmptySelection);
        }

        if self.current_token() == &Token::Eof {
            return Err(SelectionError::EmptySelection);
        }

        let expr = self.parse_or()?;

        // Ensure we consumed all tokens
        if self.current_token() != &Token::Eof {
            return Err(SelectionError::UnexpectedToken {
                expected: "end of selection".to_string(),
                found: format!("{:?}", self.current_token()),
                position: self.current_position(),
            });
        }

        Ok(expr)
    }

    /// Parse OR expressions (lowest precedence).
    fn parse_or(&mut self) -> Result<SelectionExpr, SelectionError> {
        let mut left = self.parse_and()?;

        while self.current_token() == &Token::Or {
            self.advance();
            let right = self.parse_and()?;
            left = SelectionExpr::Or(Box::new(left), Box::new(right));
        }

        Ok(left)
    }

    /// Parse AND expressions (medium precedence).
    fn parse_and(&mut self) -> Result<SelectionExpr, SelectionError> {
        let mut left = self.parse_not()?;

        while self.current_token() == &Token::And {
            self.advance();
            let right = self.parse_not()?;
            left = SelectionExpr::And(Box::new(left), Box::new(right));
        }

        Ok(left)
    }

    /// Parse NOT expressions (highest precedence).
    fn parse_not(&mut self) -> Result<SelectionExpr, SelectionError> {
        if self.current_token() == &Token::Not {
            self.advance();
            let expr = self.parse_not()?; // Right-associative
            return Ok(SelectionExpr::Not(Box::new(expr)));
        }

        self.parse_primary()
    }

    /// Parse primary expressions (atoms, parenthesized, keywords).
    fn parse_primary(&mut self) -> Result<SelectionExpr, SelectionError> {
        let token = self.current_token().clone();
        let position = self.current_position();

        match token {
            Token::LParen => {
                let paren_pos = position;
                self.advance();
                let expr = self.parse_or()?;
                if self.current_token() != &Token::RParen {
                    return Err(SelectionError::UnclosedParenthesis {
                        position: paren_pos,
                    });
                }
                self.advance();
                Ok(expr)
            }
            Token::Chain => {
                self.advance();
                let id = self.expect_identifier()?;
                Ok(SelectionExpr::Chain(id))
            }
            Token::Name => {
                self.advance();
                let name = self.expect_identifier()?;
                Ok(SelectionExpr::Name(name))
            }
            Token::Resname => {
                self.advance();
                let name = self.expect_identifier()?;
                Ok(SelectionExpr::Resname(name))
            }
            Token::Resid => {
                self.advance();
                self.parse_resid()
            }
            Token::Element => {
                self.advance();
                let elem = self.expect_identifier()?;
                Ok(SelectionExpr::Element(elem))
            }
            Token::Bfactor => {
                self.advance();
                let (op, val) = self.parse_comparison()?;
                Ok(SelectionExpr::Bfactor(op, val))
            }
            Token::Occupancy => {
                self.advance();
                let (op, val) = self.parse_comparison()?;
                Ok(SelectionExpr::Occupancy(op, val))
            }
            Token::Backbone => {
                self.advance();
                Ok(SelectionExpr::Backbone)
            }
            Token::Protein => {
                self.advance();
                Ok(SelectionExpr::Protein)
            }
            Token::Nucleic => {
                self.advance();
                Ok(SelectionExpr::Nucleic)
            }
            Token::Water => {
                self.advance();
                Ok(SelectionExpr::Water)
            }
            Token::Hetero => {
                self.advance();
                Ok(SelectionExpr::Hetero)
            }
            Token::Hydrogen => {
                self.advance();
                Ok(SelectionExpr::Hydrogen)
            }
            Token::All => {
                self.advance();
                Ok(SelectionExpr::All)
            }
            Token::Identifier(ref s) => {
                // Check if this is an unknown keyword
                Err(SelectionError::UnknownKeyword {
                    keyword: s.clone(),
                    position,
                })
            }
            _ => Err(SelectionError::UnexpectedToken {
                expected: "selection keyword or '('".to_string(),
                found: format!("{:?}", token),
                position,
            }),
        }
    }

    /// Parse residue ID (single value or range).
    fn parse_resid(&mut self) -> Result<SelectionExpr, SelectionError> {
        let start = self.expect_integer()?;

        if self.current_token() == &Token::Colon {
            self.advance();
            let end = self.expect_integer()?;
            Ok(SelectionExpr::ResidRange { start, end })
        } else {
            Ok(SelectionExpr::Resid(start))
        }
    }

    /// Parse numeric comparison (e.g., < 30.0).
    fn parse_comparison(&mut self) -> Result<(ComparisonOp, f64), SelectionError> {
        let position = self.current_position();
        let op = match self.current_token() {
            Token::Lt => ComparisonOp::Lt,
            Token::Gt => ComparisonOp::Gt,
            Token::Le => ComparisonOp::Le,
            Token::Ge => ComparisonOp::Ge,
            Token::Eq => ComparisonOp::Eq,
            _ => {
                return Err(SelectionError::UnexpectedToken {
                    expected: "comparison operator (<, >, <=, >=, =)".to_string(),
                    found: format!("{:?}", self.current_token()),
                    position,
                });
            }
        };
        self.advance();

        let value = self.expect_number()?;
        Ok((op, value))
    }

    fn expect_identifier(&mut self) -> Result<String, SelectionError> {
        let position = self.current_position();
        match self.current_token().clone() {
            Token::Identifier(s) => {
                self.advance();
                Ok(s)
            }
            // Also accept integers as identifiers (for chain IDs like "1")
            Token::Integer(n) => {
                self.advance();
                Ok(n.to_string())
            }
            other => Err(SelectionError::UnexpectedToken {
                expected: "identifier".to_string(),
                found: format!("{:?}", other),
                position,
            }),
        }
    }

    fn expect_integer(&mut self) -> Result<i32, SelectionError> {
        let position = self.current_position();
        match self.current_token().clone() {
            Token::Integer(n) => {
                self.advance();
                Ok(n)
            }
            other => Err(SelectionError::UnexpectedToken {
                expected: "integer".to_string(),
                found: format!("{:?}", other),
                position,
            }),
        }
    }

    fn expect_number(&mut self) -> Result<f64, SelectionError> {
        let position = self.current_position();
        match self.current_token().clone() {
            Token::Float(f) => {
                self.advance();
                Ok(f)
            }
            Token::Integer(n) => {
                self.advance();
                Ok(n as f64)
            }
            other => Err(SelectionError::UnexpectedToken {
                expected: "number".to_string(),
                found: format!("{:?}", other),
                position,
            }),
        }
    }

    fn current_token(&self) -> &Token {
        self.tokens
            .get(self.position)
            .map(|t| &t.token)
            .unwrap_or(&Token::Eof)
    }

    fn current_position(&self) -> usize {
        self.tokens.get(self.position).map(|t| t.start).unwrap_or(0)
    }

    fn advance(&mut self) {
        if self.position < self.tokens.len() {
            self.position += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::selection::lexer::Lexer;

    fn parse(s: &str) -> Result<SelectionExpr, SelectionError> {
        let tokens = Lexer::new(s).tokenize()?;
        Parser::new(tokens).parse()
    }

    #[test]
    fn test_parse_chain() {
        let expr = parse("chain A").unwrap();
        assert!(matches!(expr, SelectionExpr::Chain(s) if s == "A"));
    }

    #[test]
    fn test_parse_name() {
        let expr = parse("name CA").unwrap();
        assert!(matches!(expr, SelectionExpr::Name(s) if s == "CA"));
    }

    #[test]
    fn test_parse_resname() {
        let expr = parse("resname ALA").unwrap();
        assert!(matches!(expr, SelectionExpr::Resname(s) if s == "ALA"));
    }

    #[test]
    fn test_parse_resid() {
        let expr = parse("resid 50").unwrap();
        assert!(matches!(expr, SelectionExpr::Resid(50)));
    }

    #[test]
    fn test_parse_resid_range() {
        let expr = parse("resid 1:100").unwrap();
        assert!(matches!(
            expr,
            SelectionExpr::ResidRange { start: 1, end: 100 }
        ));
    }

    #[test]
    fn test_parse_and() {
        let expr = parse("chain A and name CA").unwrap();
        assert!(matches!(expr, SelectionExpr::And(_, _)));
    }

    #[test]
    fn test_parse_or() {
        let expr = parse("chain A or chain B").unwrap();
        assert!(matches!(expr, SelectionExpr::Or(_, _)));
    }

    #[test]
    fn test_parse_not() {
        let expr = parse("not hydrogen").unwrap();
        assert!(matches!(expr, SelectionExpr::Not(_)));
    }

    #[test]
    fn test_parse_parentheses() {
        let expr = parse("(chain A or chain B) and backbone").unwrap();
        // Should parse as: (chain A or chain B) and backbone
        if let SelectionExpr::And(left, right) = expr {
            assert!(matches!(*left, SelectionExpr::Or(_, _)));
            assert!(matches!(*right, SelectionExpr::Backbone));
        } else {
            panic!("Expected And expression");
        }
    }

    #[test]
    fn test_parse_precedence_and_or() {
        // OR has lower precedence than AND
        let expr = parse("chain A and name CA or chain B").unwrap();
        // Should parse as: (chain A and name CA) or chain B
        if let SelectionExpr::Or(left, _) = expr {
            assert!(matches!(*left, SelectionExpr::And(_, _)));
        } else {
            panic!("Expected Or expression at top level");
        }
    }

    #[test]
    fn test_parse_precedence_not_and() {
        // NOT binds tighter than AND
        let expr = parse("not hydrogen and protein").unwrap();
        // Should parse as: (not hydrogen) and protein
        if let SelectionExpr::And(left, _) = expr {
            assert!(matches!(*left, SelectionExpr::Not(_)));
        } else {
            panic!("Expected And expression at top level");
        }
    }

    #[test]
    fn test_parse_bfactor() {
        let expr = parse("bfactor < 30.0").unwrap();
        if let SelectionExpr::Bfactor(op, val) = expr {
            assert_eq!(op, ComparisonOp::Lt);
            assert!((val - 30.0).abs() < f64::EPSILON);
        } else {
            panic!("Expected Bfactor expression");
        }
    }

    #[test]
    fn test_parse_complex() {
        let expr = parse("(chain A or chain B) and backbone and not hydrogen").unwrap();
        // Should parse without errors
        assert!(matches!(expr, SelectionExpr::And(_, _)));
    }

    #[test]
    fn test_parse_error_unclosed_paren() {
        let result = parse("(chain A and name CA");
        assert!(matches!(
            result,
            Err(SelectionError::UnclosedParenthesis { .. })
        ));
    }

    #[test]
    fn test_parse_error_empty() {
        let result = parse("");
        assert!(matches!(result, Err(SelectionError::EmptySelection)));
    }

    #[test]
    fn test_parse_error_unknown_keyword() {
        let result = parse("unknown_keyword");
        assert!(matches!(result, Err(SelectionError::UnknownKeyword { .. })));
    }

    #[test]
    fn test_parse_all() {
        let expr = parse("all").unwrap();
        assert!(matches!(expr, SelectionExpr::All));

        let expr = parse("*").unwrap();
        assert!(matches!(expr, SelectionExpr::All));
    }

    #[test]
    fn test_parse_keywords() {
        assert!(matches!(
            parse("backbone").unwrap(),
            SelectionExpr::Backbone
        ));
        assert!(matches!(parse("protein").unwrap(), SelectionExpr::Protein));
        assert!(matches!(parse("nucleic").unwrap(), SelectionExpr::Nucleic));
        assert!(matches!(parse("water").unwrap(), SelectionExpr::Water));
        assert!(matches!(parse("hetero").unwrap(), SelectionExpr::Hetero));
        assert!(matches!(
            parse("hydrogen").unwrap(),
            SelectionExpr::Hydrogen
        ));
    }
}
