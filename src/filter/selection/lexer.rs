//! Lexer (tokenizer) for the selection language.

use super::error::SelectionError;

/// Tokens produced by the lexer.
#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    // Keywords
    Chain,
    Name,
    Resname,
    Resid,
    Element,
    Bfactor,
    Occupancy,
    Backbone,
    Protein,
    Nucleic,
    Water,
    Hetero,
    Hydrogen,
    All,

    // Boolean operators
    And,
    Or,
    Not,

    // Comparison operators
    Lt,  // <
    Gt,  // >
    Le,  // <=
    Ge,  // >=
    Eq,  // =

    // Grouping and ranges
    LParen, // (
    RParen, // )
    Colon,  // : for ranges

    // Values
    Identifier(String), // Chain IDs, atom names, etc.
    Integer(i32),       // Residue numbers
    Float(f64),         // B-factor, occupancy values

    // End of input
    Eof,
}

/// A token with its position in the source.
#[derive(Debug, Clone)]
pub struct SpannedToken {
    /// The token
    pub token: Token,
    /// Start position in the input
    pub start: usize,
    /// End position in the input (reserved for future error reporting)
    #[allow(dead_code)]
    pub end: usize,
}

/// Lexer for selection language.
pub struct Lexer<'a> {
    input: &'a str,
    chars: std::iter::Peekable<std::str::CharIndices<'a>>,
    position: usize,
}

impl<'a> Lexer<'a> {
    /// Create a new lexer for the given input.
    pub fn new(input: &'a str) -> Self {
        Self {
            input,
            chars: input.char_indices().peekable(),
            position: 0,
        }
    }

    /// Tokenize the entire input.
    pub fn tokenize(&mut self) -> Result<Vec<SpannedToken>, SelectionError> {
        let mut tokens = Vec::new();
        loop {
            let token = self.next_token()?;
            let is_eof = token.token == Token::Eof;
            tokens.push(token);
            if is_eof {
                break;
            }
        }
        Ok(tokens)
    }

    fn next_token(&mut self) -> Result<SpannedToken, SelectionError> {
        self.skip_whitespace();

        let start = self.position;

        let token = match self.peek_char() {
            None => Token::Eof,
            Some('(') => {
                self.advance();
                Token::LParen
            }
            Some(')') => {
                self.advance();
                Token::RParen
            }
            Some(':') => {
                self.advance();
                Token::Colon
            }
            Some('<') => self.scan_comparison_lt(),
            Some('>') => self.scan_comparison_gt(),
            Some('=') => {
                self.advance();
                // Handle == as well
                if self.peek_char() == Some('=') {
                    self.advance();
                }
                Token::Eq
            }
            Some('!') => {
                self.advance();
                // Support ! as NOT
                Token::Not
            }
            Some('&') => {
                self.advance();
                // Support && as AND
                if self.peek_char() == Some('&') {
                    self.advance();
                }
                Token::And
            }
            Some('|') => {
                self.advance();
                // Support || as OR
                if self.peek_char() == Some('|') {
                    self.advance();
                }
                Token::Or
            }
            Some(c) if c.is_ascii_digit() || c == '-' => self.scan_number()?,
            Some(c) if c.is_ascii_alphabetic() || c == '_' || c == '*' => {
                self.scan_identifier_or_keyword()
            }
            Some(c) => {
                return Err(SelectionError::SyntaxError {
                    message: format!("Unexpected character '{}'", c),
                    position: start,
                })
            }
        };

        Ok(SpannedToken {
            token,
            start,
            end: self.position,
        })
    }

    fn scan_comparison_lt(&mut self) -> Token {
        self.advance(); // consume '<'
        if self.peek_char() == Some('=') {
            self.advance();
            Token::Le
        } else {
            Token::Lt
        }
    }

    fn scan_comparison_gt(&mut self) -> Token {
        self.advance(); // consume '>'
        if self.peek_char() == Some('=') {
            self.advance();
            Token::Ge
        } else {
            Token::Gt
        }
    }

    fn scan_number(&mut self) -> Result<Token, SelectionError> {
        let start = self.position;
        let mut has_dot = false;

        // Handle optional negative sign
        if self.peek_char() == Some('-') {
            self.advance();
        }

        // Scan digits and optional decimal point
        while let Some(c) = self.peek_char() {
            if c.is_ascii_digit() {
                self.advance();
            } else if c == '.' && !has_dot {
                has_dot = true;
                self.advance();
            } else {
                break;
            }
        }

        let text = &self.input[start..self.position];

        if has_dot {
            text.parse::<f64>()
                .map(Token::Float)
                .map_err(|_| SelectionError::InvalidValue {
                    field: "number".to_string(),
                    value: text.to_string(),
                    reason: "invalid float".to_string(),
                })
        } else {
            text.parse::<i32>()
                .map(Token::Integer)
                .map_err(|_| SelectionError::InvalidValue {
                    field: "number".to_string(),
                    value: text.to_string(),
                    reason: "invalid integer".to_string(),
                })
        }
    }

    fn scan_identifier_or_keyword(&mut self) -> Token {
        let start = self.position;

        // Handle * for "all"
        if self.peek_char() == Some('*') {
            self.advance();
            return Token::All;
        }

        while let Some(c) = self.peek_char() {
            if c.is_ascii_alphanumeric() || c == '_' {
                self.advance();
            } else {
                break;
            }
        }

        let text = &self.input[start..self.position];

        // Match keywords (case-insensitive)
        // Note: We avoid single-letter aliases (b, h) that could conflict with chain IDs
        match text.to_lowercase().as_str() {
            "chain" | "chainid" => Token::Chain,
            "name" | "atomname" => Token::Name,
            "resname" | "resn" => Token::Resname,
            "resid" | "resi" | "resnum" => Token::Resid,
            "element" | "elem" => Token::Element,
            "bfactor" | "tempfactor" => Token::Bfactor,
            "occupancy" | "occ" => Token::Occupancy,
            "backbone" | "bb" => Token::Backbone,
            "protein" => Token::Protein,
            "nucleic" => Token::Nucleic,
            "water" | "waters" | "solvent" => Token::Water,
            "hetero" | "hetatm" | "het" => Token::Hetero,
            "hydrogen" | "hydrogens" => Token::Hydrogen,
            "and" => Token::And,
            "or" => Token::Or,
            "not" => Token::Not,
            "all" => Token::All,
            _ => Token::Identifier(text.to_string()),
        }
    }

    fn skip_whitespace(&mut self) {
        while let Some(c) = self.peek_char() {
            if c.is_whitespace() {
                self.advance();
            } else {
                break;
            }
        }
    }

    fn peek_char(&mut self) -> Option<char> {
        self.chars.peek().map(|(_, c)| *c)
    }

    fn advance(&mut self) -> Option<char> {
        if let Some((pos, c)) = self.chars.next() {
            self.position = pos + c.len_utf8();
            Some(c)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokenize(s: &str) -> Vec<Token> {
        Lexer::new(s)
            .tokenize()
            .unwrap()
            .into_iter()
            .map(|t| t.token)
            .collect()
    }

    #[test]
    fn test_simple_chain() {
        let tokens = tokenize("chain A");
        assert_eq!(tokens, vec![Token::Chain, Token::Identifier("A".to_string()), Token::Eof]);
    }

    #[test]
    fn test_boolean_operators() {
        let tokens = tokenize("chain A and name CA");
        assert_eq!(
            tokens,
            vec![
                Token::Chain,
                Token::Identifier("A".to_string()),
                Token::And,
                Token::Name,
                Token::Identifier("CA".to_string()),
                Token::Eof
            ]
        );
    }

    #[test]
    fn test_comparison() {
        let tokens = tokenize("bfactor < 30.0");
        assert_eq!(
            tokens,
            vec![Token::Bfactor, Token::Lt, Token::Float(30.0), Token::Eof]
        );
    }

    #[test]
    fn test_range() {
        let tokens = tokenize("resid 1:100");
        assert_eq!(
            tokens,
            vec![Token::Resid, Token::Integer(1), Token::Colon, Token::Integer(100), Token::Eof]
        );
    }

    #[test]
    fn test_parentheses() {
        let tokens = tokenize("(chain A or chain B)");
        assert_eq!(
            tokens,
            vec![
                Token::LParen,
                Token::Chain,
                Token::Identifier("A".to_string()),
                Token::Or,
                Token::Chain,
                Token::Identifier("B".to_string()),
                Token::RParen,
                Token::Eof
            ]
        );
    }

    #[test]
    fn test_not_operator() {
        let tokens = tokenize("not hydrogen");
        assert_eq!(tokens, vec![Token::Not, Token::Hydrogen, Token::Eof]);
    }

    #[test]
    fn test_keywords_case_insensitive() {
        let tokens = tokenize("CHAIN A AND NAME CA");
        assert_eq!(
            tokens,
            vec![
                Token::Chain,
                Token::Identifier("A".to_string()),
                Token::And,
                Token::Name,
                Token::Identifier("CA".to_string()),
                Token::Eof
            ]
        );
    }

    #[test]
    fn test_all_selector() {
        let tokens = tokenize("all");
        assert_eq!(tokens, vec![Token::All, Token::Eof]);

        let tokens = tokenize("*");
        assert_eq!(tokens, vec![Token::All, Token::Eof]);
    }

    #[test]
    fn test_alternative_keywords() {
        // Test alternate keyword forms
        let tokens = tokenize("chainid A");
        assert_eq!(tokens[0], Token::Chain);

        let tokens = tokenize("resn ALA");
        assert_eq!(tokens[0], Token::Resname);

        let tokens = tokenize("resi 50");
        assert_eq!(tokens[0], Token::Resid);
    }

    #[test]
    fn test_negative_number() {
        let tokens = tokenize("resid -5");
        assert_eq!(tokens, vec![Token::Resid, Token::Integer(-5), Token::Eof]);
    }

    #[test]
    fn test_c_style_operators() {
        let tokens = tokenize("chain A && name CA");
        assert_eq!(
            tokens,
            vec![
                Token::Chain,
                Token::Identifier("A".to_string()),
                Token::And,
                Token::Name,
                Token::Identifier("CA".to_string()),
                Token::Eof
            ]
        );
    }
}
