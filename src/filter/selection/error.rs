//! Error types for selection language parsing and evaluation.

use std::fmt;

/// Errors that can occur during selection parsing or evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum SelectionError {
    /// Unexpected token during parsing.
    UnexpectedToken {
        /// What was expected
        expected: String,
        /// What was found
        found: String,
        /// Position in the input string
        position: usize,
    },
    /// Unknown keyword in selection.
    UnknownKeyword {
        /// The unrecognized keyword
        keyword: String,
        /// Position in the input string
        position: usize,
    },
    /// Invalid value for a field.
    InvalidValue {
        /// The field being parsed
        field: String,
        /// The invalid value
        value: String,
        /// Reason for invalidity
        reason: String,
    },
    /// Unclosed parenthesis.
    UnclosedParenthesis {
        /// Position of the opening parenthesis
        position: usize,
    },
    /// Empty selection string.
    EmptySelection,
    /// General syntax error with helpful message.
    SyntaxError {
        /// Error message
        message: String,
        /// Position in the input string
        position: usize,
    },
}

impl fmt::Display for SelectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SelectionError::UnexpectedToken {
                expected,
                found,
                position,
            } => {
                write!(
                    f,
                    "Expected {} at position {}, found '{}'",
                    expected, position, found
                )
            }
            SelectionError::UnknownKeyword { keyword, position } => {
                write!(
                    f,
                    "Unknown selection keyword '{}' at position {}",
                    keyword, position
                )
            }
            SelectionError::InvalidValue {
                field,
                value,
                reason,
            } => {
                write!(
                    f,
                    "Invalid value '{}' for {}: {}",
                    value, field, reason
                )
            }
            SelectionError::UnclosedParenthesis { position } => {
                write!(f, "Unclosed parenthesis at position {}", position)
            }
            SelectionError::EmptySelection => {
                write!(f, "Selection string is empty")
            }
            SelectionError::SyntaxError { message, position } => {
                write!(f, "Syntax error at position {}: {}", position, message)
            }
        }
    }
}

impl std::error::Error for SelectionError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = SelectionError::UnexpectedToken {
            expected: "identifier".to_string(),
            found: "and".to_string(),
            position: 5,
        };
        assert!(err.to_string().contains("Expected identifier"));
        assert!(err.to_string().contains("position 5"));
    }

    #[test]
    fn test_empty_selection_error() {
        let err = SelectionError::EmptySelection;
        assert_eq!(err.to_string(), "Selection string is empty");
    }
}
