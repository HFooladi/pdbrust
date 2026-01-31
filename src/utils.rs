//! Utility functions and helper methods.

use std::num::ParseFloatError;
use std::num::ParseIntError;

/// Converts a string slice to a float, handling whitespace.
pub fn parse_float(s: &str) -> Result<f64, ParseFloatError> {
    s.trim().parse()
}

/// Converts a string slice to an integer, handling whitespace.
pub fn parse_int<T: std::str::FromStr<Err = ParseIntError>>(s: &str) -> Result<T, ParseIntError> {
    s.trim().parse()
}
