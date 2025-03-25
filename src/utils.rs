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

/// Formats a float to a fixed-width string with 3 decimal places.
pub fn format_float(f: f64) -> String {
    format!("{:.3}", f)
}

/// Formats an integer to a fixed-width string.
pub fn format_int(i: i32) -> String {
    format!("{:>5}", i)
}

/// Formats a string to a fixed width, padding with spaces if needed.
pub fn format_string(s: &str, width: usize) -> String {
    format!("{:width$}", s, width = width)
} 