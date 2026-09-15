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

/// Decodes a hybrid-36 number from a fixed-width PDB field (e.g. atom serial
/// numbers above 99,999 or residue numbers above 9,999).
///
/// Decimal values (including blank-padded and negative ones) are parsed as
/// usual. Values beyond the decimal range use base 36: `A0000` = 100,000 for a
/// 5-character field, continuing with lowercase `a0000` after `ZZZZZ`.
pub(crate) fn decode_hybrid36(field: &str, width: u32) -> Result<i32, crate::PdbError> {
    let invalid = || crate::PdbError::ParseError(format!("Invalid hybrid-36 number: {field:?}"));

    let first = field.chars().next().ok_or_else(invalid)?;
    if first == ' ' || first == '-' || first == '+' || first.is_ascii_digit() {
        return field.trim().parse::<i32>().map_err(|_| invalid());
    }
    if field.len() != width as usize {
        return Err(invalid());
    }

    let uppercase = first.is_ascii_uppercase();
    let mut value: i64 = 0;
    for c in field.chars() {
        let digit = match c {
            '0'..='9' => c as i64 - '0' as i64,
            'A'..='Z' if uppercase => c as i64 - 'A' as i64 + 10,
            'a'..='z' if !uppercase => c as i64 - 'a' as i64 + 10,
            _ => return Err(invalid()),
        };
        value = value * 36 + digit;
    }

    let power36 = 36_i64.pow(width - 1);
    let decimal_limit = 10_i64.pow(width);
    let decoded = if uppercase {
        value - 10 * power36 + decimal_limit
    } else {
        value + 16 * power36 + decimal_limit
    };
    i32::try_from(decoded).map_err(|_| invalid())
}

#[cfg(test)]
mod tests {
    use super::decode_hybrid36;

    #[test]
    fn decimal_numbers_decode_unchanged() {
        for (field, width, expected) in [
            ("    1", 5, 1),
            ("99999", 5, 99_999),
            ("-9999", 5, -9_999),
            ("+1234", 5, 1_234),
            ("   1", 4, 1),
            ("9999", 4, 9_999),
            ("-999", 4, -999),
        ] {
            assert_eq!(
                decode_hybrid36(field, width).unwrap(),
                expected,
                "{field:?}"
            );
        }
    }

    #[test]
    fn hybrid36_numbers_continue_after_the_decimal_range() {
        for (field, width, expected) in [
            ("A0000", 5, 100_000),
            ("A0001", 5, 100_001),
            ("A000Z", 5, 100_035),
            ("A0010", 5, 100_036),
            ("ZZZZZ", 5, 43_770_015),
            ("a0000", 5, 43_770_016),
            ("A000", 4, 10_000),
            ("ZZZZ", 4, 1_223_055),
            ("a000", 4, 1_223_056),
        ] {
            assert_eq!(
                decode_hybrid36(field, width).unwrap(),
                expected,
                "{field:?}"
            );
        }
    }

    #[test]
    fn invalid_fields_are_rejected() {
        for (field, width) in [("A00*", 4), ("Aa000", 5), ("", 5), ("A000", 5), ("**", 4)] {
            assert!(decode_hybrid36(field, width).is_err(), "{field:?}");
        }
    }
}
