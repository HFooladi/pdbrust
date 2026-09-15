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

/// Encodes a number into a fixed-width PDB field, using hybrid-36 when it does
/// not fit in decimal. Returns `None` if the value is out of range for `width`.
pub(crate) fn encode_hybrid36(value: i32, width: u32) -> Option<String> {
    let field_width = width as usize;
    let value = i64::from(value);
    let decimal_limit = 10_i64.pow(width);
    if value > -10_i64.pow(width - 1) && value < decimal_limit {
        return Some(format!("{value:>field_width$}"));
    }
    if value < 0 {
        return None;
    }

    // Uppercase block first (A000.. to ZZZZ..), then the lowercase block.
    let block = 26 * 36_i64.pow(width - 1);
    let mut offset = value - decimal_limit;
    let digits: &[u8] = if offset < block {
        b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    } else {
        offset -= block;
        if offset >= block {
            return None;
        }
        b"0123456789abcdefghijklmnopqrstuvwxyz"
    };

    let mut n = offset + 10 * 36_i64.pow(width - 1);
    let mut encoded = vec![b'0'; field_width];
    for slot in encoded.iter_mut().rev() {
        *slot = digits[(n % 36) as usize];
        n /= 36;
    }
    Some(encoded.iter().map(|&b| b as char).collect())
}

#[cfg(test)]
mod tests {
    use super::{decode_hybrid36, encode_hybrid36};
    use proptest::prelude::*;

    #[test]
    fn numbers_in_the_decimal_range_are_encoded_in_decimal() {
        for (value, width, expected) in [
            (1, 5, "    1"),
            (99_999, 5, "99999"),
            (-9_999, 5, "-9999"),
            (1, 4, "   1"),
            (9_999, 4, "9999"),
            (-999, 4, "-999"),
        ] {
            assert_eq!(encode_hybrid36(value, width).as_deref(), Some(expected));
        }
    }

    #[test]
    fn larger_numbers_are_encoded_in_hybrid36() {
        for (value, width, expected) in [
            (100_000, 5, "A0000"),
            (100_035, 5, "A000Z"),
            (100_036, 5, "A0010"),
            (43_770_015, 5, "ZZZZZ"),
            (43_770_016, 5, "a0000"),
            (87_440_031, 5, "zzzzz"),
            (10_000, 4, "A000"),
            (1_223_055, 4, "ZZZZ"),
            (1_223_056, 4, "a000"),
            (2_436_111, 4, "zzzz"),
        ] {
            assert_eq!(encode_hybrid36(value, width).as_deref(), Some(expected));
        }
    }

    #[test]
    fn numbers_outside_the_hybrid36_range_are_not_encoded() {
        for (value, width) in [(87_440_032, 5), (-10_000, 5), (2_436_112, 4), (-1_000, 4)] {
            assert_eq!(encode_hybrid36(value, width), None, "{value}");
        }
    }

    proptest! {
        #[test]
        fn encoded_numbers_decode_to_the_original(value in -999..2_436_112i32, wide in any::<bool>()) {
            let width = if wide { 5 } else { 4 };
            let field = encode_hybrid36(value, width).unwrap();
            prop_assert_eq!(field.len(), width as usize);
            prop_assert_eq!(decode_hybrid36(&field, width).unwrap(), value);
        }
    }

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
