//! Gzip decompression support for PDB/mmCIF files.
//!
//! This module provides transparent decompression of gzip-compressed
//! structure files, commonly used in the PDB archive (`.ent.gz` files).
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::parser::gzip::{parse_gzip_pdb_file, parse_gzip_structure_file};
//!
//! // Parse a gzip-compressed PDB file from the PDB archive
//! let structure = parse_gzip_pdb_file("pdb1ubq.ent.gz")?;
//!
//! // Auto-detect format within gzip file
//! let structure = parse_gzip_structure_file("structure.pdb.gz")?;
//! ```

use crate::core::PdbStructure;
use crate::error::PdbError;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// Parse a gzip-compressed PDB file.
///
/// Supports files with `.gz` extension or `.ent.gz` naming convention
/// commonly used in the PDB archive.
///
/// # Arguments
/// * `path` - Path to the gzip-compressed PDB file
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use pdbrust::parser::gzip::parse_gzip_pdb_file;
///
/// let structure = parse_gzip_pdb_file("pdb1ubq.ent.gz")?;
/// println!("Loaded {} atoms", structure.atoms.len());
/// ```
pub fn parse_gzip_pdb_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let file = File::open(&path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    crate::parser::pdb::parse_pdb_reader(reader)
}

/// Parse a gzip-compressed mmCIF file.
///
/// # Arguments
/// * `path` - Path to the gzip-compressed mmCIF file
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use pdbrust::parser::gzip::parse_gzip_mmcif_file;
///
/// let structure = parse_gzip_mmcif_file("1ubq.cif.gz")?;
/// println!("Loaded {} atoms", structure.atoms.len());
/// ```
pub fn parse_gzip_mmcif_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let file = File::open(&path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    crate::parser::mmcif::parse_mmcif_reader(reader)
}

/// Parse a gzip-compressed structure file with auto-detection.
///
/// Automatically detects whether the compressed content is PDB or mmCIF format
/// based on the filename or content inspection.
///
/// # Detection Logic
///
/// 1. First checks filename patterns:
///    - `.ent.gz`, `.pdb.gz` → PDB format
///    - `.cif.gz`, `.mmcif.gz` → mmCIF format
/// 2. If extension is ambiguous, inspects first bytes of decompressed content
///
/// # Arguments
/// * `path` - Path to the gzip-compressed structure file
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use pdbrust::parser::gzip::parse_gzip_structure_file;
///
/// // Works with any gzip-compressed structure file
/// let structure = parse_gzip_structure_file("structure.gz")?;
/// ```
pub fn parse_gzip_structure_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let path_ref = path.as_ref();

    // Try to determine format from filename
    if let Some(stem) = path_ref.file_stem() {
        let stem_str = stem.to_string_lossy().to_lowercase();

        // Check for PDB archive convention: pdbXXXX.ent.gz
        if stem_str.ends_with(".ent") || stem_str.ends_with(".pdb") {
            return parse_gzip_pdb_file(path);
        }

        // Check for mmCIF format
        if stem_str.ends_with(".cif") || stem_str.ends_with(".mmcif") {
            return parse_gzip_mmcif_file(path);
        }

        // Check if filename starts with "pdb" (PDB archive convention)
        if stem_str.starts_with("pdb") {
            return parse_gzip_pdb_file(path);
        }
    }

    // Content-based detection: decompress and check first bytes
    let file = File::open(&path)?;
    let mut decoder = GzDecoder::new(file);
    let mut first_bytes = vec![0u8; 512];
    let bytes_read = decoder.read(&mut first_bytes).unwrap_or(0);
    let content = String::from_utf8_lossy(&first_bytes[..bytes_read]);

    // Look for mmCIF indicators
    let first_line = content.lines().find(|line| {
        let trimmed = line.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    });

    let is_mmcif = first_line.map_or(false, |line| {
        line.starts_with("data_") || line.starts_with('_')
    });

    if is_mmcif {
        parse_gzip_mmcif_file(path)
    } else {
        // Default to PDB format
        parse_gzip_pdb_file(path)
    }
}

/// Parse gzip-compressed PDB data from a reader.
///
/// This function is useful when you have gzip-compressed data from a source
/// other than a file, such as a network stream.
///
/// # Arguments
/// * `reader` - Any type implementing Read containing gzip-compressed PDB data
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
pub fn parse_gzip_pdb_reader<R: Read>(reader: R) -> Result<PdbStructure, PdbError> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    crate::parser::pdb::parse_pdb_reader(buf_reader)
}

/// Parse gzip-compressed mmCIF data from a reader.
///
/// # Arguments
/// * `reader` - Any type implementing Read containing gzip-compressed mmCIF data
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
pub fn parse_gzip_mmcif_reader<R: Read>(reader: R) -> Result<PdbStructure, PdbError> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    crate::parser::mmcif::parse_mmcif_reader(buf_reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_gzip_file(content: &str) -> NamedTempFile {
        let file = NamedTempFile::new().unwrap();
        let mut encoder = GzEncoder::new(file.reopen().unwrap(), Compression::default());
        encoder.write_all(content.as_bytes()).unwrap();
        encoder.finish().unwrap();
        file
    }

    #[test]
    fn test_parse_gzip_pdb_simple() {
        let pdb_content = r#"HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
END
"#;
        let file = create_gzip_file(pdb_content);
        let structure = parse_gzip_pdb_file(file.path()).unwrap();
        assert_eq!(structure.atoms.len(), 2);
    }

    #[test]
    fn test_parse_gzip_structure_file_pdb_extension() {
        let pdb_content = r#"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
END
"#;
        let file = NamedTempFile::with_suffix(".pdb.gz").unwrap();
        let mut encoder = GzEncoder::new(file.reopen().unwrap(), Compression::default());
        encoder.write_all(pdb_content.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let structure = parse_gzip_structure_file(file.path()).unwrap();
        assert_eq!(structure.atoms.len(), 1);
    }

    #[test]
    fn test_parse_gzip_reader() {
        use std::io::Cursor;

        let pdb_content = r#"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
END
"#;
        let mut compressed = Vec::new();
        {
            let mut encoder = GzEncoder::new(&mut compressed, Compression::default());
            encoder.write_all(pdb_content.as_bytes()).unwrap();
        }

        let cursor = Cursor::new(compressed);
        let structure = parse_gzip_pdb_reader(cursor).unwrap();
        assert_eq!(structure.atoms.len(), 1);
    }
}
