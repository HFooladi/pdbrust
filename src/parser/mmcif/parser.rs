//! mmCIF file parsing implementation

use crate::core::PdbStructure;
use crate::core::mmcif::MmcifParser;
use crate::core::mmcif_converter::mmcif_to_pdb_structure;
use crate::error::PdbError;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parses an mmCIF file and returns a PdbStructure.
///
/// This function reads an mmCIF file and converts it to the same `PdbStructure`
/// format used by the PDB parser, allowing for unified handling of both formats.
///
/// # Arguments
/// * `path` - Path to the mmCIF file to parse
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use pdbrust::parser::parse_mmcif_file;
///
/// let structure = parse_mmcif_file("structure.cif")?;
/// println!("Loaded {} atoms from mmCIF file", structure.atoms.len());
/// ```
pub fn parse_mmcif_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let mut parser = MmcifParser::new();

    // Parse the mmCIF file
    parser
        .parse_file(
            path.as_ref()
                .to_str()
                .ok_or_else(|| PdbError::InvalidRecord("Invalid file path".to_string()))?,
        )
        .map_err(PdbError::IoError)?;

    // Convert to PdbStructure
    mmcif_to_pdb_structure(&parser)
}

/// Parses mmCIF data from a reader and returns a PdbStructure.
///
/// This function is useful when you have mmCIF data from a source other than a file,
/// such as a network stream or embedded data.
///
/// # Arguments
/// * `reader` - Any type implementing BufRead containing mmCIF data
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use std::io::Cursor;
/// use pdbrust::parser::parse_mmcif_reader;
///
/// let mmcif_data = "data_test\n_entry.id TEST\n";
/// let cursor = Cursor::new(mmcif_data);
/// let structure = parse_mmcif_reader(cursor)?;
/// ```
pub fn parse_mmcif_reader<R: BufRead>(reader: R) -> Result<PdbStructure, PdbError> {
    let mut parser = MmcifParser::new();

    // Parse the mmCIF data
    parser.parse_reader(reader).map_err(PdbError::IoError)?;

    // Convert to PdbStructure
    mmcif_to_pdb_structure(&parser)
}

/// Parses mmCIF data from a string and returns a PdbStructure.
///
/// This is a convenience function for parsing mmCIF data that's already in memory.
///
/// # Arguments
/// * `data` - String containing mmCIF data
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use pdbrust::parser::parse_mmcif_string;
///
/// let mmcif_data = r#"
/// data_test
/// _entry.id TEST
/// loop_
/// _atom_site.group_PDB
/// _atom_site.id
/// _atom_site.type_symbol
/// ATOM 1 N
/// "#;
///
/// let structure = parse_mmcif_string(mmcif_data)?;
/// ```
pub fn parse_mmcif_string(data: &str) -> Result<PdbStructure, PdbError> {
    use std::io::Cursor;
    let cursor = Cursor::new(data);
    parse_mmcif_reader(cursor)
}

/// Auto-detects file format and parses accordingly.
///
/// This function examines the file extension and/or content to determine
/// whether it's a PDB or mmCIF file, then parses it appropriately.
///
/// # Arguments
/// * `path` - Path to the structure file to parse
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
///
/// # Examples
/// ```ignore
/// use pdbrust::parser::parse_structure_file;
///
/// // Works with both .pdb and .cif files
/// let structure1 = parse_structure_file("protein.pdb")?;
/// let structure2 = parse_structure_file("protein.cif")?;
/// ```
pub fn parse_structure_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let path_ref = path.as_ref();

    // First, try to determine format from extension
    if let Some(extension) = path_ref.extension() {
        match extension.to_str() {
            Some("cif") | Some("mmcif") => {
                return parse_mmcif_file(path);
            }
            Some("pdb") | Some("ent") => {
                return crate::parser::parse_pdb_file(path);
            }
            _ => {
                // Unknown extension, try content detection
            }
        }
    }

    // Try content-based detection
    detect_and_parse_file(path)
}

/// Detects file format by examining content and parses accordingly.
fn detect_and_parse_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let file = File::open(&path)?;
    let reader = BufReader::new(file);

    // Read the first non-empty line to detect format
    let mut first_meaningful_line = String::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            first_meaningful_line = trimmed.to_string();
            break;
        }
    }

    // Detect format based on content
    if first_meaningful_line.starts_with("data_") || first_meaningful_line.starts_with("_") {
        // Looks like mmCIF format
        parse_mmcif_file(path)
    } else if first_meaningful_line.starts_with("HEADER")
        || first_meaningful_line.starts_with("ATOM")
        || first_meaningful_line.starts_with("HETATM")
        || first_meaningful_line.starts_with("TITLE")
        || first_meaningful_line.starts_with("REMARK")
    {
        // Looks like PDB format
        crate::parser::parse_pdb_file(path)
    } else {
        // Default to PDB format and let the parser handle errors
        crate::parser::parse_pdb_file(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_mmcif_string() {
        let mmcif_data = r#"
data_test
_entry.id TEST_ENTRY
_struct.title "Test Structure"
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . MET A 1 20.154 10.000 5.000 1.00 25.00
ATOM 2 C CA . MET A 1 21.500 10.500 5.500 1.00 24.50
"#;

        let structure = parse_mmcif_string(mmcif_data).unwrap();

        assert_eq!(structure.atoms.len(), 2);
        assert!(structure.header.is_some());
        assert_eq!(structure.title, Some("Test Structure".to_string()));

        let atom1 = &structure.atoms[0];
        assert_eq!(atom1.serial, 1);
        assert_eq!(atom1.name, "N");
        assert_eq!(atom1.element, "N");
        assert_eq!(atom1.residue_name, "MET");
        assert_eq!(atom1.chain_id, "A");
        assert_eq!(atom1.residue_seq, 1);
    }

    #[test]
    fn test_parse_mmcif_reader() {
        let mmcif_data = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 0.000 0.000 0.000 1.00 20.00
"#;

        let cursor = Cursor::new(mmcif_data);
        let structure = parse_mmcif_reader(cursor).unwrap();

        assert_eq!(structure.atoms.len(), 1);
        assert_eq!(structure.atoms[0].residue_name, "ALA");
    }

    #[test]
    fn test_parse_mmcif_file() {
        let mmcif_data = r#"
data_test
_entry.id TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA . VAL B 2 1.000 2.000 3.000 0.90 30.00
"#;

        // Create a temporary file
        let mut temp_file = NamedTempFile::new().unwrap();
        temp_file.write_all(mmcif_data.as_bytes()).unwrap();

        let structure = parse_mmcif_file(temp_file.path()).unwrap();

        assert_eq!(structure.atoms.len(), 1);
        assert_eq!(structure.atoms[0].residue_name, "VAL");
        assert_eq!(structure.atoms[0].chain_id, "B");
        assert_eq!(structure.atoms[0].residue_seq, 2);
    }

    #[test]
    fn test_auto_detect_mmcif() {
        let mmcif_data = r#"data_test
_entry.id TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . ALA A 1 0.0 0.0 0.0 1.0 20.0"#;

        let mut temp_file = NamedTempFile::with_suffix(".unknown").unwrap();
        temp_file.write_all(mmcif_data.as_bytes()).unwrap();

        let structure = parse_structure_file(temp_file.path()).unwrap();
        assert_eq!(structure.atoms.len(), 1);
    }

    #[test]
    fn test_auto_detect_pdb() {
        let pdb_data =
            "ATOM      1  N   ALA A   1      20.154  16.967  23.486  1.00 25.00           N  \n";

        let mut temp_file = NamedTempFile::with_suffix(".unknown").unwrap();
        temp_file.write_all(pdb_data.as_bytes()).unwrap();

        let structure = parse_structure_file(temp_file.path()).unwrap();
        assert_eq!(structure.atoms.len(), 1);
    }

    #[test]
    fn test_extension_detection() {
        // Test .cif extension
        let mmcif_data = "data_test\n_entry.id TEST\n";
        let mut cif_file = NamedTempFile::with_suffix(".cif").unwrap();
        cif_file.write_all(mmcif_data.as_bytes()).unwrap();

        let structure = parse_structure_file(cif_file.path()).unwrap();
        assert!(structure.header.is_some());

        // Test .pdb extension
        let pdb_data = "HEADER    TEST STRUCTURE\n";
        let mut pdb_file = NamedTempFile::with_suffix(".pdb").unwrap();
        pdb_file.write_all(pdb_data.as_bytes()).unwrap();

        let structure = parse_structure_file(pdb_file.path()).unwrap();
        assert!(structure.header.is_some());
    }
}
