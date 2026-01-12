//! Python bindings for PDB/mmCIF parsing functions

use crate::error::convert_error;
use crate::structure::PyPdbStructure;
use pyo3::prelude::*;

/// Parse a PDB format file
///
/// Args:
///     path: Path to the PDB file
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     IOError: If file cannot be read
///     ValueError: If file format is invalid
#[pyfunction]
pub fn parse_pdb_file(path: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_pdb_file(path)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse an mmCIF format file
///
/// Args:
///     path: Path to the mmCIF file
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     IOError: If file cannot be read
///     ValueError: If file format is invalid
#[pyfunction]
pub fn parse_mmcif_file(path: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_mmcif_file(path)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse a structure file with automatic format detection
///
/// Automatically detects whether the file is PDB or mmCIF format
/// based on file content.
///
/// Args:
///     path: Path to the structure file
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     IOError: If file cannot be read
///     ValueError: If file format is invalid or unrecognized
#[pyfunction]
pub fn parse_structure_file(path: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_structure_file(path)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse PDB format from a string
///
/// Args:
///     content: PDB format string content
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     ValueError: If format is invalid
#[pyfunction]
pub fn parse_pdb_string(content: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_pdb_string(content)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse mmCIF format from a string
///
/// Args:
///     content: mmCIF format string content
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     ValueError: If format is invalid
#[pyfunction]
pub fn parse_mmcif_string(content: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_mmcif_string(content)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse a gzip-compressed PDB file
///
/// Args:
///     path: Path to the gzip-compressed PDB file (.pdb.gz, .ent.gz)
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     IOError: If file cannot be read
///     ValueError: If file format is invalid
#[cfg(feature = "gzip")]
#[pyfunction]
pub fn parse_gzip_pdb_file(path: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_gzip_pdb_file(path)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse a gzip-compressed mmCIF file
///
/// Args:
///     path: Path to the gzip-compressed mmCIF file (.cif.gz)
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     IOError: If file cannot be read
///     ValueError: If file format is invalid
#[cfg(feature = "gzip")]
#[pyfunction]
pub fn parse_gzip_mmcif_file(path: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_gzip_mmcif_file(path)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Parse a gzip-compressed structure file with automatic format detection
///
/// Args:
///     path: Path to the gzip-compressed file
///
/// Returns:
///     PdbStructure object
///
/// Raises:
///     IOError: If file cannot be read
///     ValueError: If file format is invalid
#[cfg(feature = "gzip")]
#[pyfunction]
pub fn parse_gzip_structure_file(path: &str) -> PyResult<PyPdbStructure> {
    pdbrust::parse_gzip_structure_file(path)
        .map(PyPdbStructure::from)
        .map_err(convert_error)
}

/// Write a structure to a PDB file
///
/// Args:
///     structure: PdbStructure object to write
///     path: Output file path
///
/// Raises:
///     IOError: If file cannot be written
#[pyfunction]
pub fn write_pdb_file(structure: &PyPdbStructure, path: &str) -> PyResult<()> {
    pdbrust::write_pdb_file(&structure.inner, path).map_err(convert_error)
}

/// Write a structure to an mmCIF file
///
/// Args:
///     structure: PdbStructure object to write
///     path: Output file path
///
/// Raises:
///     IOError: If file cannot be written
///
/// Example:
///     >>> import pdbrust
///     >>> structure = pdbrust.parse_pdb_file("protein.pdb")
///     >>> pdbrust.write_mmcif_file(structure, "protein.cif")
#[pyfunction]
pub fn write_mmcif_file(structure: &PyPdbStructure, path: &str) -> PyResult<()> {
    pdbrust::write_mmcif_file(&structure.inner, path).map_err(convert_error)
}

/// Write a structure to mmCIF format string
///
/// Args:
///     structure: PdbStructure object to write
///
/// Returns:
///     String containing mmCIF formatted structure
///
/// Raises:
///     IOError: If write fails
///
/// Example:
///     >>> import pdbrust
///     >>> structure = pdbrust.parse_pdb_file("protein.pdb")
///     >>> mmcif_content = pdbrust.write_mmcif_string(structure)
///     >>> print(mmcif_content[:100])
#[pyfunction]
pub fn write_mmcif_string(structure: &PyPdbStructure) -> PyResult<String> {
    pdbrust::write_mmcif_string(&structure.inner).map_err(convert_error)
}

/// Write a structure to a gzip-compressed mmCIF file
///
/// Args:
///     structure: PdbStructure object to write
///     path: Output file path (typically ending in .cif.gz)
///
/// Raises:
///     IOError: If file cannot be written
///
/// Example:
///     >>> import pdbrust
///     >>> structure = pdbrust.parse_pdb_file("protein.pdb")
///     >>> pdbrust.write_gzip_mmcif_file(structure, "protein.cif.gz")
#[cfg(feature = "gzip")]
#[pyfunction]
pub fn write_gzip_mmcif_file(structure: &PyPdbStructure, path: &str) -> PyResult<()> {
    pdbrust::write_gzip_mmcif_file(&structure.inner, path).map_err(convert_error)
}
