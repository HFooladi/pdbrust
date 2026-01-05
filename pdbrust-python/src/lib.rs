//! Python bindings for PDBRust
//!
//! This crate provides Python bindings for the PDBRust library using PyO3.

use pyo3::prelude::*;

mod atom;
mod error;
mod parsing;
mod records;
mod structure;

#[cfg(feature = "descriptors")]
mod descriptors;

#[cfg(feature = "quality")]
mod quality;

#[cfg(feature = "summary")]
mod summary;

#[cfg(feature = "rcsb")]
mod rcsb;

#[cfg(feature = "numpy")]
pub mod numpy_support;

/// PDBRust Python module
///
/// High-performance PDB/mmCIF parsing and analysis library.
#[pymodule]
fn _pdbrust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Core types
    m.add_class::<structure::PyPdbStructure>()?;
    m.add_class::<atom::PyAtom>()?;
    m.add_class::<records::PySSBond>()?;
    m.add_class::<records::PySeqRes>()?;
    m.add_class::<records::PyConect>()?;
    m.add_class::<records::PyRemark>()?;
    m.add_class::<records::PyModel>()?;

    // Parsing functions
    m.add_function(wrap_pyfunction!(parsing::parse_pdb_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::parse_mmcif_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::parse_structure_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::parse_pdb_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::parse_mmcif_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::write_pdb_file, m)?)?;

    // Gzip parsing (feature-gated)
    #[cfg(feature = "gzip")]
    {
        m.add_function(wrap_pyfunction!(parsing::parse_gzip_pdb_file, m)?)?;
        m.add_function(wrap_pyfunction!(parsing::parse_gzip_mmcif_file, m)?)?;
        m.add_function(wrap_pyfunction!(parsing::parse_gzip_structure_file, m)?)?;
    }

    // Descriptors (feature-gated)
    #[cfg(feature = "descriptors")]
    {
        m.add_class::<descriptors::PyStructureDescriptors>()?;
    }

    // Quality (feature-gated)
    #[cfg(feature = "quality")]
    {
        m.add_class::<quality::PyQualityReport>()?;
    }

    // Summary (feature-gated)
    #[cfg(feature = "summary")]
    {
        m.add_class::<summary::PyStructureSummary>()?;
    }

    // RCSB (feature-gated)
    #[cfg(feature = "rcsb")]
    {
        m.add_class::<rcsb::PyFileFormat>()?;
        m.add_class::<rcsb::PyExperimentalMethod>()?;
        m.add_class::<rcsb::PyPolymerType>()?;
        m.add_class::<rcsb::PySearchQuery>()?;
        m.add_class::<rcsb::PySearchResult>()?;
        m.add_function(wrap_pyfunction!(rcsb::rcsb_search, m)?)?;
        m.add_function(wrap_pyfunction!(rcsb::download_structure, m)?)?;
        m.add_function(wrap_pyfunction!(rcsb::download_pdb_string, m)?)?;
        m.add_function(wrap_pyfunction!(rcsb::download_to_file, m)?)?;
    }

    Ok(())
}
