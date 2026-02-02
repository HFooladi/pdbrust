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

#[cfg(feature = "geometry")]
mod geometry;

#[cfg(feature = "dssp")]
mod dssp;

#[cfg(feature = "ligand-quality")]
mod ligand_quality;

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
    m.add_function(wrap_pyfunction!(parsing::write_mmcif_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::write_mmcif_string, m)?)?;

    // Gzip parsing and writing (feature-gated)
    #[cfg(feature = "gzip")]
    {
        m.add_function(wrap_pyfunction!(parsing::parse_gzip_pdb_file, m)?)?;
        m.add_function(wrap_pyfunction!(parsing::parse_gzip_mmcif_file, m)?)?;
        m.add_function(wrap_pyfunction!(parsing::parse_gzip_structure_file, m)?)?;
        m.add_function(wrap_pyfunction!(parsing::write_gzip_mmcif_file, m)?)?;
    }

    // Descriptors (feature-gated)
    #[cfg(feature = "descriptors")]
    {
        m.add_class::<descriptors::PyStructureDescriptors>()?;
        m.add_class::<descriptors::PyResidueBFactor>()?;
        // AlphaFold/pLDDT types
        m.add_class::<descriptors::PyConfidenceCategory>()?;
        m.add_class::<descriptors::PyResiduePlddt>()?;
        // Protein-ligand interaction types
        m.add_class::<descriptors::PyContactResidue>()?;
        m.add_class::<descriptors::PyProteinLigandHBond>()?;
        m.add_class::<descriptors::PySaltBridge>()?;
        m.add_class::<descriptors::PyHydrophobicContact>()?;
        m.add_class::<descriptors::PyLigandInteractionProfile>()?;
        m.add_class::<descriptors::PyBindingSite>()?;
    }

    // Dihedral/H-bond types (requires both descriptors and dssp)
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    {
        m.add_class::<descriptors::PyRamachandranRegion>()?;
        m.add_class::<descriptors::PyResidueDihedrals>()?;
        m.add_class::<descriptors::PyRamachandranStats>()?;
        m.add_class::<descriptors::PyHBondType>()?;
        m.add_class::<descriptors::PyMainchainHBond>()?;
        m.add_class::<descriptors::PyResidueHBonds>()?;
        m.add_class::<descriptors::PyHBondStats>()?;
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

    // RCSB Async (feature-gated)
    #[cfg(feature = "rcsb-async")]
    {
        m.add_class::<rcsb::PyAsyncDownloadOptions>()?;
        m.add_class::<rcsb::PyDownloadResult>()?;
        m.add_function(wrap_pyfunction!(rcsb::download_multiple, m)?)?;
    }

    // Geometry (feature-gated)
    #[cfg(feature = "geometry")]
    {
        m.add_class::<geometry::PyAtomSelection>()?;
        m.add_class::<geometry::PyAlignmentResult>()?;
        m.add_class::<geometry::PyPerResidueRmsd>()?;
        // LDDT types
        m.add_class::<geometry::PyLddtOptions>()?;
        m.add_class::<geometry::PyLddtResult>()?;
        m.add_class::<geometry::PyPerResidueLddt>()?;
    }

    // DSSP (feature-gated)
    #[cfg(feature = "dssp")]
    {
        m.add_class::<dssp::PySecondaryStructure>()?;
        m.add_class::<dssp::PyResidueSSAssignment>()?;
        m.add_class::<dssp::PySecondaryStructureAssignment>()?;
    }

    // Ligand Quality (feature-gated)
    #[cfg(feature = "ligand-quality")]
    {
        m.add_class::<ligand_quality::PyAtomClash>()?;
        m.add_class::<ligand_quality::PyLigandPoseReport>()?;
        m.add_function(wrap_pyfunction!(ligand_quality::vdw_radius, m)?)?;
        m.add_function(wrap_pyfunction!(ligand_quality::covalent_radius, m)?)?;
    }

    Ok(())
}
