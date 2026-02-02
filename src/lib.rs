//! # PDBRust
//!
//! A high-performance Rust library for parsing and analyzing Protein Data Bank (PDB)
//! and mmCIF structure files.
//!
//! PDBRust provides comprehensive tools for working with molecular structure data,
//! from simple file parsing to advanced structural analysis. It's designed for
//! bioinformatics pipelines, structural biology research, and machine learning
//! applications that work with protein structures.
//!
//! For detailed documentation, examples, and best practices, see the [guide] module.
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use pdbrust::{parse_pdb_file, PdbStructure};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Parse a PDB file
//!     let structure = parse_pdb_file("protein.pdb")?;
//!
//!     println!("Atoms: {}", structure.atoms.len());
//!     println!("Chains: {:?}", structure.get_chain_ids());
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Feature Flags
//!
//! PDBRust uses feature flags to keep the core library lightweight while offering
//! extensive optional functionality:
//!
//! | Feature | Description | Default |
//! |---------|-------------|---------|
//! | `filter` | Structure filtering, extraction, and cleaning | No |
//! | `descriptors` | Structural descriptors (Rg, composition, B-factor, pLDDT) | No |
//! | `quality` | Quality assessment and reports | No |
//! | `summary` | Unified summaries (requires descriptors + quality) | No |
//! | `rcsb` | RCSB PDB search and download | No |
//! | `rcsb-async` | Async/concurrent bulk downloads with rate limiting | No |
//! | `parallel` | Parallel processing with Rayon | No |
//! | `geometry` | RMSD calculation and structure alignment (Kabsch) | No |
//! | `dssp` | DSSP-like secondary structure assignment | No |
//! | `gzip` | Parse gzip-compressed files (.ent.gz, .pdb.gz) | No |
//! | `analysis` | All analysis features combined | No |
//! | `full` | Everything | No |
//!
//! Enable features in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.6", features = ["filter", "descriptors"] }
//! ```
//!
//! ## Core Features
//!
//! ### Parsing
//!
//! - Parse both PDB and mmCIF files with comprehensive error handling
//! - Automatic format detection based on file content
//! - Support for multiple models (NMR ensembles)
//! - Handle alternate conformations (altlocs)
//!
//! ### Structure Data
//!
//! - ATOM/HETATM records with full coordinate data
//! - SEQRES sequence information
//! - CONECT connectivity records
//! - SSBOND disulfide bond definitions
//! - Header, title, and remark metadata
//!
//! ## Optional Features
//!
//! ### Filtering (`filter` feature)
//!
//! ```rust,ignore
//! // Remove ligands and keep only chain A
//! let cleaned = structure
//!     .remove_ligands()
//!     .keep_only_chain("A")
//!     .keep_only_backbone();
//!
//! // Extract CA coordinates
//! let ca_coords = structure.get_ca_coords(None);
//! ```
//!
//! ### Structural Descriptors (`descriptors` feature)
//!
//! ```rust,ignore
//! // Compute structural properties
//! let rg = structure.radius_of_gyration();
//! let composition = structure.aa_composition();
//! let hydrophobic = structure.hydrophobic_ratio();
//! ```
//!
//! ### Quality Assessment (`quality` feature)
//!
//! ```rust,ignore
//! // Get comprehensive quality report
//! let report = structure.quality_report();
//!
//! if report.is_analysis_ready() {
//!     println!("Structure is ready for analysis");
//! }
//! ```
//!
//! ### RCSB Integration (`rcsb` feature)
//!
//! ```rust,ignore
//! use pdbrust::rcsb::{download_structure, rcsb_search, SearchQuery, FileFormat};
//!
//! // Download from RCSB PDB
//! let structure = download_structure("1UBQ", FileFormat::Pdb)?;
//!
//! // Search RCSB
//! let query = SearchQuery::new()
//!     .with_text("kinase")
//!     .with_resolution_max(2.0);
//! let results = rcsb_search(&query, 10)?;
//! ```
//!
//! ### Async RCSB Downloads (`rcsb-async` feature)
//!
//! ```rust,ignore
//! use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};
//!
//! // Download multiple structures concurrently
//! let pdb_ids = vec!["1UBQ", "8HM2", "4INS"];
//! let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;
//!
//! // With rate limiting options
//! let options = AsyncDownloadOptions::default()
//!     .with_max_concurrent(10)
//!     .with_rate_limit_ms(50);
//! let results = download_multiple_async(&pdb_ids, FileFormat::Cif, Some(options)).await;
//! ```
//!
//! ### Geometry (`geometry` feature)
//!
//! ```rust,ignore
//! use pdbrust::geometry::AtomSelection;
//!
//! // Calculate RMSD between structures
//! let rmsd = structure1.rmsd_to(&structure2)?;  // CA atoms by default
//!
//! // Align structures using Kabsch algorithm
//! let (aligned, result) = mobile.align_to(&target)?;
//! println!("RMSD: {:.4} Å ({} atoms)", result.rmsd, result.num_atoms);
//!
//! // Per-residue RMSD for flexibility analysis
//! let per_res = mobile.per_residue_rmsd_to(&target)?;
//! ```
//!
//! ### Secondary Structure (`dssp` feature)
//!
//! ```rust,ignore
//! // DSSP-like secondary structure assignment
//! let ss = structure.assign_secondary_structure();
//! println!("Helix: {:.1}%", ss.helix_fraction * 100.0);
//! println!("Sheet: {:.1}%", ss.sheet_fraction * 100.0);
//!
//! // Compact string representation (e.g., "HHHHEEEECCCC")
//! let ss_string = structure.secondary_structure_string();
//! ```
//!
//! ### Gzip Support (`gzip` feature)
//!
//! ```rust,ignore
//! use pdbrust::{parse_gzip_structure_file, parse_gzip_pdb_file};
//!
//! // Parse gzip-compressed files directly
//! let structure = parse_gzip_structure_file("pdb1ubq.ent.gz")?;
//! let structure = parse_gzip_pdb_file("protein.pdb.gz")?;
//! ```
//!
//! ## Format Support
//!
//! ### PDB Format
//!
//! Traditional fixed-width text format with support for:
//! - ATOM/HETATM records
//! - HEADER, TITLE, REMARK records
//! - SEQRES, CONECT, SSBOND records
//! - MODEL/ENDMDL for multi-model structures
//!
//! ### mmCIF Format
//!
//! Modern dictionary-based format with support for:
//! - `_atom_site` category (converted to ATOM records)
//! - `_entity_poly_seq` category (converted to SEQRES records)
//! - `_struct_disulfid` category (converted to SSBOND records)
//! - Header and metadata information
//!
//! ## Examples
//!
//! ### Auto-detect Format
//!
//! ```rust,no_run
//! use pdbrust::parse_structure_file;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Works with both .pdb and .cif files
//!     let structure = parse_structure_file("example.cif")?;
//!
//!     // Get all chain IDs
//!     let chains = structure.get_chain_ids();
//!
//!     // Get sequence for a specific chain
//!     if let Some(chain_id) = chains.first() {
//!         let sequence = structure.get_sequence(chain_id);
//!         println!("Sequence for chain {}: {:?}", chain_id, sequence);
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Format-Specific Parsing
//!
//! ```rust,no_run
//! use pdbrust::{parse_pdb_file, parse_mmcif_file};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Parse PDB format explicitly
//!     let pdb_structure = parse_pdb_file("structure.pdb")?;
//!
//!     // Parse mmCIF format explicitly
//!     let mmcif_structure = parse_mmcif_file("structure.cif")?;
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Performance
//!
//! PDBRust is designed for high performance. Benchmarks against the Python
//! `libraryPDB` library show 40-260x speedups for common operations:
//!
//! | Operation | Python | Rust | Speedup |
//! |-----------|--------|------|---------|
//! | Parse PDB file | 15ms | 0.36ms | 42x |
//! | Remove ligands | 8ms | 0.03ms | 267x |
//! | Radius of gyration | 2ms | 0.05ms | 40x |
//!
//! ## Error Handling
//!
//! All parsing functions return `Result<T, PdbError>` with detailed error context:
//!
//! ```rust,no_run
//! use pdbrust::{parse_pdb_file, PdbError};
//!
//! match parse_pdb_file("structure.pdb") {
//!     Ok(structure) => println!("Loaded {} atoms", structure.atoms.len()),
//!     Err(PdbError::IoError(e)) => eprintln!("File error: {}", e),
//!     Err(e) => eprintln!("Parse error: {}", e),
//! }
//! ```

// Module exports
pub mod core;
pub mod error;
pub mod guide;
pub mod parser;
pub mod records;
mod utils;
pub mod writer;

// Conditional modules
#[cfg(feature = "filter")]
pub mod filter;

// Re-export selection error type for convenience
#[cfg(feature = "filter")]
pub use filter::selection::SelectionError;

#[cfg(feature = "descriptors")]
pub mod descriptors;

// Re-export descriptor types for convenience
#[cfg(feature = "descriptors")]
pub use descriptors::{ResidueBFactor, StructureDescriptors};

// Re-export AlphaFold/pLDDT types
#[cfg(feature = "descriptors")]
pub use descriptors::{ConfidenceCategory, ResiduePlddt};

// Re-export protein-ligand interaction types
#[cfg(feature = "descriptors")]
pub use descriptors::{
    BindingSite, ContactResidue, HydrophobicContact, LigandInteractionProfile, ProteinLigandHBond,
    SaltBridge,
};

// Re-export dihedral types (requires dssp + descriptors)
#[cfg(all(feature = "descriptors", feature = "dssp"))]
pub use descriptors::{RamachandranRegion, RamachandranStats, ResidueDihedrals, ResidueRef};

// Re-export hydrogen bond types (requires dssp + descriptors)
#[cfg(all(feature = "descriptors", feature = "dssp"))]
pub use descriptors::{HBondStats, HBondType, MainchainHBond, ResidueHBonds};

#[cfg(feature = "quality")]
pub mod quality;

#[cfg(feature = "summary")]
pub mod summary;

#[cfg(feature = "rcsb")]
pub mod rcsb;

#[cfg(feature = "geometry")]
pub mod geometry;

#[cfg(feature = "dssp")]
pub mod dssp;

#[cfg(feature = "ligand-quality")]
pub mod ligand_quality;

// Re-export ligand quality types for convenience
#[cfg(feature = "ligand-quality")]
pub use ligand_quality::{AtomClash, LigandPoseReport};

// Re-exports for convenience
pub use core::PdbStructure;
pub use error::PdbError;
pub use parser::{
    parse_mmcif_file, parse_mmcif_string, parse_pdb_file, parse_pdb_reader, parse_pdb_string,
    parse_structure_file,
};
pub use records::{Atom, Conect, Model, Remark, Residue, SSBond, SeqRes};
pub use writer::{write_mmcif, write_mmcif_file, write_mmcif_string, write_pdb, write_pdb_file};

// Gzip parsing functions (requires "gzip" feature)
#[cfg(feature = "gzip")]
pub use parser::{
    parse_gzip_mmcif_file, parse_gzip_mmcif_reader, parse_gzip_pdb_file, parse_gzip_pdb_reader,
    parse_gzip_structure_file,
};

// Gzip writing functions (requires "gzip" feature)
#[cfg(feature = "gzip")]
pub use writer::write_gzip_mmcif_file;
