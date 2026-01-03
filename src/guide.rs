//! # PDBRust User Guide
//!
//! This guide provides comprehensive documentation for using the PDBRust library
//! to work with Protein Data Bank (PDB) and mmCIF files.
//!
//! For a quick start guide, see [`docs/GETTING_STARTED.md`](https://github.com/HFooladi/pdbrust/blob/main/docs/GETTING_STARTED.md).
//!
//! ## Runnable Examples
//!
//! The [`examples/`](https://github.com/HFooladi/pdbrust/tree/main/examples) directory contains
//! complete, runnable examples:
//!
//! | Example | Features | Description |
//! |---------|----------|-------------|
//! | `analysis_workflow.rs` | filter, descriptors, quality, summary | Complete load→clean→analyze→export pipeline |
//! | `filtering_demo.rs` | filter | Fluent filtering API demonstration |
//! | `rcsb_workflow.rs` | rcsb, descriptors | RCSB search and download workflows |
//! | `batch_processing.rs` | descriptors, summary | Multi-file processing with CSV export |
//!
//! Run examples with:
//! ```bash
//! cargo run --example analysis_workflow --features "filter,descriptors,quality,summary"
//! cargo run --example filtering_demo --features "filter"
//! ```
//!
//! ## Getting Started
//!
//! Add PDBRust to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! pdbrust = "0.2"
//! ```
//!
//! For additional features, enable them explicitly:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.2", features = ["filter", "descriptors", "quality", "summary", "rcsb"] }
//! ```
//!
//! ## Feature Flags
//!
//! | Feature | Description |
//! |---------|-------------|
//! | `filter` | Filtering, extraction, and structure cleaning |
//! | `descriptors` | Structural descriptors (Rg, composition, geometry) |
//! | `quality` | Quality assessment and reports |
//! | `summary` | Unified summaries (requires `descriptors` + `quality`) |
//! | `rcsb` | RCSB PDB search and download |
//! | `parallel` | Parallel processing with Rayon |
//! | `geometry` | Geometric analysis with nalgebra |
//! | `analysis` | All analysis features combined |
//! | `full` | Everything |
//!
//! ## Basic Usage
//!
//! ### Parsing Structures
//!
//! ```ignore
//! use pdbrust::{parse_pdb_file, parse_mmcif_file, parse_structure_file};
//!
//! // Parse PDB format
//! let structure = parse_pdb_file("protein.pdb")?;
//!
//! // Parse mmCIF format
//! let structure = parse_mmcif_file("protein.cif")?;
//!
//! // Auto-detect format
//! let structure = parse_structure_file("protein.ent")?;
//!
//! // Basic information
//! println!("Atoms: {}", structure.atoms.len());
//! println!("Chains: {:?}", structure.get_chain_ids());
//! println!("Residues: {}", structure.get_residues().len());
//! ```
//!
//! ## Filtering and Cleaning (feature: `filter`)
//!
//! The filter module provides functions to extract subsets of atoms and clean structures.
//!
//! ### Extracting Coordinates
//!
//! ```ignore
//! use pdbrust::parse_pdb_file;
//!
//! let structure = parse_pdb_file("protein.pdb")?;
//!
//! // Get all CA coordinates
//! let ca_coords = structure.get_ca_coords(None);
//!
//! // Get CA coordinates for chain A only
//! let chain_a_ca = structure.get_ca_coords(Some("A"));
//!
//! // Get CA atoms (full Atom structs)
//! let ca_atoms = structure.get_ca_atoms(None);
//! ```
//!
//! ### Filtering Structures
//!
//! ```ignore
//! // Remove ligands and water molecules
//! let protein_only = structure.remove_ligands();
//!
//! // Keep only a specific chain
//! let chain_a = structure.keep_only_chain("A");
//!
//! // Keep only CA atoms
//! let ca_only = structure.keep_only_ca();
//!
//! // Keep only backbone atoms (N, CA, C, O)
//! let backbone = structure.keep_only_backbone();
//!
//! // Remove hydrogen atoms
//! let no_h = structure.remove_hydrogens();
//!
//! // Chain operations using fluent API
//! let cleaned = structure
//!     .remove_ligands()
//!     .keep_only_chain("A")
//!     .keep_only_backbone();
//! ```
//!
//! ### Cleaning and Renumbering
//!
//! ```ignore
//! let mut structure = parse_pdb_file("protein.pdb")?;
//!
//! // Normalize chain IDs to A, B, C, ...
//! structure.normalize_chain_ids();
//!
//! // Reindex residues starting from 1
//! structure.reindex_residues();
//!
//! // Renumber atoms sequentially
//! structure.renumber_atoms();
//!
//! // Center structure at origin
//! structure.center_structure();
//!
//! // Get centroid
//! let (cx, cy, cz) = structure.get_centroid();
//! ```
//!
//! ## Structural Descriptors (feature: `descriptors`)
//!
//! Compute structural properties and metrics.
//!
//! ### Composition Analysis
//!
//! ```ignore
//! use pdbrust::parse_pdb_file;
//!
//! let structure = parse_pdb_file("protein.pdb")?;
//!
//! // Amino acid composition (fractions sum to 1.0)
//! let composition = structure.aa_composition();
//! for (aa, fraction) in &composition {
//!     println!("{}: {:.1}%", aa, fraction * 100.0);
//! }
//!
//! // Specific ratios
//! let gly_ratio = structure.glycine_ratio();
//! let hydrophobic = structure.hydrophobic_ratio();
//! let polar = structure.polar_ratio();
//! let charged = structure.charged_ratio();
//!
//! // Number of residues (based on CA count)
//! let n_residues = structure.count_ca_residues();
//!
//! // Missing residue ratio (sequence gaps)
//! let missing = structure.missing_residue_ratio();
//! ```
//!
//! ### Geometric Descriptors
//!
//! ```ignore
//! // Radius of gyration
//! let rg = structure.radius_of_gyration();
//! println!("Rg: {:.2} Å", rg);
//!
//! // Maximum CA-CA distance
//! let max_dist = structure.max_ca_distance();
//! println!("Max distance: {:.2} Å", max_dist);
//!
//! // Secondary structure content (heuristic)
//! let ss_ratio = structure.secondary_structure_ratio();
//!
//! // Compactness index (Rg / n^(1/3))
//! let compactness = structure.compactness_index();
//!
//! // CA density (atoms per Å³)
//! let density = structure.ca_density();
//!
//! // Bounding box
//! let ((xmin, ymin, zmin), (xmax, ymax, zmax)) = structure.ca_bounding_box();
//! ```
//!
//! ### All Descriptors at Once
//!
//! ```ignore
//! use pdbrust::descriptors::StructureDescriptors;
//!
//! let descriptors = structure.structure_descriptors();
//!
//! println!("Residues: {}", descriptors.num_residues);
//! println!("Rg: {:.2} Å", descriptors.radius_of_gyration);
//! println!("Hydrophobic: {:.1}%", descriptors.hydrophobic_ratio * 100.0);
//! ```
//!
//! ## Quality Assessment (feature: `quality`)
//!
//! Assess structure quality and suitability for analysis.
//!
//! ```ignore
//! use pdbrust::parse_pdb_file;
//!
//! let structure = parse_pdb_file("protein.pdb")?;
//!
//! // Individual checks
//! let is_ca_only = structure.has_ca_only();
//! let has_models = structure.has_multiple_models();
//! let has_altlocs = structure.has_altlocs();
//!
//! // Full quality report
//! let report = structure.quality_report();
//!
//! println!("Chains: {}", report.num_chains);
//! println!("Models: {}", report.num_models);
//! println!("Has HETATM: {}", report.has_hetatm);
//! println!("Has hydrogens: {}", report.has_hydrogens);
//! println!("Has SS bonds: {}", report.has_ssbonds);
//!
//! // Check if suitable for analysis
//! if report.is_analysis_ready() {
//!     println!("Structure is ready for analysis");
//! }
//!
//! if report.is_clean() {
//!     println!("Structure is clean (full atoms, no altlocs)");
//! }
//!
//! // Get resolution if available
//! if let Some(res) = structure.get_resolution() {
//!     println!("Resolution: {:.2} Å", res);
//! }
//! ```
//!
//! ## Unified Summaries (feature: `summary`)
//!
//! Combine quality and descriptors into a single summary.
//!
//! ```ignore
//! use pdbrust::parse_pdb_file;
//! use pdbrust::summary::{batch_summarize, summaries_to_csv};
//!
//! let structure = parse_pdb_file("protein.pdb")?;
//!
//! // Get comprehensive summary
//! let summary = structure.summary();
//!
//! println!("Atoms: {}", summary.num_atoms);
//! println!("Residues: {}", summary.num_residues);
//! println!("Rg: {:.2} Å", summary.radius_of_gyration);
//! println!("Analysis ready: {}", summary.is_analysis_ready());
//!
//! // Batch processing
//! let structures = vec![structure1, structure2, structure3];
//! let summaries = batch_summarize(&structures);
//!
//! // Export to CSV
//! let csv = summaries_to_csv(&summaries, true); // with header
//! std::fs::write("summaries.csv", csv)?;
//! ```
//!
//! ## RCSB PDB Integration (feature: `rcsb`)
//!
//! Search and download structures from RCSB PDB.
//!
//! ### Downloading Structures
//!
//! ```ignore
//! use pdbrust::rcsb::{download_structure, download_to_file, FileFormat};
//!
//! // Download and parse directly
//! let structure = download_structure("1UBQ", FileFormat::Pdb)?;
//! println!("Downloaded {} atoms", structure.atoms.len());
//!
//! // Download mmCIF format
//! let structure = download_structure("1UBQ", FileFormat::Cif)?;
//!
//! // Download to file
//! download_to_file("1UBQ", "1UBQ.pdb", FileFormat::Pdb)?;
//! ```
//!
//! ### Searching RCSB
//!
//! ```ignore
//! use pdbrust::rcsb::{rcsb_search, SearchQuery, ExperimentalMethod};
//!
//! // Build a search query
//! let query = SearchQuery::new()
//!     .with_text("kinase")
//!     .with_organism("Homo sapiens")
//!     .with_resolution_max(2.0)
//!     .with_experimental_method(ExperimentalMethod::XRay);
//!
//! // Execute search (returns up to 10 results)
//! let result = rcsb_search(&query, 10)?;
//!
//! println!("Found {} structures", result.total_count);
//! for pdb_id in &result.pdb_ids {
//!     println!("  {}", pdb_id);
//! }
//! ```
//!
//! ### Search Query Options
//!
//! ```ignore
//! let query = SearchQuery::new()
//!     .with_text("ubiquitin")              // Full-text search
//!     .with_organism("Homo sapiens")       // Source organism
//!     .with_resolution_min(1.0)            // Min resolution (Å)
//!     .with_resolution_max(2.5)            // Max resolution (Å)
//!     .with_experimental_method(ExperimentalMethod::XRay)
//!     .with_polymer_type(PolymerType::Protein)
//!     .with_sequence_length_min(50)        // Min sequence length
//!     .with_sequence_length_max(200)       // Max sequence length
//!     .with_release_date_min("2020-01-01") // Released after
//!     .with_release_date_max("2024-12-31") // Released before
//!     .with_ec_number("2.7.11.1");         // Enzyme classification
//! ```
//!
//! ## Working with Multi-Model Structures
//!
//! ```ignore
//! let structure = parse_pdb_file("nmr_ensemble.pdb")?;
//!
//! // Check for multiple models
//! if structure.has_multiple_models() {
//!     println!("NMR ensemble with {} models", structure.models.len());
//!
//!     for model in &structure.models {
//!         println!("Model {}: {} atoms", model.serial, model.atoms.len());
//!     }
//! }
//! ```
//!
//! ## Connectivity Analysis
//!
//! ```ignore
//! // Disulfide bonds
//! for bond in &structure.ssbonds {
//!     println!("SS-bond: {}{} - {}{}",
//!         bond.chain1_id, bond.residue1_seq,
//!         bond.chain2_id, bond.residue2_seq);
//!     println!("  Distance: {:.2} Å", bond.length);
//! }
//!
//! // Atom connectivity (CONECT records)
//! let connected = structure.get_connected_atoms(atom_serial);
//! ```
//!
//! ## Writing Structures
//!
//! ```ignore
//! use pdbrust::write_pdb_file;
//!
//! // Write to PDB format
//! write_pdb_file(&structure, "output.pdb")?;
//!
//! // Or use the structure method
//! structure.to_file("output.pdb")?;
//! ```
//!
//! ## Error Handling
//!
//! ```ignore
//! use pdbrust::{parse_pdb_file, PdbError};
//!
//! match parse_pdb_file("structure.pdb") {
//!     Ok(structure) => {
//!         println!("Loaded {} atoms", structure.atoms.len());
//!     }
//!     Err(PdbError::Io(e)) => {
//!         eprintln!("File error: {}", e);
//!     }
//!     Err(PdbError::InvalidRecord(msg)) => {
//!         eprintln!("Invalid record: {}", msg);
//!     }
//!     Err(e) => {
//!         eprintln!("Error: {}", e);
//!     }
//! }
//! ```
//!
//! ## Common Workflows
//!
//! ### Workflow 1: Structure Cleaning for MD Simulations
//!
//! ```ignore
//! use pdbrust::{parse_pdb_file, write_pdb_file};
//!
//! let structure = parse_pdb_file("raw_structure.pdb")?;
//!
//! // Clean: remove ligands, waters, hydrogens
//! let mut cleaned = structure
//!     .remove_ligands()
//!     .remove_hydrogens()
//!     .keep_only_chain("A");
//!
//! // Prepare for simulation
//! cleaned.center_structure();
//! cleaned.renumber_atoms();
//! cleaned.reindex_residues();
//!
//! write_pdb_file(&cleaned, "cleaned_for_md.pdb")?;
//! ```
//!
//! ### Workflow 2: Dataset Characterization
//!
//! ```ignore
//! use pdbrust::{parse_structure_file, PdbStructure};
//! use pdbrust::summary::{batch_summarize, summaries_to_csv};
//! use std::fs;
//!
//! // Load multiple structures
//! let mut structures: Vec<PdbStructure> = Vec::new();
//! for entry in fs::read_dir("pdb_files")? {
//!     let path = entry?.path();
//!     if let Ok(s) = parse_structure_file(&path) {
//!         structures.push(s);
//!     }
//! }
//!
//! // Compute summaries and export
//! let summaries = batch_summarize(&structures);
//! let csv = summaries_to_csv(&summaries, true);
//! fs::write("dataset_summary.csv", csv)?;
//! ```
//!
//! ### Workflow 3: RCSB Search and Analysis Pipeline
//!
//! ```ignore
//! use pdbrust::rcsb::{rcsb_search, download_structure, SearchQuery, FileFormat};
//!
//! // Search for structures
//! let query = SearchQuery::new()
//!     .with_text("kinase")
//!     .with_organism("Homo sapiens")
//!     .with_resolution_max(2.0);
//!
//! let results = rcsb_search(&query, 10)?;
//!
//! // Download and analyze each
//! for pdb_id in &results.pdb_ids {
//!     let structure = download_structure(pdb_id, FileFormat::Pdb)?;
//!     let rg = structure.radius_of_gyration();
//!     let composition = structure.aa_composition();
//!     println!("{}: Rg={:.1}Å, {} residues", pdb_id, rg, structure.count_ca_residues());
//! }
//! ```
//!
//! ### Workflow 4: Quality Filtering
//!
//! ```ignore
//! use pdbrust::parse_pdb_file;
//!
//! let structure = parse_pdb_file("structure.pdb")?;
//!
//! // Check quality before analysis
//! let report = structure.quality_report();
//!
//! if !report.is_analysis_ready() {
//!     println!("Warning: Structure may need preprocessing");
//!     if report.has_multiple_models {
//!         println!("  - Multiple models (NMR ensemble)");
//!     }
//!     if report.has_altlocs {
//!         println!("  - Alternate conformations present");
//!     }
//!     if report.has_ca_only {
//!         println!("  - CA-only structure");
//!     }
//! }
//!
//! // Get unified summary for clean structures
//! if report.is_clean() {
//!     let summary = structure.summary();
//!     println!("Rg: {:.2}Å, Hydrophobic: {:.1}%",
//!         summary.radius_of_gyration,
//!         summary.hydrophobic_ratio * 100.0);
//! }
//! ```
//!
//! ## Performance Tips
//!
//! 1. **Parse once, reuse**: Parse the structure once and perform multiple analyses
//! 2. **Use iterators**: Prefer `.iter()` over collecting to vectors
//! 3. **Filter early**: Apply filters before expensive computations
//! 4. **Batch operations**: Use `batch_summarize` for multiple structures
//!
//! ```ignore
//! // Efficient: filter then compute
//! let chain_a = structure.keep_only_chain("A");
//! let rg = chain_a.radius_of_gyration();
//!
//! // Efficient: use iterators
//! let ca_count = structure.atoms.iter()
//!     .filter(|a| a.name.trim() == "CA")
//!     .count();
//! ```
//!
//! ## See Also
//!
//! - [Getting Started Guide](https://github.com/HFooladi/pdbrust/blob/main/docs/GETTING_STARTED.md)
//! - [Examples Directory](https://github.com/HFooladi/pdbrust/tree/main/examples)
//! - [API Documentation](https://docs.rs/pdbrust)
