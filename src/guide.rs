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
//! | `filtering_demo.rs` | filter | Fluent filtering API and method chaining |
//! | `selection_demo.rs` | filter | PyMOL/VMD-style selection language |
//! | `geometry_demo.rs` | geometry | RMSD calculation and Kabsch alignment |
//! | `lddt_demo.rs` | geometry | LDDT calculation (superposition-free) |
//! | `rcsb_workflow.rs` | rcsb, descriptors | RCSB search and download workflows |
//! | `async_download_demo.rs` | rcsb-async, descriptors | Concurrent bulk downloads with rate limiting |
//! | `batch_processing.rs` | descriptors, summary | Multi-file processing with CSV export |
//! | `b_factor_demo.rs` | descriptors | B-factor analysis and flexibility detection |
//! | `secondary_structure_demo.rs` | dssp | DSSP-like secondary structure assignment |
//! | `dockq_demo.rs` | dockq | DockQ v2 interface quality assessment |
//! | `full_pdb_benchmark.rs` | gzip, parallel, descriptors, quality, summary | Full PDB archive benchmark |
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
//! pdbrust = "0.6"
//! ```
//!
//! For additional features, enable them explicitly:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.6", features = ["filter", "descriptors", "quality", "summary", "rcsb"] }
//! ```
//!
//! ## Feature Flags
//!
//! | Feature | Description |
//! |---------|-------------|
//! | `filter` | Filtering, extraction, cleaning, and selection language |
//! | `descriptors` | Structural descriptors (Rg, composition, B-factor, pLDDT) |
//! | `quality` | Quality assessment and reports |
//! | `summary` | Unified summaries (requires `descriptors` + `quality`) |
//! | `rcsb` | RCSB PDB search and download |
//! | `rcsb-async` | Async/concurrent bulk downloads with rate limiting |
//! | `parallel` | Parallel processing with Rayon |
//! | `geometry` | RMSD, LDDT (superposition-free), structure alignment (Kabsch) |
//! | `dssp` | DSSP-like secondary structure assignment |
//! | `dockq` | DockQ v2 interface quality for protein-protein complexes |
//! | `gzip` | Parse gzip-compressed files (.ent.gz, .pdb.gz) |
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
//! ### Selection Language
//!
//! PDBRust supports PyMOL/VMD-style selection syntax:
//!
//! ```ignore
//! // Basic selectors
//! let chain_a = structure.select("chain A")?;
//! let ca_atoms = structure.select("name CA")?;
//! let alanines = structure.select("resname ALA")?;
//! let residue_range = structure.select("resid 1:100")?;
//!
//! // Keywords
//! let backbone = structure.select("backbone")?;       // N, CA, C, O
//! let protein = structure.select("protein")?;         // Standard amino acids
//! let water = structure.select("water")?;             // HOH, WAT
//! let ligands = structure.select("hetero and not water")?;
//!
//! // Boolean operators
//! let complex = structure.select("chain A and name CA")?;
//! let either = structure.select("chain A or chain B")?;
//! let exclude = structure.select("protein and not hydrogen")?;
//!
//! // Numeric comparisons
//! let low_b = structure.select("bfactor < 30.0")?;
//! let high_occ = structure.select("occupancy >= 0.5")?;
//!
//! // Validate syntax without executing
//! PdbStructure::validate_selection("chain A and name CA")?;
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
//! let aromatic = structure.aromatic_ratio();   // PHE, TYR, TRP
//! let small = structure.small_ratio();         // GLY, ALA, SER, PRO
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
//! ### B-Factor Analysis
//!
//! Analyze crystallographic B-factors (temperature factors):
//!
//! ```ignore
//! // Global statistics
//! let mean_b = structure.b_factor_mean();
//! let mean_ca = structure.b_factor_mean_ca();  // CA atoms only
//! let min_b = structure.b_factor_min();
//! let max_b = structure.b_factor_max();
//! let std_b = structure.b_factor_std();
//!
//! // Per-residue B-factor profile
//! let profile = structure.b_factor_profile();
//! for res in &profile {
//!     println!("{}{}: mean={:.1} min={:.1} max={:.1}",
//!         res.chain_id, res.residue_seq, res.mean, res.min, res.max);
//! }
//!
//! // Identify flexible/rigid regions
//! let flexible = structure.flexible_residues(50.0);  // B > 50 Å²
//! let rigid = structure.rigid_residues(20.0);        // B < 20 Å²
//!
//! // Normalize B-factors (z-score)
//! let normalized = structure.normalize_b_factors();
//! ```
//!
//! ### AlphaFold/pLDDT Analysis
//!
//! Analyze AlphaFold predicted structures and confidence scores:
//!
//! ```ignore
//! // Detect if structure is from AlphaFold
//! if structure.is_predicted() {
//!     // pLDDT scores are stored in B-factor column
//!     let mean_plddt = structure.plddt_mean();
//!     println!("Mean pLDDT: {:.1}", mean_plddt);
//!
//!     // Per-residue pLDDT scores
//!     let scores = structure.per_residue_plddt();
//!     for score in &scores {
//!         println!("{}{}: {:.1} ({})",
//!             score.chain_id, score.residue_seq,
//!             score.plddt, score.category);
//!     }
//!
//!     // Find disordered/confident regions
//!     let low = structure.low_confidence_regions(50.0);   // pLDDT < 50
//!     let high = structure.high_confidence_regions(90.0); // pLDDT > 90
//! }
//! ```
//!
//! pLDDT confidence categories:
//!
//! | Range | Category | Interpretation |
//! |-------|----------|----------------|
//! | ≥90 | Very High | High accuracy, suitable for atomic details |
//! | 70-90 | Confident | Good backbone prediction |
//! | 50-70 | Low | Caution advised |
//! | <50 | Very Low | May be disordered/unstructured |
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
//! ### Async Downloads (feature: `rcsb-async`)
//!
//! Download multiple structures concurrently with rate limiting:
//!
//! ```ignore
//! use pdbrust::rcsb::{download_multiple_async, AsyncDownloadOptions, FileFormat};
//!
//! // Download with default settings
//! let pdb_ids = vec!["1UBQ", "8HM2", "4INS"];
//! let results = download_multiple_async(&pdb_ids, FileFormat::Pdb, None).await;
//!
//! // With custom options
//! let options = AsyncDownloadOptions::default()
//!     .with_max_concurrent(10)   // Max concurrent downloads
//!     .with_rate_limit_ms(50);   // Delay between requests
//!
//! let results = download_multiple_async(&pdb_ids, FileFormat::Cif, Some(options)).await;
//!
//! // Preset configurations
//! let conservative = AsyncDownloadOptions::conservative(); // 2 concurrent, 500ms
//! let fast = AsyncDownloadOptions::fast();                 // 20 concurrent, 25ms
//!
//! // Handle results
//! for (pdb_id, result) in results {
//!     match result {
//!         Ok(structure) => println!("{}: {} atoms", pdb_id, structure.atoms.len()),
//!         Err(e) => eprintln!("{}: {}", pdb_id, e),
//!     }
//! }
//! ```
//!
//! ## Geometry: RMSD, LDDT, and Alignment (feature: `geometry`)
//!
//! Calculate RMSD, LDDT, and align structures using the Kabsch algorithm.
//!
//! ### RMSD Calculation
//!
//! ```ignore
//! use pdbrust::geometry::AtomSelection;
//!
//! // RMSD using CA atoms (default)
//! let rmsd = structure1.rmsd_to(&structure2)?;
//!
//! // RMSD with different atom selections
//! let rmsd_bb = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;
//! let rmsd_all = structure1.rmsd_to_with_selection(&structure2, AtomSelection::AllAtoms)?;
//! ```
//!
//! ### LDDT Calculation (Superposition-Free)
//!
//! LDDT (Local Distance Difference Test) is a superposition-free metric widely used
//! in AlphaFold (pLDDT) and CASP evaluations. It measures the fraction of inter-atomic
//! distances preserved within specified thresholds.
//!
//! ```ignore
//! use pdbrust::geometry::{AtomSelection, LddtOptions};
//!
//! // LDDT using CA atoms (default)
//! let result = model.lddt_to(&reference)?;
//! println!("LDDT: {:.4}", result.score);  // 0.0 (poor) to 1.0 (perfect)
//!
//! // LDDT with custom options
//! let options = LddtOptions::default()
//!     .with_inclusion_radius(10.0)
//!     .with_thresholds(vec![0.5, 1.0, 2.0, 4.0]);
//! let result = model.lddt_to_with_options(&reference, AtomSelection::Backbone, options)?;
//!
//! // Per-residue LDDT for quality analysis
//! let per_res = model.per_residue_lddt_to(&reference)?;
//! for r in per_res.iter().filter(|r| r.score < 0.7) {
//!     println!("Low LDDT: {}{} {:.2}", r.residue_id.0, r.residue_id.1, r.score);
//! }
//! ```
//!
//! Key properties of LDDT:
//! - Superposition-free: Invariant to rotation and translation
//! - Local: Focuses on distance preservation within 15Å (configurable)
//! - Multi-threshold: Uses 0.5Å, 1.0Å, 2.0Å, 4.0Å by default
//! - Range: 0.0 (poor) to 1.0 (perfect)
//!
//! ### Structure Alignment
//!
//! ```ignore
//! // Align mobile to target (returns aligned structure and result)
//! let (aligned, result) = mobile.align_to(&target)?;
//! println!("RMSD after alignment: {:.4} Å", result.rmsd);
//! println!("Atoms used: {}", result.num_atoms);
//!
//! // Align with specific atom selection
//! let (aligned, result) = mobile.align_to_with_selection(&target, AtomSelection::Backbone)?;
//! ```
//!
//! ### Per-Residue RMSD
//!
//! Identify flexible regions by computing per-residue RMSD:
//!
//! ```ignore
//! let per_res = mobile.per_residue_rmsd_to(&target)?;
//! for r in &per_res {
//!     if r.rmsd > 2.0 {
//!         println!("Flexible: {}{} RMSD={:.2} Å", r.residue_id.0, r.residue_id.1, r.rmsd);
//!     }
//! }
//! ```
//!
//! ## Secondary Structure (feature: `dssp`)
//!
//! DSSP-like secondary structure assignment from hydrogen bond patterns.
//!
//! ```ignore
//! // Complete secondary structure assignment
//! let ss = structure.assign_secondary_structure();
//!
//! // Summary statistics
//! println!("Helix: {:.1}%", ss.helix_fraction * 100.0);
//! println!("Sheet: {:.1}%", ss.sheet_fraction * 100.0);
//! println!("Coil:  {:.1}%", ss.coil_fraction * 100.0);
//!
//! // Compact string representation (e.g., "HHHHEEEECCCC")
//! let ss_string = structure.secondary_structure_string();
//!
//! // Composition tuple (helix, sheet, coil)
//! let (helix, sheet, coil) = structure.secondary_structure_composition();
//!
//! // Per-residue assignments
//! for res in &ss.residue_assignments {
//!     println!("{}{}: {} ({})",
//!         res.chain_id, res.residue_seq, res.residue_name, res.ss.code());
//! }
//! ```
//!
//! Secondary structure codes (DSSP):
//!
//! | Code | Name | Description |
//! |------|------|-------------|
//! | H | α-helix | i → i+4 hydrogen bond pattern |
//! | G | 3₁₀-helix | i → i+3 hydrogen bond pattern |
//! | I | π-helix | i → i+5 hydrogen bond pattern |
//! | P | κ-helix | PPII helix (polyproline II) |
//! | E | Extended strand | Part of β-sheet |
//! | B | Beta bridge | Isolated β-bridge |
//! | T | Turn | Hydrogen-bonded turn |
//! | S | Bend | High backbone curvature |
//! | C | Coil | None of the above |
//!
//! ## Dihedral Angles (features: `descriptors` + `dssp`)
//!
//! Compute backbone dihedral angles for Ramachandran analysis.
//!
//! ```ignore
//! // Compute phi/psi angles
//! let dihedrals = structure.compute_dihedrals();
//! for d in &dihedrals {
//!     if let (Some(phi), Some(psi)) = (d.phi, d.psi) {
//!         println!("{}{}: φ={:.1}° ψ={:.1}° ({})",
//!             d.chain_id, d.residue_seq, phi, psi, d.region);
//!     }
//! }
//!
//! // Ramachandran region classification
//! let stats = structure.ramachandran_stats();
//! println!("Favored: {:.1}%", stats.favored_fraction * 100.0);
//! println!("Allowed: {:.1}%", stats.allowed_fraction * 100.0);
//! println!("Outliers: {:.1}%", stats.outlier_fraction * 100.0);
//!
//! // Find outliers
//! let outliers = structure.ramachandran_outliers();
//! for o in &outliers {
//!     println!("Outlier: {}{} ({}) φ={:.1}° ψ={:.1}°",
//!         o.chain_id, o.residue_seq, o.residue_name,
//!         o.phi.unwrap_or(0.0), o.psi.unwrap_or(0.0));
//! }
//!
//! // Detect cis peptide bonds
//! let cis_peptides = structure.cis_peptides();
//! ```
//!
//! ## Hydrogen Bond Network (features: `descriptors` + `dssp`)
//!
//! Detect mainchain hydrogen bonds.
//!
//! ```ignore
//! // Get all mainchain H-bonds
//! let hbonds = structure.mainchain_hbonds();
//! for hb in &hbonds {
//!     println!("H-bond: {}{} -> {}{} ({:.2} Å, {:.1}°)",
//!         hb.donor_chain, hb.donor_seq,
//!         hb.acceptor_chain, hb.acceptor_seq,
//!         hb.distance, hb.angle);
//! }
//!
//! // H-bonds for a specific residue
//! let res_hbonds = structure.hbonds_for_residue("A", 50);
//!
//! // H-bond statistics
//! let stats = structure.hbond_stats();
//! println!("Total H-bonds: {}", stats.total_hbonds);
//! println!("Helical: {}", stats.helical_hbonds);
//! println!("Sheet: {}", stats.sheet_hbonds);
//! ```
//!
//! ## Protein-Ligand Interactions (feature: `descriptors`)
//!
//! Analyze interactions between proteins and bound ligands.
//!
//! ```ignore
//! // Find binding site residues around a ligand
//! let site = structure.binding_site("ATP", 5.0);  // 5 Å cutoff
//! for res in &site.residues {
//!     println!("{}{} ({}) - distance: {:.2} Å",
//!         res.chain_id, res.residue_seq, res.residue_name, res.distance);
//! }
//!
//! // Detailed interaction profile
//! let interactions = structure.ligand_interactions("ATP", 4.0);
//!
//! // Hydrogen bonds
//! for hb in &interactions.hbonds {
//!     println!("H-bond: {}{} {} - {} ({:.2} Å)",
//!         hb.residue_chain, hb.residue_seq, hb.residue_name,
//!         hb.ligand_atom, hb.distance);
//! }
//!
//! // Salt bridges
//! for sb in &interactions.salt_bridges {
//!     println!("Salt bridge: {}{} - {} ({:.2} Å)",
//!         sb.residue_chain, sb.residue_seq, sb.ligand_atom, sb.distance);
//! }
//!
//! // Hydrophobic contacts
//! for hc in &interactions.hydrophobic_contacts {
//!     println!("Hydrophobic: {}{} - {} ({:.2} Å)",
//!         hc.residue_chain, hc.residue_seq, hc.ligand_atom, hc.distance);
//! }
//!
//! // Analyze all ligands at once
//! let all_interactions = structure.all_ligand_interactions(4.0);
//! ```
//!
//! ## Gzip Support (feature: `gzip`)
//!
//! Parse gzip-compressed structure files directly.
//!
//! ```ignore
//! use pdbrust::{parse_gzip_structure_file, parse_gzip_pdb_file, parse_gzip_mmcif_file};
//!
//! // Auto-detect format within gzip
//! let structure = parse_gzip_structure_file("pdb1ubq.ent.gz")?;
//!
//! // Explicit format
//! let structure = parse_gzip_pdb_file("protein.pdb.gz")?;
//! let structure = parse_gzip_mmcif_file("protein.cif.gz")?;
//!
//! // Write gzip-compressed mmCIF
//! use pdbrust::write_gzip_mmcif_file;
//! write_gzip_mmcif_file(&structure, "output.cif.gz")?;
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
