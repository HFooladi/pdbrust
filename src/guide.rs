//! # PDBRust User Guide
//!
//! This guide provides comprehensive documentation for using the PDBRust library
//! to work with Protein Data Bank (PDB) files.
//!
//! ## Getting Started
//!
//! First, add PDBRust to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! pdbrust = "0.1.0"
//! ```
//!
//! Then, import and use the library:
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Load a PDB file
//!     let structure = PdbStructure::from_file("example.pdb")?;
//!     
//!     // Get basic information
//!     println!("Number of atoms: {}", structure.atoms.len());
//!     println!("Chains: {:?}", structure.get_chain_ids());
//!     
//!     Ok(())
//! }
//! ```
//!
//! ## Common Use Cases
//!
//! ### 1. Analyzing Protein Structure
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! fn analyze_structure(pdb_file: &str) -> Result<(), Box<dyn std::error::Error>> {
//!     let structure = PdbStructure::from_file(pdb_file)?;
//!     
//!     // Get all chains
//!     let chains = structure.get_chain_ids();
//!     
//!     for chain_id in chains {
//!         // Get residues in this chain
//!         let residues = structure.get_residues_for_chain(&chain_id);
//!         println!("Chain {} has {} residues", chain_id, residues.len());
//!         
//!         // Get sequence from SEQRES records
//!         let sequence = structure.get_sequence(&chain_id);
//!         println!("Chain {} sequence length: {}", chain_id, sequence.len());
//!     }
//!     
//!     Ok(())
//! }
//! ```
//!
//! ### 2. Working with Multiple Models (e.g., NMR Structures)
//!
//! ```ignore
//! fn analyze_models(structure: &PdbStructure) {
//!     // Check number of models
//!     println!("Structure contains {} models", structure.models.len());
//!     
//!     // Analyze each model
//!     for model in &structure.models {
//!         println!("Model {} has {} atoms", model.serial, model.atoms.len());
//!         
//!         // Process model-specific remarks
//!         for remark in &model.remarks {
//!             println!("Model {} remark {}: {}",
//!                     model.serial, remark.number, remark.content);
//!         }
//!     }
//! }
//! ```
//!
//! ### 3. Analyzing Connectivity
//!
//! ```ignore
//! fn analyze_connectivity(structure: &PdbStructure) {
//!     // Analyze disulfide bonds
//!     for bond in &structure.ssbonds {
//!         println!("SS-bond between:");
//!         println!("  Residue: {} {} {}",
//!                 bond.residue1_name, bond.chain1_id, bond.residue1_seq);
//!         println!("  Residue: {} {} {}",
//!                 bond.residue2_name, bond.chain2_id, bond.residue2_seq);
//!         if let Some(distance) = bond.distance {
//!             println!("  Distance: {:.2} Å", distance);
//!         }
//!     }
//!     
//!     // Analyze atom connectivity
//!     for atom in &structure.atoms {
//!         let connected = structure.get_connected_atoms(atom.serial);
//!         if !connected.is_empty() {
//!             println!("Atom {} ({}) is connected to:",
//!                     atom.serial, atom.name);
//!             for connected_atom in connected {
//!                 println!("  Atom {} ({})",
//!                         connected_atom.serial, connected_atom.name);
//!             }
//!         }
//!     }
//! }
//! ```
//!
//! ## Best Practices
//!
//! 1. **Error Handling**
//!    - Always use the `Result` type to handle potential errors
//!    - Check for specific error types using pattern matching
//!    ```ignore
//!    match PdbStructure::from_file("structure.pdb") {
//!        Ok(structure) => {
//!            // Process the structure
//!        }
//!        Err(PdbError::Io(e)) => {
//!            eprintln!("File error: {}", e);
//!        }
//!        Err(PdbError::InvalidRecord(msg)) => {
//!            eprintln!("Invalid record: {}", msg);
//!        }
//!        Err(e) => {
//!            eprintln!("Other error: {}", e);
//!        }
//!    }
//!    ```
//!
//! 2. **Memory Management**
//!    - Use references when passing `PdbStructure` to functions
//!    - Consider using iterators instead of collecting vectors for large structures
//!    ```ignore
//!    // Efficient way to process atoms
//!    structure.atoms
//!        .iter()
//!        .filter(|atom| atom.chain_id == "A")
//!        .for_each(|atom| {
//!            // Process each atom
//!        });
//!    ```
//!
//! 3. **Data Validation**
//!    - Validate PDB files before processing
//!    - Check for required records
//!    ```ignore
//!    if structure.atoms.is_empty() {
//!        eprintln!("Warning: No atom records found");
//!    }
//!    if structure.seqres.is_empty() {
//!        eprintln!("Warning: No sequence information available");
//!    }
//!    ```
//!
//! ## Performance Considerations
//!
//! 1. **Memory Usage**
//!    - The library stores all atoms in memory
//!    - For very large structures, consider processing atoms in chunks
//!    - Use iterators instead of collecting vectors when possible
//!
//! 2. **Processing Speed**
//!    - Parsing is done in a single pass through the file
//!    - Atom lookups are O(n) - consider building indices for frequent lookups
//!    - Chain and residue operations are optimized with sorting and deduplication
//!
//! 3. **Optimization Tips**
//!    ```ignore
//!    // Create indices for frequent lookups
//!    use std::collections::HashMap;
//!    
//!    let atom_map: HashMap<i32, &Atom> = structure.atoms
//!        .iter()
//!        .map(|atom| (atom.serial, atom))
//!        .collect();
//!    
//!    // Now lookups are O(1)
//!    if let Some(atom) = atom_map.get(&serial_number) {
//!        // Process atom
//!    }
//!    ```
//!
//! ## Examples of Handling Different PDB Record Types
//!
//! ### 1. ATOM Records
//! ```ignore
//! // Processing ATOM records
//! for atom in &structure.atoms {
//!     if atom.name == "CA" {  // Alpha carbons
//!         println!("CA atom at position ({}, {}, {})",
//!                 atom.x, atom.y, atom.z);
//!     }
//! }
//! ```
//!
//! ### 2. SEQRES Records
//! ```ignore
//! // Processing sequence information
//! for seqres in &structure.seqres {
//!     println!("Chain {} sequence:", seqres.chain_id);
//!     for residue in &seqres.residues {
//!         print!("{} ", residue);
//!     }
//!     println!();
//! }
//! ```
//!
//! ### 3. SSBOND Records
//! ```ignore
//! // Processing disulfide bonds
//! for bond in &structure.ssbonds {
//!     println!("Disulfide bond between:");
//!     println!("{} {} {} - {} {} {}",
//!             bond.residue1_name, bond.chain1_id, bond.residue1_seq,
//!             bond.residue2_name, bond.chain2_id, bond.residue2_seq);
//! }
//! ```
//!
//! ### 4. REMARK Records
//! ```ignore
//! // Processing specific remarks
//! let resolution_remarks = structure.get_remarks_by_number(2);
//! for remark in resolution_remarks {
//!     if remark.content.contains("RESOLUTION") {
//!         println!("Structure resolution: {}", remark.content);
//!     }
//! }
//! ```
//!
//! ### 5. MODEL Records
//! ```ignore
//! // Processing multiple models
//! for model in &structure.models {
//!     println!("Model {}:", model.serial);
//!     println!("  Number of atoms: {}", model.atoms.len());
//!     println!("  Number of remarks: {}", model.remarks.len());
//! }
//! ```
