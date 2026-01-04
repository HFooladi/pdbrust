//! Filtering and Cleaning Demo
//!
//! This example demonstrates the fluent filtering API for
//! structure manipulation and cleaning. Common operations include:
//! - Extracting specific chains
//! - Keeping only CA or backbone atoms
//! - Removing ligands, waters, and hydrogens
//! - Centering and renumbering structures
//!
//! Run with:
//! ```bash
//! cargo run --example filtering_demo --features "filter"
//! ```

use pdbrust::parse_pdb_file;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== PDBRust Filtering Demo ===\n");

    // Load structure
    let structure = parse_pdb_file("examples/pdb_files/1UBQ.pdb")?;
    println!("Original structure: {} atoms", structure.atoms.len());
    println!("Chains: {:?}", structure.get_chain_ids());

    // ========== Chain Extraction ==========
    println!("\n--- Chain Extraction ---");

    let chain_a = structure.keep_only_chain("A");
    println!("Chain A only: {} atoms", chain_a.atoms.len());

    // ========== Atom Type Filtering ==========
    println!("\n--- Atom Type Filtering ---");

    // CA atoms only (useful for distance matrices, contact maps)
    let ca_only = structure.keep_only_ca();
    println!("CA atoms only: {} atoms", ca_only.atoms.len());

    // Backbone atoms (N, CA, C, O) - useful for secondary structure
    let backbone = structure.keep_only_backbone();
    println!(
        "Backbone only (N, CA, C, O): {} atoms",
        backbone.atoms.len()
    );

    // Extract CA coordinates as vectors
    let ca_coords = structure.get_ca_coords(None);
    println!("CA coordinates extracted: {} points", ca_coords.len());

    // CA coords for specific chain
    let chain_a_ca_coords = structure.get_ca_coords(Some("A"));
    println!("Chain A CA coordinates: {} points", chain_a_ca_coords.len());

    // ========== Remove Unwanted Atoms ==========
    println!("\n--- Removing Unwanted Atoms ---");

    // Remove ligands and water molecules (HETATM records)
    let no_ligands = structure.remove_ligands();
    println!(
        "After remove_ligands(): {} atoms (removed {})",
        no_ligands.atoms.len(),
        structure.atoms.len() - no_ligands.atoms.len()
    );

    // Remove hydrogen atoms (for lighter representation)
    let no_hydrogens = structure.remove_hydrogens();
    println!(
        "After remove_hydrogens(): {} atoms (removed {})",
        no_hydrogens.atoms.len(),
        structure.atoms.len() - no_hydrogens.atoms.len()
    );

    // ========== Method Chaining (Fluent API) ==========
    println!("\n--- Fluent API (Method Chaining) ---");

    // Chain multiple operations for clean, readable code
    let cleaned = structure
        .remove_ligands()
        .remove_hydrogens()
        .keep_only_chain("A")
        .keep_only_backbone();

    println!("Chained: remove_ligands -> remove_hydrogens -> chain A -> backbone");
    println!("Result: {} atoms", cleaned.atoms.len());

    // Another common workflow: CA-only representation
    let ca_repr = structure.remove_ligands().keep_only_ca();
    println!(
        "\nCA representation (no ligands): {} atoms",
        ca_repr.atoms.len()
    );

    // ========== In-Place Modifications ==========
    println!("\n--- In-Place Modifications ---");

    // Clone for modification (original is immutable)
    let mut mutable = structure.clone();

    // Get original centroid
    let (cx, cy, cz) = mutable.get_centroid();
    println!("Original centroid: ({:.2}, {:.2}, {:.2})", cx, cy, cz);

    // Center structure at origin
    mutable.center_structure();
    let (cx, cy, cz) = mutable.get_centroid();
    println!("After centering: ({:.2}, {:.2}, {:.2})", cx, cy, cz);

    // Renumber atoms sequentially from 1
    mutable.renumber_atoms();
    println!(
        "Atoms renumbered: first={}, last={}",
        mutable.atoms.first().map(|a| a.serial).unwrap_or(0),
        mutable.atoms.last().map(|a| a.serial).unwrap_or(0)
    );

    // Normalize chain IDs to A, B, C, ...
    let mut multi_chain = structure.clone();
    multi_chain.normalize_chain_ids();
    println!("Chains normalized: {:?}", multi_chain.get_chain_ids());

    // Reindex residues starting from 1
    let mut reindexed = structure.clone();
    reindexed.reindex_residues();
    if let Some(first_atom) = reindexed.atoms.first() {
        println!("First residue after reindex: {}", first_atom.residue_seq);
    }

    // ========== Working with Residues ==========
    println!("\n--- Working with Residues ---");

    let residues = structure.get_residues();
    println!("Total residues: {}", residues.len());

    let chain_a_residues = structure.get_residues_for_chain("A");
    println!("Chain A residues: {}", chain_a_residues.len());

    // Get sequence (returns Vec<String> of residue names)
    let sequence = structure.get_sequence("A");
    if !sequence.is_empty() {
        let seq_preview: Vec<_> = sequence.iter().take(10).collect();
        println!(
            "Chain A sequence (first 10 residues): {:?}{}",
            seq_preview,
            if sequence.len() > 10 { "..." } else { "" }
        );
    }

    // ========== Summary ==========
    println!("\n=== Summary ===");
    println!("Filtering operations are non-destructive (return new structures)");
    println!("Use method chaining for clean, readable pipelines");
    println!("In-place modifications require &mut and modify the structure");

    Ok(())
}
