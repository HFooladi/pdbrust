//! Demonstration of the selection language feature.
//!
//! Run with: cargo run --example selection_demo --features filter

use pdbrust::parse_pdb_file;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join("1UBQ.pdb");

    println!("Loading structure from {:?}", pdb_path);
    let structure = parse_pdb_file(&pdb_path)?;
    println!("Total atoms: {}", structure.get_num_atoms());
    println!();

    // ====================
    // Basic Selectors
    // ====================
    println!("=== Basic Selectors ===");

    // Chain selection
    let chain_a = structure.select("chain A")?;
    println!("chain A: {} atoms", chain_a.get_num_atoms());

    // Atom name selection
    let ca_atoms = structure.select("name CA")?;
    println!("name CA: {} atoms", ca_atoms.get_num_atoms());

    // Residue name selection
    let met_residues = structure.select("resname MET")?;
    println!("resname MET: {} atoms", met_residues.get_num_atoms());

    // Single residue selection
    let res1 = structure.select("resid 1")?;
    println!("resid 1: {} atoms", res1.get_num_atoms());

    // Residue range selection
    let res_range = structure.select("resid 1:10")?;
    println!("resid 1:10: {} atoms", res_range.get_num_atoms());

    // Element selection
    let nitrogens = structure.select("element N")?;
    println!("element N: {} atoms", nitrogens.get_num_atoms());
    println!();

    // ====================
    // Keyword Selectors
    // ====================
    println!("=== Keyword Selectors ===");

    let backbone = structure.select("backbone")?;
    println!("backbone: {} atoms", backbone.get_num_atoms());

    let protein = structure.select("protein")?;
    println!("protein: {} atoms", protein.get_num_atoms());

    let all_atoms = structure.select("all")?;
    println!("all: {} atoms", all_atoms.get_num_atoms());
    println!();

    // ====================
    // Boolean Operators
    // ====================
    println!("=== Boolean Operators ===");

    // AND operator
    let chain_a_ca = structure.select("chain A and name CA")?;
    println!("chain A and name CA: {} atoms", chain_a_ca.get_num_atoms());

    // OR operator
    let met_or_gln = structure.select("resname MET or resname GLN")?;
    println!("resname MET or resname GLN: {} atoms", met_or_gln.get_num_atoms());

    // NOT operator
    let not_hydrogen = structure.select("not hydrogen")?;
    println!("not hydrogen: {} atoms", not_hydrogen.get_num_atoms());

    // Parentheses for grouping
    let complex1 = structure.select("(resid 1 or resid 2) and name CA")?;
    println!("(resid 1 or resid 2) and name CA: {} atoms", complex1.get_num_atoms());
    println!();

    // ====================
    // Numeric Comparisons
    // ====================
    println!("=== Numeric Comparisons ===");

    let low_bfactor = structure.select("bfactor < 20.0")?;
    println!("bfactor < 20.0: {} atoms", low_bfactor.get_num_atoms());

    let high_occupancy = structure.select("occupancy >= 0.5")?;
    println!("occupancy >= 0.5: {} atoms", high_occupancy.get_num_atoms());
    println!();

    // ====================
    // Complex Selections
    // ====================
    println!("=== Complex Selections ===");

    // Chain A backbone with low B-factor
    let complex2 = structure.select("chain A and backbone and bfactor < 30.0")?;
    println!(
        "chain A and backbone and bfactor < 30.0: {} atoms",
        complex2.get_num_atoms()
    );

    // First and last 10 residues, CA only
    let termini = structure.select("(resid 1:10 or resid 67:76) and name CA")?;
    println!(
        "(resid 1:10 or resid 67:76) and name CA: {} atoms",
        termini.get_num_atoms()
    );

    // Protein without hydrogen
    let heavy_atoms = structure.select("protein and not hydrogen")?;
    println!("protein and not hydrogen: {} atoms", heavy_atoms.get_num_atoms());
    println!();

    // ====================
    // Equivalence with Existing Methods
    // ====================
    println!("=== Equivalence with Existing Methods ===");

    // select("name CA") == keep_only_ca()
    let ca_via_select = structure.select("name CA")?.get_num_atoms();
    let ca_via_method = structure.keep_only_ca().get_num_atoms();
    println!(
        "select(\"name CA\"): {}, keep_only_ca(): {} - Match: {}",
        ca_via_select,
        ca_via_method,
        ca_via_select == ca_via_method
    );

    // select("backbone") == keep_only_backbone()
    let bb_via_select = structure.select("backbone")?.get_num_atoms();
    let bb_via_method = structure.keep_only_backbone().get_num_atoms();
    println!(
        "select(\"backbone\"): {}, keep_only_backbone(): {} - Match: {}",
        bb_via_select,
        bb_via_method,
        bb_via_select == bb_via_method
    );

    // select("chain A") == keep_only_chain("A")
    let chain_via_select = structure.select("chain A")?.get_num_atoms();
    let chain_via_method = structure.keep_only_chain("A").get_num_atoms();
    println!(
        "select(\"chain A\"): {}, keep_only_chain(\"A\"): {} - Match: {}",
        chain_via_select,
        chain_via_method,
        chain_via_select == chain_via_method
    );
    println!();

    // ====================
    // Selection Validation
    // ====================
    println!("=== Selection Validation ===");

    // Valid selection
    match pdbrust::PdbStructure::validate_selection("chain A and name CA") {
        Ok(()) => println!("\"chain A and name CA\" - Valid"),
        Err(e) => println!("\"chain A and name CA\" - Error: {}", e),
    }

    // Invalid selection (empty)
    match pdbrust::PdbStructure::validate_selection("") {
        Ok(()) => println!("\"\" - Valid"),
        Err(e) => println!("\"\" - Error: {}", e),
    }

    // Invalid selection (unclosed parenthesis)
    match pdbrust::PdbStructure::validate_selection("(chain A and name CA") {
        Ok(()) => println!("\"(chain A and name CA\" - Valid"),
        Err(e) => println!("\"(chain A and name CA\" - Error: {}", e),
    }

    // Invalid selection (unknown keyword)
    match pdbrust::PdbStructure::validate_selection("invalid_keyword") {
        Ok(()) => println!("\"invalid_keyword\" - Valid"),
        Err(e) => println!("\"invalid_keyword\" - Error: {}", e),
    }

    println!("\nDone!");
    Ok(())
}
