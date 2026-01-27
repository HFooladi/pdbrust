//! Validation Example for New Test Files
//!
//! This example validates the functionality of PDBRust against the newly
//! downloaded test files:
//! - AF-P62987-F1.pdb: AlphaFold structure for pLDDT testing
//! - 1HSG.pdb: HIV-1 protease with indinavir (MK1) for ligand interaction testing
//! - 2GB1.cif: Protein G B1 domain for mmCIF parsing
//! - 1L2Y.pdb: Trp-cage miniprotein for multi-model NMR testing
//!
//! Run with:
//! ```bash
//! cargo run --example validate_new_files --features "full"
//! ```

use pdbrust::{parse_mmcif_file, parse_pdb_file, parse_structure_file};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== PDBRust Test File Validation ===\n");

    let mut all_passed = true;

    // Test 1: AlphaFold pLDDT functionality
    all_passed &= validate_alphafold()?;

    // Test 2: Protein-ligand interactions
    all_passed &= validate_ligand_interactions()?;

    // Test 3: mmCIF parsing
    all_passed &= validate_mmcif()?;

    // Test 4: Multi-model NMR structure
    all_passed &= validate_multimodel()?;

    // Summary
    println!("\n=== Validation Summary ===");
    if all_passed {
        println!("All tests PASSED!");
    } else {
        println!("Some tests FAILED. Please check the output above.");
    }

    Ok(())
}

fn validate_alphafold() -> Result<bool, Box<dyn Error>> {
    println!("--- Test 1: AlphaFold Structure (AF-P62987-F1.pdb) ---");

    let structure = parse_pdb_file("examples/pdb_files/AF-P62987-F1.pdb")?;
    let mut passed = true;

    // Check basic parsing
    println!("  Atoms: {}", structure.atoms.len());
    if structure.atoms.is_empty() {
        println!("  [FAIL] No atoms parsed");
        passed = false;
    } else {
        println!("  [PASS] Structure parsed successfully");
    }

    // Check if detected as predicted structure
    let is_predicted = structure.is_predicted();
    println!("  is_predicted(): {}", is_predicted);
    if !is_predicted {
        println!("  [FAIL] Structure not detected as predicted");
        passed = false;
    } else {
        println!("  [PASS] Detected as predicted structure");
    }

    // Check pLDDT mean
    let plddt_mean = structure.plddt_mean();
    println!("  plddt_mean(): {:.2}", plddt_mean);
    if !(30.0..=100.0).contains(&plddt_mean) {
        println!("  [FAIL] pLDDT mean out of expected range");
        passed = false;
    } else {
        println!("  [PASS] pLDDT mean in expected range");
    }

    // Check per-residue pLDDT
    let per_residue = structure.per_residue_plddt();
    println!("  per_residue_plddt(): {} residues", per_residue.len());
    if per_residue.is_empty() {
        println!("  [FAIL] No per-residue pLDDT data");
        passed = false;
    } else {
        println!("  [PASS] Per-residue pLDDT computed");
    }

    // Check confidence regions
    let low_conf = structure.low_confidence_regions(70.0);
    let high_conf = structure.high_confidence_regions(90.0);
    println!("  low_confidence_regions(<70): {} residues", low_conf.len());
    println!(
        "  high_confidence_regions(>=90): {} residues",
        high_conf.len()
    );
    println!("  [PASS] Confidence regions computed");

    // Check pLDDT distribution
    let (very_high, confident, low, very_low) = structure.plddt_distribution();
    println!("  pLDDT distribution:");
    println!("    Very high (>90): {:.1}%", very_high * 100.0);
    println!("    Confident (70-90): {:.1}%", confident * 100.0);
    println!("    Low (50-70): {:.1}%", low * 100.0);
    println!("    Very low (<50): {:.1}%", very_low * 100.0);
    println!("  [PASS] pLDDT distribution computed");

    println!(
        "\n  AlphaFold validation: {}",
        if passed { "PASSED" } else { "FAILED" }
    );
    Ok(passed)
}

fn validate_ligand_interactions() -> Result<bool, Box<dyn Error>> {
    println!("\n--- Test 2: Protein-Ligand Complex (1HSG.pdb) ---");

    let structure = parse_pdb_file("examples/pdb_files/1HSG.pdb")?;
    let mut passed = true;

    // Check basic parsing
    println!("  Atoms: {}", structure.atoms.len());
    println!("  Chains: {:?}", structure.get_chain_ids());

    // Check if ligand MK1 (indinavir) is present
    let ligand_atoms: Vec<_> = structure
        .atoms
        .iter()
        .filter(|a| a.residue_name == "MK1")
        .collect();
    println!("  Ligand MK1 atoms: {}", ligand_atoms.len());
    if ligand_atoms.is_empty() {
        println!("  [FAIL] Ligand MK1 not found");
        passed = false;
    } else {
        println!("  [PASS] Ligand MK1 found");
    }

    // Check binding site
    if let Some(site) = structure.binding_site("MK1", 5.0) {
        println!("  Binding site residues: {}", site.num_residues());
        if site.num_residues() == 0 {
            println!("  [FAIL] No binding site residues found");
            passed = false;
        } else {
            println!("  [PASS] Binding site detected");
            // Show closest residues
            println!("  Closest residues:");
            for res in site.residues_by_distance().iter().take(5) {
                println!(
                    "    {}{} {}: {:.2} Å",
                    res.chain_id, res.residue_seq, res.residue_name, res.min_distance
                );
            }
        }
    } else {
        println!("  [FAIL] binding_site() returned None");
        passed = false;
    }

    // Check ligand interactions
    if let Some(profile) = structure.ligand_interactions("MK1") {
        println!("  Ligand interactions:");
        println!("    H-bonds: {}", profile.hydrogen_bonds.len());
        println!("    Salt bridges: {}", profile.salt_bridges.len());
        println!(
            "    Hydrophobic contacts: {}",
            profile.hydrophobic_contacts.len()
        );
        println!("    Total interactions: {}", profile.total_interactions());

        if profile.total_interactions() == 0 {
            println!("  [WARN] No interactions detected (may need parameter tuning)");
        } else {
            println!("  [PASS] Ligand interactions detected");
        }

        // Show some H-bond details
        if !profile.hydrogen_bonds.is_empty() {
            println!("  Sample H-bonds:");
            for hb in profile.hydrogen_bonds.iter().take(3) {
                println!(
                    "    {}{} {}:{} - {}:{} ({:.2} Å)",
                    hb.protein_chain,
                    hb.protein_resid,
                    hb.protein_resname,
                    hb.protein_atom,
                    hb.ligand_name,
                    hb.ligand_atom,
                    hb.distance
                );
            }
        }
    } else {
        println!("  [FAIL] ligand_interactions() returned None");
        passed = false;
    }

    // Check all_ligand_interactions
    let all_lig = structure.all_ligand_interactions();
    println!("  all_ligand_interactions(): {} ligands", all_lig.len());
    for profile in &all_lig {
        println!(
            "    {} (chain {}): {} contacts",
            profile.ligand_name,
            profile.ligand_chain,
            profile.contact_residues.len()
        );
    }

    println!(
        "\n  Protein-ligand validation: {}",
        if passed { "PASSED" } else { "FAILED" }
    );
    Ok(passed)
}

fn validate_mmcif() -> Result<bool, Box<dyn Error>> {
    println!("\n--- Test 3: mmCIF Parsing (2GB1.cif) ---");

    let mut passed = true;

    // Test explicit mmCIF parsing
    let structure = parse_mmcif_file("examples/pdb_files/2GB1.cif")?;
    println!("  parse_mmcif_file():");
    println!("    Atoms: {}", structure.atoms.len());
    println!("    Chains: {:?}", structure.get_chain_ids());

    if structure.atoms.is_empty() {
        println!("  [FAIL] No atoms parsed from mmCIF");
        passed = false;
    } else {
        println!("  [PASS] mmCIF parsed successfully");
    }

    // Test auto-detection
    let structure2 = parse_structure_file("examples/pdb_files/2GB1.cif")?;
    println!("  parse_structure_file() (auto-detect):");
    println!("    Atoms: {}", structure2.atoms.len());

    if structure2.atoms.len() != structure.atoms.len() {
        println!("  [FAIL] Auto-detection gave different result");
        passed = false;
    } else {
        println!("  [PASS] Auto-detection works correctly");
    }

    // Check chain IDs
    let chains = structure.get_chain_ids();
    if chains.is_empty() {
        println!("  [FAIL] No chains extracted");
        passed = false;
    } else {
        println!("  [PASS] Chains extracted: {:?}", chains);
    }

    // Check residue count
    let residues = structure.get_residues();
    println!("  Residues: {}", residues.len());

    println!(
        "\n  mmCIF validation: {}",
        if passed { "PASSED" } else { "FAILED" }
    );
    Ok(passed)
}

fn validate_multimodel() -> Result<bool, Box<dyn Error>> {
    println!("\n--- Test 4: Multi-model NMR Structure (1L2Y.pdb) ---");

    let structure = parse_pdb_file("examples/pdb_files/1L2Y.pdb")?;
    let mut passed = true;

    // Check basic parsing
    println!("  Total atoms: {}", structure.atoms.len());

    // Check models
    let num_models = structure.models.len();
    println!("  Number of models: {}", num_models);
    if num_models <= 1 {
        println!("  [FAIL] Expected multiple models");
        passed = false;
    } else {
        println!("  [PASS] Multiple models detected");
    }

    // Check has_multiple_models
    let has_multiple = structure.has_multiple_models();
    println!("  has_multiple_models(): {}", has_multiple);
    if !has_multiple && num_models > 1 {
        println!("  [FAIL] has_multiple_models() disagrees with models count");
        passed = false;
    } else {
        println!("  [PASS] has_multiple_models() correct");
    }

    // Show model info
    println!("  Model details:");
    for model in structure.models.iter().take(5) {
        println!("    Model {}: {} atoms", model.serial, model.atoms.len());
    }
    if num_models > 5 {
        println!("    ... and {} more models", num_models - 5);
    }

    // Check if it's an NMR structure (typically 20+ models)
    if num_models >= 10 {
        println!("  [INFO] Appears to be an NMR ensemble");
    }

    // Check residue count in first model (should be ~20 for Trp-cage)
    let chains = structure.get_chain_ids();
    if !chains.is_empty() {
        let residues = structure.get_residues_for_chain(&chains[0]);
        println!("  Residues in chain {}: {}", chains[0], residues.len());
    }

    println!(
        "\n  Multi-model validation: {}",
        if passed { "PASSED" } else { "FAILED" }
    );
    Ok(passed)
}
