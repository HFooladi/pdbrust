//! DockQ validation against the original DockQ tool's test data.
//!
//! Downloads and compares results from:
//! https://github.com/wallnerlab/DockQ
//!
//! Run with:
//! ```
//! cargo run --example dockq_validation --features dockq
//! ```

use pdbrust::dockq::{ChainMappingStrategy, DockQOptions};
use pdbrust::{PdbStructure, parse_pdb_file};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let base = "tests/dockq_validation";

    println!("=== DockQ Validation Against Reference Implementation ===\n");

    // ========================================================================
    // Test 1: dimer_dimer (perfect match, expected DockQ = 1.000)
    // ========================================================================
    println!("--- Test 1: dimer_dimer (perfect match) ---");
    println!("Expected: Total DockQ = 1.000, 4 interfaces, all DockQ = 1.000\n");

    let model = parse_pdb_file(&format!("{}/dimer_model.pdb", base))?;
    let native = parse_pdb_file(&format!("{}/dimer_native.pdb", base))?;

    println!(
        "Model:  {} atoms, chains: {:?}",
        model.atoms.len(),
        model.get_chain_ids()
    );
    println!(
        "Native: {} atoms, chains: {:?}",
        native.atoms.len(),
        native.get_chain_ids()
    );

    let result = model.dockq_to(&native)?;
    println!("Chain mapping: {:?}", result.chain_mapping);
    println!(
        "Total DockQ: {:.3} over {} interfaces",
        result.total_dockq, result.num_interfaces
    );
    for iface in &result.interfaces {
        println!(
            "  DockQ {:.3} iRMSD {:.3} LRMSD {:.3} fnat {:.3} fnonnat {:.3} F1 {:.3} mapping {}{} -> {}{}",
            iface.dockq,
            iface.irmsd,
            iface.lrmsd,
            iface.fnat,
            iface.fnonnat,
            iface.f1,
            iface.model_chains.0,
            iface.model_chains.1,
            iface.native_chains.0,
            iface.native_chains.1,
        );
    }
    let pass1 = (result.total_dockq - 1.0).abs() < 0.01;
    println!("RESULT: {}\n", if pass1 { "PASS" } else { "FAIL" });

    // ========================================================================
    // Test 2: model.pdb vs native.pdb (expected DockQ = 0.700)
    // Ref: fnat=0.533, iRMSD=1.232, LRMSD=1.516
    // ========================================================================
    println!("--- Test 2: model.pdb vs native.pdb ---");
    println!("Expected: DockQ=0.700, fnat=0.533, iRMSD=1.232, LRMSD=1.516\n");

    let model = parse_pdb_file(&format!("{}/model.pdb", base))?;
    let native = parse_pdb_file(&format!("{}/native.pdb", base))?;

    println!(
        "Model:  {} atoms, chains: {:?}",
        model.atoms.len(),
        model.get_chain_ids()
    );
    println!(
        "Native: {} atoms, chains: {:?}",
        native.atoms.len(),
        native.get_chain_ids()
    );

    // The original DockQ uses --allowed_mismatches 1 for this case
    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };

    let result = model.dockq_to_with_options(&native, options)?;
    println!(
        "Total DockQ: {:.3} over {} interfaces",
        result.total_dockq, result.num_interfaces
    );
    for iface in &result.interfaces {
        println!(
            "  DockQ {:.3} iRMSD {:.3} LRMSD {:.3} fnat {:.3} fnonnat {:.3} F1 {:.3} mapping {},{} -> {},{}",
            iface.dockq,
            iface.irmsd,
            iface.lrmsd,
            iface.fnat,
            iface.fnonnat,
            iface.f1,
            iface.model_chains.0,
            iface.model_chains.1,
            iface.native_chains.0,
            iface.native_chains.1,
        );
    }
    println!("  Ref:  DockQ 0.700 iRMSD 1.232 LRMSD 1.516 fnat 0.533 fnonnat 0.238 F1 0.627");
    let pass2 = (result.total_dockq - 0.700).abs() < 0.10; // 10% tolerance for now
    println!(
        "RESULT: {} (tolerance ±0.10)\n",
        if pass2 { "PASS" } else { "FAIL" }
    );

    // ========================================================================
    // Test 3: 1A2K trimeric complex (expected total DockQ = 0.653)
    // Interface A,B: DockQ=0.994, fnat=0.983
    // Interface A,C: DockQ=0.511, fnat=0.333
    // Interface B,C: DockQ=0.453, fnat=0.500
    // ========================================================================
    println!("--- Test 3: 1A2K trimeric complex ---");
    println!("Expected: Total DockQ=0.653, 3 interfaces\n");

    let model = parse_pdb_file(&format!("{}/1A2K_model.pdb", base))?;
    let native = parse_pdb_file(&format!("{}/1A2K_native.pdb", base))?;

    println!(
        "Model:  {} atoms, chains: {:?}",
        model.atoms.len(),
        model.get_chain_ids()
    );
    println!(
        "Native: {} atoms, chains: {:?}",
        native.atoms.len(),
        native.get_chain_ids()
    );

    let result = model.dockq_to(&native)?;
    println!("Chain mapping: {:?}", result.chain_mapping);
    println!(
        "Total DockQ: {:.3} over {} interfaces",
        result.total_dockq, result.num_interfaces
    );
    for iface in &result.interfaces {
        println!(
            "  DockQ {:.3} iRMSD {:.3} LRMSD {:.3} fnat {:.3} fnonnat {:.3} F1 {:.3} mapping {},{} -> {},{}",
            iface.dockq,
            iface.irmsd,
            iface.lrmsd,
            iface.fnat,
            iface.fnonnat,
            iface.f1,
            iface.model_chains.0,
            iface.model_chains.1,
            iface.native_chains.0,
            iface.native_chains.1,
        );
    }
    println!("  Ref interfaces:");
    println!("    A,B -> B,A: DockQ 0.994 iRMSD 0.000 LRMSD 0.000 fnat 0.983");
    println!("    A,C -> B,C: DockQ 0.511 iRMSD 1.237 LRMSD 6.864 fnat 0.333");
    println!("    B,C -> A,C: DockQ 0.453 iRMSD 2.104 LRMSD 8.131 fnat 0.500");
    let pass3 = (result.total_dockq - 0.653).abs() < 0.15;
    println!(
        "RESULT: {} (tolerance ±0.15)\n",
        if pass3 { "PASS" } else { "FAIL" }
    );

    // ========================================================================
    // Test 4: 1EXB octameric complex (expected total DockQ = 0.852)
    // ========================================================================
    println!("--- Test 4: 1EXB octameric complex ---");
    println!("Expected: Total DockQ=0.852, 16 interfaces\n");

    let model = parse_pdb_file(&format!("{}/1EXB_model.pdb", base))?;
    let native = parse_pdb_file(&format!("{}/1EXB_native.pdb", base))?;

    println!(
        "Model:  {} atoms, chains: {:?}",
        model.atoms.len(),
        model.get_chain_ids()
    );
    println!(
        "Native: {} atoms, chains: {:?}",
        native.atoms.len(),
        native.get_chain_ids()
    );

    // This is a big complex - use explicit mapping from the reference
    // Reference mapping: BACDFHEG:ABDCEGFH (model:native)
    // model B -> native A, model A -> native B, model C -> native D, etc.
    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("B".to_string(), "A".to_string()),
            ("A".to_string(), "B".to_string()),
            ("C".to_string(), "D".to_string()),
            ("D".to_string(), "C".to_string()),
            ("F".to_string(), "E".to_string()),
            ("H".to_string(), "G".to_string()),
            ("E".to_string(), "F".to_string()),
            ("G".to_string(), "H".to_string()),
        ]),
        ..Default::default()
    };

    let result = model.dockq_to_with_options(&native, options)?;
    println!(
        "Total DockQ: {:.3} over {} interfaces",
        result.total_dockq, result.num_interfaces
    );
    for iface in &result.interfaces {
        println!(
            "  DockQ {:.3} iRMSD {:.3} LRMSD {:.3} fnat {:.3} mapping {}{} -> {}{}",
            iface.dockq,
            iface.irmsd,
            iface.lrmsd,
            iface.fnat,
            iface.model_chains.0,
            iface.model_chains.1,
            iface.native_chains.0,
            iface.native_chains.1,
        );
    }
    let pass4 = (result.total_dockq - 0.852).abs() < 0.15;
    println!(
        "RESULT: {} (tolerance ±0.15)\n",
        if pass4 { "PASS" } else { "FAIL" }
    );

    // ========================================================================
    // Summary
    // ========================================================================
    println!("=== Summary ===");
    println!(
        "Test 1 (dimer perfect match): {}",
        if pass1 { "PASS" } else { "FAIL" }
    );
    println!(
        "Test 2 (model/native AB):     {}",
        if pass2 { "PASS" } else { "FAIL" }
    );
    println!(
        "Test 3 (1A2K trimer):         {}",
        if pass3 { "PASS" } else { "FAIL" }
    );
    println!(
        "Test 4 (1EXB octamer):        {}",
        if pass4 { "PASS" } else { "FAIL" }
    );
    println!(
        "\n{}/{} tests passed",
        [pass1, pass2, pass3, pass4].iter().filter(|x| **x).count(),
        4
    );

    Ok(())
}
