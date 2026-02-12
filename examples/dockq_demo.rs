//! DockQ interface quality assessment demo.
//!
//! Demonstrates DockQ scoring for protein-protein complex quality evaluation.
//!
//! Run with:
//! ```
//! cargo run --example dockq_demo --features dockq
//! ```

use pdbrust::PdbStructure;
use pdbrust::dockq::{ChainMappingStrategy, DockQOptions};
use pdbrust::records::Atom;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== DockQ Interface Quality Assessment Demo ===\n");

    // Create a synthetic two-chain complex (native)
    let native = create_complex();
    println!(
        "Native structure: {} atoms, chains: {:?}",
        native.atoms.len(),
        native.get_chain_ids()
    );

    // ========================================================================
    // 1. Self-comparison (perfect score)
    // ========================================================================
    println!("\n--- 1. Self-comparison ---");
    let result = native.dockq_to(&native)?;
    print_result(&result);

    // ========================================================================
    // 2. Slightly perturbed model
    // ========================================================================
    println!("\n--- 2. Slightly perturbed model (chain B shifted 1 A) ---");
    let mut perturbed = native.clone();
    for atom in &mut perturbed.atoms {
        if atom.chain_id == "B" {
            atom.y += 1.0;
        }
    }

    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };
    let result = perturbed.dockq_to_with_options(&native, options)?;
    print_result(&result);

    // ========================================================================
    // 3. Badly docked model
    // ========================================================================
    println!("\n--- 3. Badly docked model (chain B shifted 15 A) ---");
    let mut bad_model = native.clone();
    for atom in &mut bad_model.atoms {
        if atom.chain_id == "B" {
            atom.y += 15.0;
        }
    }

    let options = DockQOptions {
        chain_mapping: ChainMappingStrategy::Explicit(vec![
            ("A".to_string(), "A".to_string()),
            ("B".to_string(), "B".to_string()),
        ]),
        ..Default::default()
    };
    let result = bad_model.dockq_to_with_options(&native, options)?;
    print_result(&result);

    // ========================================================================
    // 4. Auto chain mapping with swapped labels
    // ========================================================================
    println!("\n--- 4. Auto chain mapping (model has swapped chain labels) ---");
    let mut swapped = PdbStructure::new();
    for atom in &native.atoms {
        let mut new_atom = atom.clone();
        new_atom.chain_id = if atom.chain_id == "A" {
            "B".to_string()
        } else {
            "A".to_string()
        };
        swapped.atoms.push(new_atom);
    }

    let result = swapped.dockq_to(&native)?;
    println!("Chain mapping: {:?}", result.chain_mapping);
    print_result(&result);

    println!("\n=== Demo complete ===");
    Ok(())
}

fn print_result(result: &pdbrust::dockq::DockQResult) {
    println!(
        "  Total DockQ: {:.4} ({} interface{})",
        result.total_dockq,
        result.num_interfaces,
        if result.num_interfaces != 1 { "s" } else { "" }
    );

    for iface in &result.interfaces {
        println!(
            "  Interface {}-{} (model {}-{}):",
            iface.native_chains.0,
            iface.native_chains.1,
            iface.model_chains.0,
            iface.model_chains.1,
        );
        println!("    DockQ:  {:.4} ({})", iface.dockq, iface.quality);
        println!("    fnat:   {:.4}", iface.fnat);
        println!("    iRMSD:  {:.4} A", iface.irmsd);
        println!("    LRMSD:  {:.4} A", iface.lrmsd);
        println!("    F1:     {:.4}", iface.f1);
        println!(
            "    Contacts: {} native, {} model",
            iface.num_native_contacts, iface.num_model_contacts
        );
    }
}

fn create_complex() -> PdbStructure {
    let mut structure = PdbStructure::new();
    let mut serial = 1;

    let chain_a_residues = ["ALA", "GLY", "VAL", "LEU", "ILE", "PHE", "TRP"];
    let chain_b_residues = ["SER", "THR", "ASP", "GLU", "LYS", "ARG", "HIS"];

    // Chain A along x-axis at y=0
    for (i, resname) in chain_a_residues.iter().enumerate() {
        let x = i as f64 * 3.8;
        let seq = (i + 1) as i32;
        for &(name, dx, dy, element) in &[
            ("N", -0.5, -0.5, "N"),
            ("CA", 0.0, 0.0, "C"),
            ("C", 0.5, 0.0, "C"),
            ("O", 0.5, 0.5, "O"),
        ] {
            structure.atoms.push(Atom {
                serial,
                name: name.to_string(),
                alt_loc: None,
                residue_name: resname.to_string(),
                chain_id: "A".to_string(),
                residue_seq: seq,
                ins_code: None,
                is_hetatm: false,
                x: x + dx,
                y: dy,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: element.to_string(),
            });
            serial += 1;
        }
    }

    // Chain B along x-axis at y=4.0 (within contact distance)
    for (i, resname) in chain_b_residues.iter().enumerate() {
        let x = i as f64 * 3.8;
        let seq = (i + 1) as i32;
        for &(name, dx, dy, element) in &[
            ("N", -0.5, 3.5, "N"),
            ("CA", 0.0, 4.0, "C"),
            ("C", 0.5, 4.0, "C"),
            ("O", 0.5, 4.5, "O"),
        ] {
            structure.atoms.push(Atom {
                serial,
                name: name.to_string(),
                alt_loc: None,
                residue_name: resname.to_string(),
                chain_id: "B".to_string(),
                residue_seq: seq,
                ins_code: None,
                is_hetatm: false,
                x: x + dx,
                y: dy,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: element.to_string(),
            });
            serial += 1;
        }
    }

    structure
}
