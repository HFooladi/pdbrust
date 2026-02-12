//! Detailed per-interface comparison of our DockQ vs reference implementation.
//!
//! Run with:
//! ```
//! cargo run --example dockq_detailed_comparison --features dockq
//! ```

use pdbrust::dockq::{ChainMappingStrategy, DockQOptions};
use pdbrust::parse_pdb_file;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let base = "tests/dockq_validation";

    println!("=== Detailed DockQ Comparison (per-interface) ===\n");

    // ========================================================================
    // Test 1: dimer_dimer - exact match expected
    // ========================================================================
    println!("=== Test 1: dimer_dimer (self-match) ===");
    {
        let model = parse_pdb_file(&format!("{}/dimer_model.pdb", base))?;
        let native = parse_pdb_file(&format!("{}/dimer_native.pdb", base))?;
        let result = model.dockq_to(&native)?;

        // Reference: 4 interfaces, all DockQ=1.000
        let expected = [
            ("AB", "AB", 1.000, 0.000, 0.000, 1.000),
            ("AL", "AL", 1.000, 0.000, 0.000, 1.000),
            ("AH", "AH", 1.000, 0.000, 0.000, 1.000),
            ("LH", "LH", 1.000, 0.000, 0.000, 1.000),
        ];

        println!(
            "Total DockQ: ours={:.3}, ref=1.000, diff={:.4}",
            result.total_dockq,
            (result.total_dockq - 1.000).abs()
        );
        println!("Interfaces found: {}, expected: 4\n", result.num_interfaces);

        println!(
            "{:<12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
            "Interface", "DockQ", "ref", "iRMSD", "ref", "LRMSD", "ref", "fnat", "ref"
        );
        println!("{}", "-".repeat(84));
        for (exp_mc, exp_nc, exp_dq, exp_ir, exp_lr, exp_fn) in &expected {
            // Find matching interface (check both chain orderings)
            let found = result.interfaces.iter().find(|i| {
                let mc = format!("{}{}", i.model_chains.0, i.model_chains.1);
                let nc = format!("{}{}", i.native_chains.0, i.native_chains.1);
                let mc_rev = format!("{}{}", i.model_chains.1, i.model_chains.0);
                let nc_rev = format!("{}{}", i.native_chains.1, i.native_chains.0);
                (mc == *exp_mc && nc == *exp_nc) || (mc_rev == *exp_mc && nc_rev == *exp_nc)
            });
            if let Some(iface) = found {
                println!(
                    "{}->{:<6} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3}",
                    exp_mc,
                    exp_nc,
                    iface.dockq,
                    exp_dq,
                    iface.irmsd,
                    exp_ir,
                    iface.lrmsd,
                    exp_lr,
                    iface.fnat,
                    exp_fn,
                );
            } else {
                println!("{}->{}: NOT FOUND", exp_mc, exp_nc);
            }
        }
    }

    // ========================================================================
    // Test 2: model vs native (single interface)
    // ========================================================================
    println!("\n=== Test 2: model.pdb vs native.pdb ===");
    {
        let model = parse_pdb_file(&format!("{}/model.pdb", base))?;
        let native = parse_pdb_file(&format!("{}/native.pdb", base))?;
        let options = DockQOptions {
            chain_mapping: ChainMappingStrategy::Explicit(vec![
                ("A".to_string(), "A".to_string()),
                ("B".to_string(), "B".to_string()),
            ]),
            ..Default::default()
        };
        let result = model.dockq_to_with_options(&native, options)?;

        let iface = &result.interfaces[0];
        println!(
            "{:<12} {:>10} {:>10} {:>10}",
            "", "Ours", "Reference", "Diff"
        );
        println!("{}", "-".repeat(45));
        println!(
            "{:<12} {:>10.3} {:>10.3} {:>10.4}",
            "DockQ",
            iface.dockq,
            0.700,
            (iface.dockq - 0.700).abs()
        );
        println!(
            "{:<12} {:>10.3} {:>10.3} {:>10.4}",
            "iRMSD",
            iface.irmsd,
            1.232,
            (iface.irmsd - 1.232).abs()
        );
        println!(
            "{:<12} {:>10.3} {:>10.3} {:>10.4}",
            "LRMSD",
            iface.lrmsd,
            1.516,
            (iface.lrmsd - 1.516).abs()
        );
        println!(
            "{:<12} {:>10.3} {:>10.3} {:>10.4}",
            "fnat",
            iface.fnat,
            0.533,
            (iface.fnat - 0.533).abs()
        );
        println!(
            "{:<12} {:>10.3} {:>10.3} {:>10.4}",
            "fnonnat",
            iface.fnonnat,
            0.238,
            (iface.fnonnat - 0.238).abs()
        );
        println!(
            "{:<12} {:>10.3} {:>10.3} {:>10.4}",
            "F1",
            iface.f1,
            0.627,
            (iface.f1 - 0.627).abs()
        );
    }

    // ========================================================================
    // Test 3: 1EXB octamer with explicit mapping (16 interfaces)
    // ========================================================================
    println!("\n=== Test 3: 1EXB octamer (explicit mapping) ===");
    {
        let model = parse_pdb_file(&format!("{}/1EXB_model.pdb", base))?;
        let native = parse_pdb_file(&format!("{}/1EXB_native.pdb", base))?;
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

        // Reference values for each interface
        let expected = [
            ("BC", "AD", 1.000, 0.000, 0.000, 1.000, 0.000),
            ("BD", "AC", 0.995, 0.000, 0.000, 0.985, 0.000),
            ("BF", "AE", 0.552, 1.458, 3.757, 0.304, 0.000),
            ("BH", "AG", 0.694, 0.793, 2.831, 0.400, 0.000),
            ("AC", "BD", 0.995, 0.000, 0.000, 0.985, 0.000),
            ("AD", "BC", 0.995, 0.000, 0.000, 0.985, 0.000),
            ("AE", "BF", 0.742, 0.850, 3.426, 0.609, 0.176),
            ("AG", "BH", 0.824, 0.625, 3.994, 0.800, 0.000),
            ("CF", "DE", 0.687, 1.161, 3.757, 0.600, 0.000),
            ("CG", "DH", 0.477, 1.425, 3.994, 0.087, 0.000),
            ("DH", "CG", 0.775, 0.711, 2.831, 0.609, 0.125),
            ("DE", "CF", 0.914, 0.547, 3.426, 1.000, 0.000),
            ("FH", "EG", 0.992, 0.001, 0.000, 0.976, 0.000),
            ("FG", "EH", 0.992, 0.001, 0.000, 0.976, 0.000),
            ("HE", "GF", 0.992, 0.001, 0.001, 0.976, 0.000),
            ("EG", "FH", 1.000, 0.001, 0.001, 1.000, 0.024),
        ];

        println!(
            "Total DockQ: ours={:.3}, ref=0.852, diff={:.4}",
            result.total_dockq,
            (result.total_dockq - 0.852).abs()
        );
        println!("Interfaces: {}/16\n", result.num_interfaces);

        println!(
            "{:<10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
            "Map", "DockQ", "ref", "dDQ", "iRMSD", "ref", "LRMSD", "ref", "fnat", "ref", "dfnat"
        );
        println!("{}", "-".repeat(110));

        let mut total_dq_diff = 0.0;
        let mut total_ir_diff = 0.0;
        let mut total_lr_diff = 0.0;
        let mut total_fn_diff = 0.0;
        let mut n_matched = 0;

        for (exp_mc, exp_nc, exp_dq, exp_ir, exp_lr, exp_fn, _exp_fnn) in &expected {
            let found = result.interfaces.iter().find(|i| {
                let mc = format!("{}{}", i.model_chains.0, i.model_chains.1);
                let nc = format!("{}{}", i.native_chains.0, i.native_chains.1);
                let mc_rev = format!("{}{}", i.model_chains.1, i.model_chains.0);
                let nc_rev = format!("{}{}", i.native_chains.1, i.native_chains.0);
                (mc == *exp_mc && nc == *exp_nc) || (mc_rev == *exp_mc && nc_rev == *exp_nc)
            });
            if let Some(iface) = found {
                let dq_diff = (iface.dockq - exp_dq).abs();
                let ir_diff = (iface.irmsd - exp_ir).abs();
                let lr_diff = (iface.lrmsd - exp_lr).abs();
                let fn_diff = (iface.fnat - exp_fn).abs();
                total_dq_diff += dq_diff;
                total_ir_diff += ir_diff;
                total_lr_diff += lr_diff;
                total_fn_diff += fn_diff;
                n_matched += 1;
                println!(
                    "{}->{:<5} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3}",
                    exp_mc,
                    exp_nc,
                    iface.dockq,
                    exp_dq,
                    dq_diff,
                    iface.irmsd,
                    exp_ir,
                    iface.lrmsd,
                    exp_lr,
                    iface.fnat,
                    exp_fn,
                    fn_diff,
                );
            } else {
                println!("{}->{}: NOT FOUND", exp_mc, exp_nc);
            }
        }

        if n_matched > 0 {
            println!(
                "\nMean absolute differences across {} interfaces:",
                n_matched
            );
            println!("  DockQ: {:.4}", total_dq_diff / n_matched as f64);
            println!("  iRMSD: {:.4}", total_ir_diff / n_matched as f64);
            println!("  LRMSD: {:.4}", total_lr_diff / n_matched as f64);
            println!("  fnat:  {:.4}", total_fn_diff / n_matched as f64);
        }
    }

    // ========================================================================
    // Test 4: 1A2K trimer with reference mapping BAC:ABC
    // ========================================================================
    println!("\n=== Test 4: 1A2K trimer (reference mapping BAC:ABC) ===");
    {
        let model = parse_pdb_file(&format!("{}/1A2K_model.pdb", base))?;
        let native = parse_pdb_file(&format!("{}/1A2K_native.pdb", base))?;

        // Reference mapping: BAC:ABC (model B->native A, model A->native B, model C->native C)
        let options = DockQOptions {
            chain_mapping: ChainMappingStrategy::Explicit(vec![
                ("B".to_string(), "A".to_string()),
                ("A".to_string(), "B".to_string()),
                ("C".to_string(), "C".to_string()),
            ]),
            ..Default::default()
        };
        let result = model.dockq_to_with_options(&native, options)?;

        let expected = [
            ("BA", "AB", 0.994, 0.000, 0.000, 0.983, 0.008),
            ("BC", "AC", 0.511, 1.237, 6.864, 0.333, 0.000),
            ("AC", "BC", 0.453, 2.104, 8.131, 0.500, 0.107),
        ];

        println!(
            "Total DockQ: ours={:.3}, ref=0.653, diff={:.4}",
            result.total_dockq,
            (result.total_dockq - 0.653).abs()
        );
        println!("Interfaces: {}/3\n", result.num_interfaces);

        println!(
            "{:<10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
            "Map", "DockQ", "ref", "dDQ", "iRMSD", "ref", "LRMSD", "ref", "fnat", "ref"
        );
        println!("{}", "-".repeat(100));

        for (exp_mc, exp_nc, exp_dq, exp_ir, exp_lr, exp_fn, _exp_fnn) in &expected {
            let found = result.interfaces.iter().find(|i| {
                let mc = format!("{}{}", i.model_chains.0, i.model_chains.1);
                let nc = format!("{}{}", i.native_chains.0, i.native_chains.1);
                let mc_rev = format!("{}{}", i.model_chains.1, i.model_chains.0);
                let nc_rev = format!("{}{}", i.native_chains.1, i.native_chains.0);
                (mc == *exp_mc && nc == *exp_nc) || (mc_rev == *exp_mc && nc_rev == *exp_nc)
            });
            if let Some(iface) = found {
                println!(
                    "{}->{:<5} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3}",
                    exp_mc,
                    exp_nc,
                    iface.dockq,
                    exp_dq,
                    (iface.dockq - exp_dq).abs(),
                    iface.irmsd,
                    exp_ir,
                    iface.lrmsd,
                    exp_lr,
                    iface.fnat,
                    exp_fn,
                );
            } else {
                println!("{}->{}: NOT FOUND", exp_mc, exp_nc);
            }
        }
    }

    Ok(())
}
