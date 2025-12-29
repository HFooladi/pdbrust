//! Rust benchmark to compare with Python libraryPDB.
//!
//! Run with: cargo run --release --features "filter,descriptors,quality,summary" --example rust_benchmark

use std::path::PathBuf;
use std::time::Instant;

fn get_test_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join(name)
}

fn benchmark<T, F: FnMut() -> T>(name: &str, mut f: F, iterations: usize) -> T {
    // Warmup
    let _ = f();

    let mut times = Vec::with_capacity(iterations);
    let mut result = None;

    for _ in 0..iterations {
        let start = Instant::now();
        result = Some(f());
        let elapsed = start.elapsed();
        times.push(elapsed.as_secs_f64() * 1000.0); // Convert to ms
    }

    let mean: f64 = times.iter().sum::<f64>() / times.len() as f64;
    let variance: f64 = times.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / times.len() as f64;
    let std = variance.sqrt();

    print!("  {:<20} {:8.3} ms ± {:.3} ms", name, mean, std);

    result.unwrap()
}

fn main() {
    use pdbrust::parse_pdb_file;

    println!("{}", "=".repeat(70));
    println!("Rust pdbrust Benchmark");
    println!("{}", "=".repeat(70));
    println!();

    let test_files = [("1UBQ", "1UBQ.pdb"), ("8HM2", "8HM2.pdb")];
    let iterations = 100;

    for (name, filename) in test_files.iter() {
        let path = get_test_file(filename);
        if !path.exists() {
            println!("Skipping {}: file not found", name);
            continue;
        }

        println!("Benchmarking: {}", name);
        println!("{}", "-".repeat(50));

        // Parsing benchmarks
        println!("\n[Parsing Operations]");

        let structure = benchmark("parse_pdb_file:", || parse_pdb_file(&path).unwrap(), iterations);
        println!("  ({} atoms)", structure.atoms.len());

        let chains = benchmark("get_chain_ids:", || structure.get_chain_ids(), iterations);
        println!("  ({:?})", chains);

        let residues = benchmark("get_residues:", || structure.get_residues(), iterations);
        println!("  ({} residues)", residues.len());

        #[cfg(feature = "filter")]
        {
            let ca_coords = benchmark("get_ca_coords:", || structure.get_ca_coords(None), iterations);
            println!("  ({} CA atoms)", ca_coords.len());
        }

        #[cfg(not(feature = "filter"))]
        {
            println!("  get_ca_coords:      (filter feature not enabled)");
        }

        // Descriptor benchmarks
        #[cfg(feature = "descriptors")]
        {
            println!("\n[Descriptor Operations]");

            let num_atoms = benchmark("num_atoms:", || structure.atoms.len(), iterations);
            println!("  ({})", num_atoms);

            let num_residues = benchmark("num_residues:", || structure.count_ca_residues(), iterations);
            println!("  ({})", num_residues);

            let _aa_comp = benchmark("aa_composition:", || structure.aa_composition(), iterations);
            println!();

            let gly_ratio = benchmark("glycine_ratio:", || structure.glycine_ratio(), iterations);
            println!("  ({:.4})", gly_ratio);

            let hydro_ratio = benchmark("hydrophobic_ratio:", || structure.hydrophobic_ratio(), iterations);
            println!("  ({:.4})", hydro_ratio);

            let rg = benchmark("radius_of_gyration:", || structure.radius_of_gyration(), iterations);
            println!("  ({:.4} Å)", rg);

            let max_dist = benchmark("max_ca_distance:", || structure.max_ca_distance(), iterations);
            println!("  ({:.4} Å)", max_dist);

            let missing = benchmark("missing_res_ratio:", || structure.missing_residue_ratio(), iterations);
            println!("  ({:.4})", missing);

            let ss_ratio = benchmark("ss_ratio:", || structure.secondary_structure_ratio(), iterations);
            println!("  ({:.4})", ss_ratio);

            let compact = benchmark("compactness_index:", || structure.compactness_index(), iterations);
            println!("  ({:.4})", compact);

            let density = benchmark("ca_density:", || structure.ca_density(), iterations);
            println!("  ({:.6})", density);
        }

        #[cfg(not(feature = "descriptors"))]
        {
            println!("\n[Descriptor Operations]");
            println!("  (descriptors feature not enabled)");
        }

        println!();
    }

    println!("{}", "=".repeat(70));
    println!("Benchmark complete.");
    println!("{}", "=".repeat(70));
}
