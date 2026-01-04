use criterion::{Criterion, black_box, criterion_group, criterion_main};
use pdbrust::PdbStructure;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_pdb(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file
}

#[allow(clippy::too_many_arguments)]
fn generate_atom_record(
    serial: i32,
    name: &str,
    residue_name: &str,
    chain_id: &str,
    residue_seq: i32,
    x: f64,
    y: f64,
    z: f64,
) -> String {
    format!(
        "ATOM  {:5} {:4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           {:>2}  \n",
        serial,
        name,
        residue_name,
        chain_id,
        residue_seq,
        x,
        y,
        z,
        name.chars().next().unwrap()
    )
}

fn generate_large_structure(num_atoms: usize) -> String {
    let mut content = String::with_capacity(num_atoms * 80);
    for i in 0..num_atoms {
        content.push_str(&generate_atom_record(
            i as i32 + 1,
            "N",
            "ALA",
            "A",
            (i / 3) as i32 + 1,
            0.0,
            0.0,
            0.0,
        ));
    }
    content
}

fn generate_multi_model_structure(num_models: usize, atoms_per_model: usize) -> String {
    let mut content = String::with_capacity(num_models * atoms_per_model * 82);
    for model in 1..=num_models {
        content.push_str(&format!("MODEL {:8}\n", model));
        for atom in 1..=atoms_per_model {
            content.push_str(&generate_atom_record(
                atom as i32,
                "N",
                "ALA",
                "A",
                (atom / 3) as i32 + 1,
                0.0,
                0.0,
                0.0,
            ));
        }
        content.push_str("ENDMDL\n");
    }
    content
}

fn benchmark_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("parsing");

    // Benchmark parsing different sizes of structures
    for size in [100, 1000, 10000].iter() {
        group.bench_function(format!("parse_{}_atoms", size), |b| {
            let content = generate_large_structure(*size);
            let file = create_test_pdb(&content);
            b.iter(|| {
                black_box(PdbStructure::from_file(file.path()).unwrap());
            });
        });
    }

    // Benchmark parsing multiple models
    for (models, atoms) in [(2, 100), (5, 100), (10, 100)].iter() {
        group.bench_function(format!("parse_{}_models_{}_atoms", models, atoms), |b| {
            let content = generate_multi_model_structure(*models, *atoms);
            let file = create_test_pdb(&content);
            b.iter(|| {
                black_box(PdbStructure::from_file(file.path()).unwrap());
            });
        });
    }

    group.finish();
}

fn benchmark_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("operations");

    // Create a test structure
    let content = generate_large_structure(1000);
    let file = create_test_pdb(&content);
    let structure = PdbStructure::from_file(file.path()).unwrap();

    // Benchmark chain operations
    group.bench_function("get_chain_ids", |b| {
        b.iter(|| {
            black_box(structure.get_chain_ids());
        });
    });

    // Benchmark residue operations
    group.bench_function("get_residues_for_chain", |b| {
        b.iter(|| {
            black_box(structure.get_residues_for_chain("A"));
        });
    });

    // Benchmark atom connectivity operations
    group.bench_function("get_connected_atoms", |b| {
        b.iter(|| {
            black_box(structure.get_connected_atoms(1));
        });
    });

    group.finish();
}

criterion_group!(benches, benchmark_parsing, benchmark_operations);
criterion_main!(benches);
