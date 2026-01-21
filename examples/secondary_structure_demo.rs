//! DSSP secondary structure assignment example
//!
//! This example demonstrates:
//! - Computing secondary structure assignment using DSSP 4-like algorithm
//! - Displaying helix/sheet/coil fractions
//! - Generating compact SS string representation
//! - Iterating over per-residue assignments
//!
//! Run with: cargo run --example secondary_structure_demo --features dssp

use pdbrust::parse_pdb_file;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("================================================");
    println!("PDBRust Secondary Structure Demo (DSSP)");
    println!("================================================\n");

    // Load structure
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pdb_files")
        .join("1UBQ.pdb");

    println!("Loading structure from {:?}", pdb_path);
    let structure = parse_pdb_file(&pdb_path)?;
    println!("Total atoms: {}", structure.atoms.len());
    println!("Total residues: {}", structure.get_residues().len());
    println!();

    // ====================
    // Secondary Structure Assignment
    // ====================
    println!("=== Secondary Structure Assignment ===");

    let ss = structure.assign_secondary_structure();
    println!("Residues assigned: {}", ss.residue_assignments.len());
    println!();

    // ====================
    // Overall Composition
    // ====================
    println!("=== Secondary Structure Composition ===");

    println!("Helix:  {:.1}%", ss.helix_fraction * 100.0);
    println!("Sheet:  {:.1}%", ss.sheet_fraction * 100.0);
    println!("Coil:   {:.1}%", ss.coil_fraction * 100.0);
    println!();

    // Alternative: use convenience method
    let (helix, sheet, coil) = structure.secondary_structure_composition();
    println!("Via convenience method:");
    println!(
        "  Helix: {:.1}%, Sheet: {:.1}%, Coil: {:.1}%",
        helix * 100.0,
        sheet * 100.0,
        coil * 100.0
    );
    println!();

    // ====================
    // Compact String Representation
    // ====================
    println!("=== Compact String Representation ===");

    let ss_string = structure.secondary_structure_string();
    println!("SS String (first 60 chars):");
    if ss_string.len() > 60 {
        println!("  {}...", &ss_string[..60]);
    } else {
        println!("  {}", ss_string);
    }
    println!("  Length: {} characters", ss_string.len());
    println!();

    // Count each type
    let h_count = ss_string.chars().filter(|&c| c == 'H').count();
    let g_count = ss_string.chars().filter(|&c| c == 'G').count();
    let i_count = ss_string.chars().filter(|&c| c == 'I').count();
    let p_count = ss_string.chars().filter(|&c| c == 'P').count();
    let e_count = ss_string.chars().filter(|&c| c == 'E').count();
    let b_count = ss_string.chars().filter(|&c| c == 'B').count();
    let t_count = ss_string.chars().filter(|&c| c == 'T').count();
    let s_count = ss_string.chars().filter(|&c| c == 'S').count();
    let c_count = ss_string.chars().filter(|&c| c == 'C').count();

    println!("Character counts:");
    println!("  H (α-helix):     {}", h_count);
    println!("  G (3₁₀-helix):   {}", g_count);
    println!("  I (π-helix):     {}", i_count);
    println!("  P (PPII/κ):      {}", p_count);
    println!("  E (strand):      {}", e_count);
    println!("  B (β-bridge):    {}", b_count);
    println!("  T (turn):        {}", t_count);
    println!("  S (bend):        {}", s_count);
    println!("  C (coil):        {}", c_count);
    println!();

    // ====================
    // Per-Residue Assignments
    // ====================
    println!("=== Per-Residue Assignments ===");

    println!("First 15 residues:");
    println!("{:<6} {:<8} {:<5} {:<6}", "Chain", "ResSeq", "Name", "SS");
    println!("{}", "-".repeat(30));

    for res in ss.residue_assignments.iter().take(15) {
        println!(
            "{:<6} {:<8} {:<5} {:<6}",
            res.chain_id,
            res.residue_seq,
            res.residue_name,
            res.ss.code()
        );
    }
    println!();

    // ====================
    // Find Secondary Structure Elements
    // ====================
    println!("=== Secondary Structure Elements ===");

    // Find helices (consecutive H residues)
    let mut in_helix = false;
    let mut helix_start = 0;
    let mut helices: Vec<(usize, usize)> = Vec::new();

    for (i, res) in ss.residue_assignments.iter().enumerate() {
        let is_helix = res.ss.code() == 'H' || res.ss.code() == 'G';
        if is_helix && !in_helix {
            in_helix = true;
            helix_start = i;
        } else if !is_helix && in_helix {
            in_helix = false;
            if i - helix_start >= 3 {
                // Minimum 3 residues for a helix
                helices.push((helix_start, i - 1));
            }
        }
    }
    if in_helix && ss.residue_assignments.len() - helix_start >= 3 {
        helices.push((helix_start, ss.residue_assignments.len() - 1));
    }

    println!("Helices found: {}", helices.len());
    for (start, end) in helices.iter().take(5) {
        let start_res = &ss.residue_assignments[*start];
        let end_res = &ss.residue_assignments[*end];
        println!(
            "  {}{}-{}{}: {} residues",
            start_res.chain_id,
            start_res.residue_seq,
            end_res.chain_id,
            end_res.residue_seq,
            end - start + 1
        );
    }
    println!();

    // Find strands (consecutive E residues)
    let mut in_strand = false;
    let mut strand_start = 0;
    let mut strands: Vec<(usize, usize)> = Vec::new();

    for (i, res) in ss.residue_assignments.iter().enumerate() {
        let is_strand = res.ss.code() == 'E';
        if is_strand && !in_strand {
            in_strand = true;
            strand_start = i;
        } else if !is_strand && in_strand {
            in_strand = false;
            if i - strand_start >= 2 {
                // Minimum 2 residues for a strand
                strands.push((strand_start, i - 1));
            }
        }
    }
    if in_strand && ss.residue_assignments.len() - strand_start >= 2 {
        strands.push((strand_start, ss.residue_assignments.len() - 1));
    }

    println!("β-strands found: {}", strands.len());
    for (start, end) in strands.iter().take(5) {
        let start_res = &ss.residue_assignments[*start];
        let end_res = &ss.residue_assignments[*end];
        println!(
            "  {}{}-{}{}: {} residues",
            start_res.chain_id,
            start_res.residue_seq,
            end_res.chain_id,
            end_res.residue_seq,
            end - start + 1
        );
    }
    println!();

    // ====================
    // DSSP Codes Reference
    // ====================
    println!("=== DSSP Secondary Structure Codes ===");
    println!(
        "
Code  Type           Description
----  -------------  ----------------------------------------
H     α-helix        Main helix type (i → i+4 H-bond pattern)
G     3₁₀-helix      Tighter helix (i → i+3 H-bond pattern)
I     π-helix        Wider helix (i → i+5 H-bond pattern)
P     κ-helix/PPII   Polyproline II (dihedral-based)
E     Extended       β-strand in sheet
B     β-bridge       Isolated β-bridge
T     Turn           Hydrogen-bonded turn
S     Bend           High backbone curvature
C     Coil           None of the above
"
    );

    // ====================
    // Summary
    // ====================
    println!("=== Summary ===");
    println!(
        "
Secondary structure analysis for {}:
  - Total residues: {}
  - Helix content: {:.1}% ({} helices)
  - Sheet content: {:.1}% ({} strands)
  - Coil content:  {:.1}%
",
        pdb_path.file_name().unwrap().to_string_lossy(),
        ss.residue_assignments.len(),
        ss.helix_fraction * 100.0,
        helices.len(),
        ss.sheet_fraction * 100.0,
        strands.len(),
        ss.coil_fraction * 100.0
    );

    println!("================================================");
    println!("Secondary structure demo completed successfully!");
    println!("================================================");

    Ok(())
}
