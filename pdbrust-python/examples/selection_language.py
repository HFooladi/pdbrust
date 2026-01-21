#!/usr/bin/env python3
"""
Selection language example for PDBRust Python bindings.

This example demonstrates:
- PyMOL/VMD-style selection language
- Basic selections (chain, name, resname, resid)
- Range selections (resid 1:100)
- Keyword selections (backbone, protein, water, hetero)
- Numeric comparisons (bfactor, occupancy)
- Boolean operators (and, or, not, parentheses)
- Complex multi-criteria queries
"""

import pdbrust

# Sample PDB file
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust Selection Language Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms, {structure.num_residues} residues")

    # --- Basic Selections ---
    print("\n1. BASIC SELECTIONS")
    print("-" * 40)

    # Select by chain
    chain_a = structure.select("chain A")
    print(f"chain A: {chain_a.num_atoms} atoms")

    # Select by atom name
    ca_atoms = structure.select("name CA")
    print(f"name CA: {ca_atoms.num_atoms} atoms")

    # Select by residue name
    ala_residues = structure.select("resname ALA")
    print(f"resname ALA: {ala_residues.num_atoms} atoms")

    # Select by residue number
    res50 = structure.select("resid 50")
    print(f"resid 50: {res50.num_atoms} atoms")

    # Select by element
    carbon_atoms = structure.select("element C")
    print(f"element C: {carbon_atoms.num_atoms} atoms")

    # --- Range Selections ---
    print("\n2. RANGE SELECTIONS")
    print("-" * 40)

    # Select residue range (inclusive)
    first_10 = structure.select("resid 1:10")
    print(f"resid 1:10: {first_10.num_atoms} atoms")

    mid_range = structure.select("resid 20:40")
    print(f"resid 20:40: {mid_range.num_atoms} atoms")

    # --- Keyword Selections ---
    print("\n3. KEYWORD SELECTIONS")
    print("-" * 40)

    # Backbone atoms (N, CA, C, O)
    backbone = structure.select("backbone")
    print(f"backbone: {backbone.num_atoms} atoms")

    # Protein atoms (standard amino acids)
    protein = structure.select("protein")
    print(f"protein: {protein.num_atoms} atoms")

    # Select all atoms
    all_atoms = structure.select("all")
    print(f"all: {all_atoms.num_atoms} atoms")

    # Alternative: select all with wildcard
    all_atoms2 = structure.select("*")
    print(f"* (wildcard): {all_atoms2.num_atoms} atoms")

    # Hydrogen atoms
    hydrogens = structure.select("hydrogen")
    print(f"hydrogen: {hydrogens.num_atoms} atoms")

    # HETATM records (ligands, waters)
    hetero = structure.select("hetero")
    print(f"hetero: {hetero.num_atoms} atoms")

    # Water molecules
    water = structure.select("water")
    print(f"water: {water.num_atoms} atoms")

    # --- Boolean Operators ---
    print("\n4. BOOLEAN OPERATORS")
    print("-" * 40)

    # AND: combine conditions
    chain_a_ca = structure.select("chain A and name CA")
    print(f"chain A and name CA: {chain_a_ca.num_atoms} atoms")

    # OR: either condition
    ala_or_gly = structure.select("resname ALA or resname GLY")
    print(f"resname ALA or resname GLY: {ala_or_gly.num_atoms} atoms")

    # NOT: negation
    not_backbone = structure.select("not backbone")
    print(f"not backbone: {not_backbone.num_atoms} atoms")

    # Combine NOT with other operators
    non_hydrogen = structure.select("protein and not hydrogen")
    print(f"protein and not hydrogen: {non_hydrogen.num_atoms} atoms")

    # --- Parentheses for Grouping ---
    print("\n5. PARENTHESES FOR GROUPING")
    print("-" * 40)

    # Group conditions with parentheses
    grouped = structure.select("(resname ALA or resname GLY) and backbone")
    print(f"(resname ALA or resname GLY) and backbone: {grouped.num_atoms} atoms")

    # Complex grouping
    complex_sel = structure.select("(chain A or chain B) and name CA")
    print(f"(chain A or chain B) and name CA: {complex_sel.num_atoms} atoms")

    # --- Numeric Comparisons ---
    print("\n6. NUMERIC COMPARISONS")
    print("-" * 40)

    # B-factor comparisons
    low_bfactor = structure.select("bfactor < 20.0")
    print(f"bfactor < 20.0: {low_bfactor.num_atoms} atoms")

    high_bfactor = structure.select("bfactor > 30.0")
    print(f"bfactor > 30.0: {high_bfactor.num_atoms} atoms")

    # B-factor range
    mid_bfactor = structure.select("bfactor >= 15.0 and bfactor <= 25.0")
    print(f"bfactor >= 15.0 and bfactor <= 25.0: {mid_bfactor.num_atoms} atoms")

    # Occupancy comparisons
    full_occupancy = structure.select("occupancy >= 1.0")
    print(f"occupancy >= 1.0: {full_occupancy.num_atoms} atoms")

    partial_occupancy = structure.select("occupancy < 1.0")
    print(f"occupancy < 1.0: {partial_occupancy.num_atoms} atoms")

    # --- Complex Queries ---
    print("\n7. COMPLEX QUERIES")
    print("-" * 40)

    # High-quality protein backbone atoms
    hq_backbone = structure.select("protein and backbone and bfactor < 30.0")
    print(f"protein and backbone and bfactor < 30.0: {hq_backbone.num_atoms} atoms")

    # First 20 residues CA atoms
    first_20_ca = structure.select("resid 1:20 and name CA")
    print(f"resid 1:20 and name CA: {first_20_ca.num_atoms} atoms")

    # Hydrophobic core residues (example: ALA, VAL, ILE, LEU, MET, PHE)
    hydrophobic = structure.select(
        "(resname ALA or resname VAL or resname ILE or resname LEU or resname MET or resname PHE) and name CA"
    )
    print(f"Hydrophobic CA atoms: {hydrophobic.num_atoms} atoms")

    # --- Practical Use Case ---
    print("\n8. PRACTICAL USE CASE: ANALYSIS SUBSET")
    print("-" * 40)

    # Select clean subset for analysis:
    # - Protein only (no ligands)
    # - No hydrogens
    # - Backbone atoms
    # - Low B-factor (well-ordered)
    analysis_subset = structure.select("protein and backbone and not hydrogen and bfactor < 40.0")
    print(f"Analysis subset: {analysis_subset.num_atoms} atoms")

    # Get info about the selection
    print("\nAnalysis subset details:")
    print(f"  Chains: {analysis_subset.get_chain_ids()}")
    print(f"  Residues: {analysis_subset.num_residues}")

    # --- Method Chaining Alternative ---
    print("\n9. SELECTION vs METHOD CHAINING")
    print("-" * 40)

    # Using selection language
    sel_result = structure.select("chain A and name CA")
    print(f"Selection 'chain A and name CA': {sel_result.num_atoms} atoms")

    # Equivalent using method chaining
    chain_result = structure.keep_only_chain("A").keep_only_ca()
    print(f"Method chaining equivalent: {chain_result.num_atoms} atoms")

    # Selection language is more flexible for complex queries!

    # --- Selection Syntax Reference ---
    print("\n10. SELECTION SYNTAX REFERENCE")
    print("-" * 40)

    print("""
Basic selectors:
  chain A          - atoms in chain A
  name CA          - atoms named CA
  resname ALA      - atoms in alanine residues
  resid 50         - atoms in residue 50
  resid 1:100      - atoms in residues 1-100
  element C        - carbon atoms

Keywords:
  backbone         - N, CA, C, O atoms
  protein          - standard amino acids
  nucleic          - standard nucleotides
  water            - water molecules (HOH, WAT)
  hetero           - HETATM records
  hydrogen         - hydrogen atoms
  all / *          - all atoms

Numeric comparisons:
  bfactor < 30.0   - B-factor less than 30
  occupancy >= 0.5 - occupancy at least 0.5

Boolean operators:
  and              - both conditions
  or               - either condition
  not              - negation
  ()               - grouping
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
