#!/usr/bin/env python3
"""
Secondary structure (DSSP) example for PDBRust Python bindings.

This example demonstrates:
- DSSP-like secondary structure assignment
- Secondary structure composition (helix/sheet/coil fractions)
- Compact string representation
- Per-residue secondary structure assignments
- Secondary structure code interpretation
- Finding helices and strands in a structure
"""

import pdbrust
from pdbrust import SecondaryStructure

# Sample PDB file
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust Secondary Structure (DSSP) Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms, {structure.num_residues} residues")

    # --- Basic Secondary Structure Assignment ---
    print("\n1. SECONDARY STRUCTURE ASSIGNMENT")
    print("-" * 40)

    # Compute DSSP-like secondary structure assignment
    ss = structure.assign_secondary_structure()

    print(f"Total residues assigned: {len(ss)}")
    print(f"Helix residues: {ss.helix_count}")
    print(f"Sheet residues: {ss.sheet_count}")
    print(f"Coil residues: {ss.coil_count}")

    # --- Composition Fractions ---
    print("\n2. SECONDARY STRUCTURE COMPOSITION")
    print("-" * 40)

    # Using properties
    print(f"Helix fraction: {ss.helix_fraction * 100:.1f}%")
    print(f"Sheet fraction: {ss.sheet_fraction * 100:.1f}%")
    print(f"Coil fraction:  {ss.coil_fraction * 100:.1f}%")

    # Using the composition method (returns tuple)
    helix, sheet, coil = ss.composition()
    print(f"\nUsing composition() method:")
    print(f"  ({helix:.3f}, {sheet:.3f}, {coil:.3f})")

    # Using structure method directly
    h2, s2, c2 = structure.secondary_structure_composition()
    print(f"\nDirect structure method:")
    print(f"  Helix: {h2:.3f}, Sheet: {s2:.3f}, Coil: {c2:.3f}")

    # --- String Representation ---
    print("\n3. SECONDARY STRUCTURE STRING")
    print("-" * 40)

    # Get compact string representation
    ss_string = ss.to_string()
    print(f"Length: {len(ss_string)} characters")
    print(f"\nSecondary structure sequence:")
    # Print in chunks for readability
    chunk_size = 50
    for i in range(0, len(ss_string), chunk_size):
        print(f"  {i+1:3d}: {ss_string[i:i+chunk_size]}")

    # Alternative: use structure method directly
    ss_direct = structure.secondary_structure_string()
    print(f"\nDirect method: {ss_direct[:50]}...")

    # --- DSSP Codes Reference ---
    print("\n4. DSSP CODE REFERENCE")
    print("-" * 40)

    print("""
DSSP Codes (single character per residue):
  H - Alpha helix (i -> i+4 H-bond)
  G - 3-10 helix (i -> i+3 H-bond)
  I - Pi helix (i -> i+5 H-bond)
  P - Kappa helix / PPII (polyproline II)
  E - Extended strand (beta-sheet)
  B - Isolated beta-bridge
  T - Hydrogen-bonded turn
  S - Bend (high backbone curvature)
  C - Coil (none of the above)

Categories:
  Helix: H, G, I, P
  Sheet: E, B
  Coil:  T, S, C
""")

    # --- Per-Residue Assignments ---
    print("\n5. PER-RESIDUE ASSIGNMENTS")
    print("-" * 40)

    print(f"{'Chain':<6} {'ResSeq':<8} {'Name':<5} {'SS Code':<8}")
    print("-" * 30)

    # Print first 15 residues
    for res in ss.residue_assignments[:15]:
        print(f"{res.chain_id:<6} {res.residue_seq:<8} {res.residue_name:<5} {res.code():<8}")

    print(f"... ({len(ss) - 15} more residues)")

    # --- Iterating Over Assignments ---
    print("\n6. ITERATING OVER ASSIGNMENTS")
    print("-" * 40)

    # Count each secondary structure type
    ss_counts = {}
    for res in ss:  # Uses __iter__
        code = res.code()
        ss_counts[code] = ss_counts.get(code, 0) + 1

    print("Detailed SS breakdown:")
    for code in sorted(ss_counts.keys()):
        count = ss_counts[code]
        pct = count / len(ss) * 100
        print(f"  {code}: {count:3d} residues ({pct:5.1f}%)")

    # --- Finding Helices and Strands ---
    print("\n7. FINDING HELICES AND STRANDS")
    print("-" * 40)

    # Find all helical regions
    helices = []
    current_helix = []

    for res in ss:
        if res.ss.is_helix():
            current_helix.append(res)
        else:
            if len(current_helix) >= 3:  # Minimum helix length
                helices.append(current_helix)
            current_helix = []

    if len(current_helix) >= 3:
        helices.append(current_helix)

    print(f"Found {len(helices)} helix segments:")
    for i, helix in enumerate(helices, 1):
        start = helix[0]
        end = helix[-1]
        print(f"  Helix {i}: {start.chain_id}{start.residue_seq}-{end.residue_seq} "
              f"({len(helix)} residues)")

    # Find all strand regions
    strands = []
    current_strand = []

    for res in ss:
        if res.ss.is_sheet():
            current_strand.append(res)
        else:
            if len(current_strand) >= 2:  # Minimum strand length
                strands.append(current_strand)
            current_strand = []

    if len(current_strand) >= 2:
        strands.append(current_strand)

    print(f"\nFound {len(strands)} strand segments:")
    for i, strand in enumerate(strands, 1):
        start = strand[0]
        end = strand[-1]
        print(f"  Strand {i}: {start.chain_id}{start.residue_seq}-{end.residue_seq} "
              f"({len(strand)} residues)")

    # --- Using SecondaryStructure Class Attributes ---
    print("\n8. SECONDARY STRUCTURE CLASS")
    print("-" * 40)

    # Access predefined SS types
    print("Predefined SecondaryStructure types:")
    print(f"  ALPHA_HELIX: code = '{SecondaryStructure.ALPHA_HELIX.code()}'")
    print(f"  HELIX_310:   code = '{SecondaryStructure.HELIX_310.code()}'")
    print(f"  PI_HELIX:    code = '{SecondaryStructure.PI_HELIX.code()}'")
    print(f"  KAPPA_HELIX: code = '{SecondaryStructure.KAPPA_HELIX.code()}'")
    print(f"  EXTENDED_STRAND: code = '{SecondaryStructure.EXTENDED_STRAND.code()}'")
    print(f"  BETA_BRIDGE: code = '{SecondaryStructure.BETA_BRIDGE.code()}'")
    print(f"  TURN:        code = '{SecondaryStructure.TURN.code()}'")
    print(f"  BEND:        code = '{SecondaryStructure.BEND.code()}'")
    print(f"  COIL:        code = '{SecondaryStructure.COIL.code()}'")

    # Check properties
    alpha = SecondaryStructure.ALPHA_HELIX
    print(f"\nALPHA_HELIX properties:")
    print(f"  is_helix(): {alpha.is_helix()}")
    print(f"  is_sheet(): {alpha.is_sheet()}")
    print(f"  is_coil():  {alpha.is_coil()}")

    # --- Accessing by Index ---
    print("\n9. INDEXING ASSIGNMENTS")
    print("-" * 40)

    # Access by positive index
    first_res = ss[0]
    print(f"First residue: {first_res}")

    # Access by negative index
    last_res = ss[-1]
    print(f"Last residue: {last_res}")

    # Accessing specific index
    mid_index = len(ss) // 2
    mid_res = ss[mid_index]
    print(f"Middle residue (index {mid_index}): {mid_res}")

    # --- Warnings ---
    print("\n10. ASSIGNMENT WARNINGS")
    print("-" * 40)

    if ss.warnings:
        print(f"Warnings generated during assignment ({len(ss.warnings)}):")
        for warning in ss.warnings[:5]:
            print(f"  - {warning}")
        if len(ss.warnings) > 5:
            print(f"  ... and {len(ss.warnings) - 5} more")
    else:
        print("No warnings generated during assignment.")

    # --- Summary ---
    print("\n11. SUMMARY")
    print("-" * 40)

    print(f"""
Structure: {PDB_FILE}
Total residues: {len(ss)}

Secondary Structure Composition:
  Helix: {ss.helix_fraction * 100:5.1f}% ({ss.helix_count} residues)
  Sheet: {ss.sheet_fraction * 100:5.1f}% ({ss.sheet_count} residues)
  Coil:  {ss.coil_fraction * 100:5.1f}% ({ss.coil_count} residues)

Structural Elements:
  Alpha helices: {len(helices)}
  Beta strands:  {len(strands)}

SS String (first 50 chars): {ss_string[:50]}
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
