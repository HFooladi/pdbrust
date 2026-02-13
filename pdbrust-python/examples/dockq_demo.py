"""DockQ interface quality assessment demo.

Demonstrates DockQ scoring for evaluating protein-protein docking quality.

Usage:
    cd pdbrust-python
    maturin develop --release
    python examples/dockq_demo.py
"""

import os
import sys

import pdbrust
from pdbrust import (
    ChainMappingStrategy,
    DockQOptions,
    DockQQuality,
)


def main():
    print("=== DockQ Interface Quality Assessment (Python) ===\n")

    # Find the 1HSG PDB file (HIV protease dimer with chains A and B)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(script_dir, "..", "..", "examples", "pdb_files", "1HSG.pdb")

    if not os.path.exists(pdb_path):
        print(f"PDB file not found: {pdb_path}")
        sys.exit(1)

    native = pdbrust.parse_pdb_file(pdb_path)
    print(f"Loaded: {native}")
    print(f"Chains: {native.get_chain_ids()}")
    print(f"Atoms:  {native.num_atoms}")

    # ================================================================
    # 1. Self-comparison (perfect DockQ = 1.0)
    # ================================================================
    print("\n--- 1. Self-comparison (expect DockQ = 1.0) ---")
    result = native.dockq_to(native)
    print_result(result)

    # ================================================================
    # 2. Using explicit chain mapping
    # ================================================================
    print("\n--- 2. Explicit chain mapping ---")
    options = DockQOptions(
        chain_mapping=ChainMappingStrategy.explicit([("A", "A"), ("B", "B")])
    )
    result = native.dockq_to_with_options(native, options)
    print(f"  Chain mapping: {result.chain_mapping}")
    print_result(result)

    # ================================================================
    # 3. Inspecting per-interface results
    # ================================================================
    print("\n--- 3. Per-interface details ---")
    result = native.dockq_to(native)
    for i, iface in enumerate(result.interfaces):
        print(f"  Interface {i + 1}:")
        print(f"    Native chains: {iface.native_receptor_chain}-{iface.native_ligand_chain}")
        print(f"    Model chains:  {iface.model_receptor_chain}-{iface.model_ligand_chain}")
        print(f"    DockQ:  {iface.dockq:.4f}")
        print(f"    fnat:   {iface.fnat:.4f}")
        print(f"    fnonnat:{iface.fnonnat:.4f}")
        print(f"    iRMSD:  {iface.irmsd:.4f} A")
        print(f"    LRMSD:  {iface.lrmsd:.4f} A")
        print(f"    F1:     {iface.f1:.4f}")
        print(f"    Native contacts: {iface.num_native_contacts}")
        print(f"    Model contacts:  {iface.num_model_contacts}")
        print(f"    Clashes: {iface.num_clashes}")

    # ================================================================
    # 4. Quality classification
    # ================================================================
    print("\n--- 4. Quality classification ---")
    for score in [0.10, 0.30, 0.60, 0.90]:
        quality = DockQQuality.from_score(score)
        print(f"  DockQ = {score:.2f} -> {quality}")

    # ================================================================
    # 5. Indexing into result
    # ================================================================
    print("\n--- 5. DockQResult iteration ---")
    result = native.dockq_to(native)
    print(f"  Number of interfaces: {len(result)}")
    if len(result) > 0:
        first = result[0]
        print(f"  First interface: {first}")

    print("\n=== Demo complete ===")


def print_result(result):
    """Print a DockQ result summary."""
    print(f"  Total DockQ: {result.total_dockq:.4f} ({result.num_interfaces} interface(s))")
    for iface in result.interfaces:
        print(
            f"    {iface.native_receptor_chain}-{iface.native_ligand_chain}: "
            f"DockQ={iface.dockq:.4f} fnat={iface.fnat:.4f} "
            f"iRMSD={iface.irmsd:.2f} LRMSD={iface.lrmsd:.2f} "
            f"({iface.quality})"
        )


if __name__ == "__main__":
    main()
