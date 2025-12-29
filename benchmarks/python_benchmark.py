#!/usr/bin/env python3
"""
Benchmark script for Python libraryPDB.

This script benchmarks parsing and analysis operations to compare
with the Rust pdbrust library performance.
"""

import sys
import time
import os
import statistics

# Add the libraryPDB to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'libraryPDB'))

from libraryPDB.PDBparser import parse_atoms, get_chains, get_residues, get_ca_coords
from libraryPDB.PDBdescriptors import (
    num_residues, num_atoms, aa_composition, glycine_ratio, hydrophobic_ratio,
    radius_of_gyration, max_ca_distance, missing_residue_ratio,
    secondary_structure_ratio, compactness_index, ca_density
)

# Test files
TEST_FILES = {
    "1UBQ": "examples/pdb_files/1UBQ.pdb",
    "8HM2": "examples/pdb_files/8HM2.pdb",
}

def benchmark(func, *args, iterations=100):
    """Run a function multiple times and return timing statistics."""
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        result = func(*args)
        end = time.perf_counter()
        times.append((end - start) * 1000)  # Convert to milliseconds

    return {
        'mean': statistics.mean(times),
        'std': statistics.stdev(times) if len(times) > 1 else 0,
        'min': min(times),
        'max': max(times),
        'result': result
    }

def run_benchmarks():
    """Run all benchmarks and print results."""
    print("=" * 70)
    print("Python libraryPDB Benchmark")
    print("=" * 70)
    print()

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    for name, rel_path in TEST_FILES.items():
        pdb_file = os.path.join(base_dir, rel_path)
        if not os.path.exists(pdb_file):
            print(f"Skipping {name}: file not found at {pdb_file}")
            continue

        print(f"Benchmarking: {name}")
        print("-" * 50)

        # Parsing benchmarks
        print("\n[Parsing Operations]")

        result = benchmark(parse_atoms, pdb_file, iterations=100)
        print(f"  parse_atoms:        {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({len(result['result'])} atoms)")

        result = benchmark(get_chains, pdb_file, iterations=100)
        print(f"  get_chains:         {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']})")

        result = benchmark(get_residues, pdb_file, iterations=100)
        print(f"  get_residues:       {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({len(result['result'])} residues)")

        result = benchmark(get_ca_coords, pdb_file, iterations=100)
        print(f"  get_ca_coords:      {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({len(result['result'])} CA atoms)")

        # Descriptor benchmarks
        print("\n[Descriptor Operations]")

        result = benchmark(num_atoms, pdb_file, iterations=100)
        print(f"  num_atoms:          {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']})")

        result = benchmark(num_residues, pdb_file, iterations=100)
        print(f"  num_residues:       {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']})")

        result = benchmark(aa_composition, pdb_file, iterations=100)
        print(f"  aa_composition:     {result['mean']:8.3f} ms ± {result['std']:.3f} ms")

        result = benchmark(glycine_ratio, pdb_file, iterations=100)
        print(f"  glycine_ratio:      {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f})")

        result = benchmark(hydrophobic_ratio, pdb_file, iterations=100)
        print(f"  hydrophobic_ratio:  {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f})")

        result = benchmark(radius_of_gyration, pdb_file, iterations=100)
        print(f"  radius_of_gyration: {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f} Å)")

        result = benchmark(max_ca_distance, pdb_file, iterations=100)
        print(f"  max_ca_distance:    {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f} Å)")

        result = benchmark(missing_residue_ratio, pdb_file, iterations=100)
        print(f"  missing_res_ratio:  {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f})")

        result = benchmark(secondary_structure_ratio, pdb_file, iterations=100)
        print(f"  ss_ratio:           {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f})")

        result = benchmark(compactness_index, pdb_file, iterations=100)
        print(f"  compactness_index:  {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.4f})")

        result = benchmark(ca_density, pdb_file, iterations=100)
        print(f"  ca_density:         {result['mean']:8.3f} ms ± {result['std']:.3f} ms  ({result['result']:.6f})")

        print()

    print("=" * 70)
    print("Benchmark complete.")
    print("=" * 70)

if __name__ == "__main__":
    run_benchmarks()
