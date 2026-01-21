"""
PDBRust: High-performance PDB/mmCIF parsing and analysis library.

This package provides Python bindings for the PDBRust Rust library,
offering 40-260x speedups over pure Python alternatives.

Basic usage:
    >>> import pdbrust
    >>> structure = pdbrust.parse_pdb_file("protein.pdb")
    >>> print(f"Loaded {structure.num_atoms} atoms")

    >>> # Filter operations
    >>> cleaned = structure.remove_ligands().keep_only_chain("A")

    >>> # Structural analysis
    >>> rg = structure.radius_of_gyration()
    >>> desc = structure.structure_descriptors()

    >>> # RCSB integration
    >>> from pdbrust import SearchQuery, rcsb_search, download_structure, FileFormat
    >>> query = SearchQuery().with_text("kinase").with_resolution_max(2.0)
    >>> results = rcsb_search(query, 10)
    >>> structure = download_structure("1UBQ", FileFormat.pdb())
"""

from pdbrust._pdbrust import (
    # Core types
    PdbStructure,
    Atom,
    SSBond,
    SeqRes,
    Conect,
    Remark,
    Model,
    # Parsing functions
    parse_pdb_file,
    parse_mmcif_file,
    parse_structure_file,
    parse_pdb_string,
    parse_mmcif_string,
    # Writing functions
    write_pdb_file,
    write_mmcif_file,
    write_mmcif_string,
)

# Optional: gzip support
try:
    from pdbrust._pdbrust import (
        parse_gzip_pdb_file,
        parse_gzip_mmcif_file,
        parse_gzip_structure_file,
        write_gzip_mmcif_file,
    )
except ImportError:
    pass

# Optional: descriptors
try:
    from pdbrust._pdbrust import StructureDescriptors, ResidueBFactor
except ImportError:
    pass

# Optional: quality
try:
    from pdbrust._pdbrust import QualityReport
except ImportError:
    pass

# Optional: summary
try:
    from pdbrust._pdbrust import StructureSummary
except ImportError:
    pass

# Optional: RCSB
try:
    from pdbrust._pdbrust import (
        FileFormat,
        ExperimentalMethod,
        PolymerType,
        SearchQuery,
        SearchResult,
        rcsb_search,
        download_structure,
        download_pdb_string,
        download_to_file,
    )
except ImportError:
    pass

# Optional: RCSB Async
try:
    from pdbrust._pdbrust import (
        AsyncDownloadOptions,
        DownloadResult,
        download_multiple,
    )
except ImportError:
    pass

# Optional: geometry (RMSD, alignment)
try:
    from pdbrust._pdbrust import (
        AtomSelection,
        AlignmentResult,
        PerResidueRmsd,
    )
except ImportError:
    pass

# Optional: DSSP secondary structure
try:
    from pdbrust._pdbrust import (
        SecondaryStructure,
        ResidueSSAssignment,
        SecondaryStructureAssignment,
    )
except ImportError:
    pass

__version__ = "0.6.0"
__all__ = [
    # Core types
    "PdbStructure",
    "Atom",
    "SSBond",
    "SeqRes",
    "Conect",
    "Remark",
    "Model",
    # Parsing
    "parse_pdb_file",
    "parse_mmcif_file",
    "parse_structure_file",
    "parse_pdb_string",
    "parse_mmcif_string",
    # Writing
    "write_pdb_file",
    "write_mmcif_file",
    "write_mmcif_string",
    # Gzip (optional)
    "parse_gzip_pdb_file",
    "parse_gzip_mmcif_file",
    "parse_gzip_structure_file",
    "write_gzip_mmcif_file",
    # Descriptors (optional)
    "StructureDescriptors",
    "ResidueBFactor",
    # Quality (optional)
    "QualityReport",
    # Summary (optional)
    "StructureSummary",
    # RCSB (optional)
    "FileFormat",
    "ExperimentalMethod",
    "PolymerType",
    "SearchQuery",
    "SearchResult",
    "rcsb_search",
    "download_structure",
    "download_pdb_string",
    "download_to_file",
    # RCSB Async (optional)
    "AsyncDownloadOptions",
    "DownloadResult",
    "download_multiple",
    # Geometry (optional)
    "AtomSelection",
    "AlignmentResult",
    "PerResidueRmsd",
    # DSSP (optional)
    "SecondaryStructure",
    "ResidueSSAssignment",
    "SecondaryStructureAssignment",
]
