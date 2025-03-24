//! PDBRust: A Rust library for parsing and analyzing Protein Data Bank (PDB) files
//! 
//! This library provides a robust and efficient way to work with PDB files, supporting various record
//! types including ATOM, SEQRES, CONECT, SSBOND, and more. It offers functionality for structural
//! analysis, sequence information retrieval, and connectivity analysis of molecular structures.
//!
//! For detailed documentation, examples, and best practices, see the [guide](guide) module.
//!
//! # Features
//!
//! - Parse PDB files with comprehensive error handling
//! - Support for multiple models in a single structure
//! - Chain and residue analysis
//! - Connectivity information through CONECT records
//! - Sequence information through SEQRES records
//! - Support for disulfide bonds through SSBOND records
//! - Remark handling for additional structural information
//!
//! # Example
//!
//! ```rust
//! use pdbrust::PdbStructure;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let structure = PdbStructure::from_file("example.pdb")?;
//!     
//!     // Get all chain IDs in the structure
//!     let chains = structure.get_chain_ids();
//!     
//!     // Get sequence for a specific chain
//!     if let Some(chain_id) = chains.first() {
//!         let sequence = structure.get_sequence(chain_id);
//!         println!("Sequence for chain {}: {:?}", chain_id, sequence);
//!     }
//!     
//!     Ok(())
//! }
//! ```

/// Documentation and usage guide
pub mod guide;

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use thiserror::Error;

/// Represents errors that can occur during PDB file parsing and analysis.
#[derive(Debug, Error)]
pub enum PdbError {
    /// Represents I/O errors that occur when reading PDB files.
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    
    /// Indicates that a record in the PDB file has an invalid format.
    #[error("Invalid record format: {0}")]
    InvalidRecord(String),
    
    /// Represents general parsing errors that occur during data extraction.
    #[error("Parse error: {0}")]
    ParseError(String),
}

/// Represents an atom record from a PDB file.
///
/// Contains all standard PDB ATOM record fields including position,
/// identification, and thermal factor information.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atom serial number.
    pub serial: i32,
    /// Atom name.
    pub name: String,
    /// Alternate location indicator (if any).
    pub alt_loc: Option<char>,
    /// Residue name.
    pub residue_name: String,
    /// Chain identifier.
    pub chain_id: String,
    /// Residue sequence number.
    pub residue_seq: i32,
    /// X coordinate in Angstroms.
    pub x: f64,
    /// Y coordinate in Angstroms.
    pub y: f64,
    /// Z coordinate in Angstroms.
    pub z: f64,
    /// Occupancy.
    pub occupancy: f64,
    /// Temperature factor.
    pub temp_factor: f64,
    /// Element symbol.
    pub element: String,
    /// Insertion code.
    pub ins_code: Option<char>, 
}

/// Represents a SEQRES record from a PDB file.
///
/// Contains sequence information for a specific chain in the structure.
#[derive(Debug, Clone)]
pub struct SeqRes {
    /// Serial number of the SEQRES record.
    pub serial: i32,
    /// Chain identifier.
    pub chain_id: String,
    /// Number of residues in the sequence.
    pub num_residues: i32,
    /// List of residue names in the sequence.
    pub residues: Vec<String>,
}

/// Represents a CONECT record from a PDB file.
///
/// Specifies connectivity between atoms in the structure.
#[derive(Debug, Clone)]
pub struct Conect {
    /// Serial number of the central atom.
    pub atom_serial: i32,
    /// Serial numbers of atoms bonded to the central atom.
    pub bonded_atoms: Vec<i32>,
}

/// Represents an SSBOND record from a PDB file.
///
/// Describes disulfide bonds between pairs of cysteine residues.
#[derive(Debug, Clone)]
pub struct SSBond {
    /// Serial number of the SSBOND record.
    pub serial: i32,
    /// Name of the first residue.
    pub residue1_name: String,
    /// Chain identifier for the first residue.
    pub chain1_id: String,
    /// Sequence number of the first residue.
    pub residue1_seq: i32,
    /// Name of the second residue.
    pub residue2_name: String,
    /// Chain identifier for the second residue.
    pub chain2_id: String,
    /// Sequence number of the second residue.
    pub residue2_seq: i32,
    /// Distance between the two sulfur atoms (if available).
    pub distance: Option<f64>,
}

/// Represents a REMARK record from a PDB file.
///
/// Contains additional information about the structure.
#[derive(Debug, Clone)]
pub struct Remark {
    /// Remark number identifying the type of information.
    pub number: i32,
    /// Content of the remark.
    pub content: String,
}

/// Represents a MODEL record from a PDB file.
///
/// Contains atoms and remarks specific to one model in a multi-model structure.
#[derive(Debug, Clone)]
pub struct Model {
    /// Serial number of the model.
    pub serial: i32,
    /// List of atoms in the model.
    pub atoms: Vec<Atom>,
    /// List of remarks specific to this model.
    pub remarks: Vec<Remark>,
}

/// Main structure representing a PDB file.
///
/// This structure contains all the information parsed from a PDB file,
/// including atomic coordinates, sequence information, connectivity,
/// and other structural features.
#[derive(Debug, Default)]
pub struct PdbStructure {
    /// List of all atoms in the structure.
    pub atoms: Vec<Atom>,
    /// Header information from the PDB file.
    pub header: Option<String>,
    /// Title of the structure.
    pub title: Option<String>,
    /// Sequence information from SEQRES records.
    pub seqres: Vec<SeqRes>,
    /// Connectivity information from CONECT records.
    pub connects: Vec<Conect>,
    /// Disulfide bond information from SSBOND records.
    pub ssbonds: Vec<SSBond>,
    /// Additional information from REMARK records.
    pub remarks: Vec<Remark>,
    /// List of models in case of multi-model structures.
    pub models: Vec<Model>,
    /// Currently active model number (if any).
    pub current_model: Option<i32>,
}


impl PdbStructure {
    /// Creates a new empty PDB structure.
    ///
    /// # Returns
    ///
    /// Returns a new `PdbStructure` instance with all fields initialized to their default empty values.
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            header: None,
            title: None,
            seqres: Vec::new(),
            connects: Vec::new(),
            ssbonds: Vec::new(),
            remarks: Vec::new(),
            models: Vec::new(),
            current_model: None,
        }
    }

    /// Reads and parses a PDB file from the specified path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the PDB file to read
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `PdbStructure` or a `PdbError`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// println!("Number of atoms: {}", structure.atoms.len());
    /// ```
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, PdbError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut structure = PdbStructure::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                let atom = Self::parse_atom_record(&line)?;
                
                // If we're in a model, add it to the current model
                // Otherwise, add it to the general atoms list
                if let Some(model_num) = structure.current_model {
                    // Find the model and add the atom to it
                    if let Some(model) = structure.models.iter_mut()
                        .find(|m| m.serial == model_num) {
                        model.atoms.push(atom.clone());
                    }
                }
                
                // Always add to the general atoms list for backward compatibility
                structure.atoms.push(atom);
            } else if line.starts_with("HEADER") {
                structure.header = Some(line[6..].trim().to_string());
            } else if line.starts_with("TITLE ") {
                structure.title = Some(line[6..].trim().to_string());
            } else if line.starts_with("MODEL ") {
                let model_num = Self::parse_model_record(&line)?;
                structure.current_model = Some(model_num);
                structure.models.push(Model {
                    serial: model_num,
                    atoms: Vec::new(),
                    remarks: Vec::new(),
                });
            } else if line.starts_with("ENDMDL") {
                structure.current_model = None;
            } else if line.starts_with("SEQRES") {
                let seqres = Self::parse_seqres_record(&line)?;
                structure.seqres.push(seqres);
            } else if line.starts_with("CONECT") {
                let conect = Self::parse_conect_record(&line)?;
                structure.connects.push(conect);
            } else if line.starts_with("SSBOND") {
                let ssbond = Self::parse_ssbond_record(&line)?;
                structure.ssbonds.push(ssbond);
            } else if line.starts_with("REMARK") {
                let remark = Self::parse_remark_record(&line)?;
                structure.remarks.push(remark);
            }
        }

        Ok(structure)
    }

    /// Parses an ATOM record from a PDB file line.
    ///
    /// This method extracts all relevant information from an ATOM record line
    /// according to the PDB file format specification.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the ATOM record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `Atom` or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if numeric values cannot be parsed.
    fn parse_atom_record(line: &str) -> Result<Atom, PdbError> {
        if line.len() < 66 {
            return Err(PdbError::InvalidRecord(
                "ATOM record line too short".to_string(),
            ));
        }

        let serial = line[6..11]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid atom serial number".to_string()))?;

        let name = line[12..16].trim().to_string();
        let alt_loc = {
            let c = line.chars().nth(16).unwrap_or(' ');
            if c == ' ' {
                None
            } else {
                Some(c)
            }
        };

        let residue_name = line[17..20].trim().to_string();
        let chain_id = line[21..22].trim().to_string();
        
        let residue_seq = line[22..26]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid residue sequence number".to_string()))?;

        let x = line[30..38]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid x coordinate".to_string()))?;
        let y = line[38..46]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid y coordinate".to_string()))?;
        let z = line[46..54]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid z coordinate".to_string()))?;

        let occupancy = line[54..60]
            .trim()
            .parse()
            .unwrap_or(1.0);
        let temp_factor = line[60..66]
            .trim()
            .parse()
            .unwrap_or(0.0);

        let element = if line.len() >= 78 {
            line[76..78].trim().to_string()
        } else {
            "".to_string()
        };

        let ins_code = {
            let c = if line.len() > 26 {
                line.chars().nth(26).unwrap_or(' ')
            } else {
                ' '
            };
            if c == ' ' {
                None
            } else {
                Some(c)
            }
        };

        Ok(Atom {
            serial,
            name,
            alt_loc,
            residue_name,
            chain_id,
            residue_seq,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element,
            ins_code,
        })
    }

    /// Parses a MODEL record from a PDB file line.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the MODEL record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed model serial number
    /// or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if the model number cannot be parsed.
    fn parse_model_record(line: &str) -> Result<i32, PdbError> {
        if line.len() < 14 {
            return Err(PdbError::InvalidRecord(
                "MODEL record line too short".to_string(),
            ));
        }

        line[10..14]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid model serial number".to_string()))
    }

    /// Parses a SEQRES record from a PDB file line.
    ///
    /// This method extracts sequence information including chain ID,
    /// residue count, and residue names from a SEQRES record.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the SEQRES record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `SeqRes` or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if numeric values cannot be parsed.
    fn parse_seqres_record(line: &str) -> Result<SeqRes, PdbError> {
        if line.len() < 70 {
            return Err(PdbError::InvalidRecord(
                "SEQRES record line too short".to_string(),
            ));
        }

        let serial = line[7..10]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid SEQRES serial number".to_string()))?;

        let chain_id = line[11..12].trim().to_string();
        
        let num_residues = line[13..17]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid number of residues".to_string()))?;

        let mut residues = Vec::new();
        let mut pos = 19;
        while pos + 3 <= line.len() && residues.len() < 13 {  // Max 13 residues per line
            let residue = line[pos..pos+3].trim();
            if !residue.is_empty() {
                residues.push(residue.to_string());
            }
            pos += 4;  // Move to next residue (3 chars + 1 space)
        }

        Ok(SeqRes {
            serial,
            chain_id,
            num_residues,
            residues,
        })
    }

    /// Parses a CONECT record from a PDB file line.
    ///
    /// This method extracts connectivity information between atoms,
    /// including the central atom and all atoms bonded to it.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the CONECT record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `Conect` or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if atom serial numbers cannot be parsed.
    fn parse_conect_record(line: &str) -> Result<Conect, PdbError> {
        if line.len() < 11 {
            return Err(PdbError::InvalidRecord(
                "CONECT record line too short".to_string(),
            ));
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            return Err(PdbError::InvalidRecord(
                "CONECT record must have at least one atom serial number".to_string(),
            ));
        }

        let atom_serial = parts[1]
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid atom serial number".to_string()))?;

        let bonded_atoms: Vec<i32> = parts[2..]
            .iter()
            .filter_map(|s| s.parse().ok())
            .collect();

        Ok(Conect {
            atom_serial,
            bonded_atoms,
        })
    }

    /// Parses an SSBOND record from a PDB file line.
    ///
    /// This method extracts information about disulfide bonds between
    /// cysteine residues, including residue names, chain IDs, and
    /// sequence numbers.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the SSBOND record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `SSBond` or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if numeric values cannot be parsed.
    fn parse_ssbond_record(line: &str) -> Result<SSBond, PdbError> {
        if line.len() < 70 {
            return Err(PdbError::InvalidRecord(
                "SSBOND record line too short".to_string(),
            ));
        }

        let serial = line[7..10]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid SSBOND serial number".to_string()))?;

        let residue1_name = line[11..14].trim().to_string();
        let chain1_id = line[15..16].trim().to_string();
        let residue1_seq = line[17..21]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid residue1 sequence number".to_string()))?;

        let residue2_name = line[25..28].trim().to_string();
        let chain2_id = line[29..30].trim().to_string();
        let residue2_seq = line[31..35]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid residue2 sequence number".to_string()))?;

        // Distance is in columns 73-78
        let distance = if line.len() >= 78 {
            Some(line[72..78]
                .trim()
                .parse()
                .unwrap_or(0.0))
        } else {
            Some(0.0)  // Default distance if not specified
        };

        Ok(SSBond {
            serial,
            residue1_name,
            chain1_id,
            residue1_seq,
            residue2_name,
            chain2_id,
            residue2_seq,
            distance,
        })
    }

    /// Parses a REMARK record from a PDB file line.
    ///
    /// This method extracts the remark number and content from a REMARK record.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the REMARK record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `Remark` or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if the remark number cannot be parsed.
    fn parse_remark_record(line: &str) -> Result<Remark, PdbError> {
        if line.len() < 11 {
            return Err(PdbError::InvalidRecord(
                "REMARK record line too short".to_string(),
            ));
        }

        // REMARK number is in columns 7-10
        let number = line[7..10]
            .trim()
            .parse()
            .unwrap_or(0); // Default to 0 if not a valid number

        // Content starts at column 11
        let content = line[11..].trim().to_string();

        Ok(Remark { number, content })
    }

    /// Returns a list of all unique chain identifiers in the structure.
    ///
    /// This method collects chain IDs from both ATOM records and SEQRES records
    /// to ensure complete coverage of all chains in the structure.
    ///
    /// # Returns
    ///
    /// Returns a `Vec<String>` containing unique chain identifiers.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let chains = structure.get_chain_ids();
    /// println!("Structure contains chains: {:?}", chains);
    /// ```
    pub fn get_chain_ids(&self) -> Vec<String> {
        let mut chain_ids: Vec<String> = self.atoms
            .iter()
            .map(|atom| atom.chain_id.clone())
            .collect();
        chain_ids.sort();
        chain_ids.dedup();
        chain_ids
    }

    /// Gets all residues for a specific chain.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - The chain identifier to get residues for
    ///
    /// # Returns
    ///
    /// Returns a `Vec` of tuples containing residue sequence numbers and names
    /// for the specified chain.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let residues = structure.get_residues_for_chain("A");
    /// for (seq, name) in residues {
    ///     println!("Residue {} at position {}", name, seq);
    /// }
    /// ```
    pub fn get_residues_for_chain(&self, chain_id: &str) -> Vec<(i32, String)> {
        let mut residues: Vec<(i32, String)> = self.atoms
            .iter()
            .filter(|atom| atom.chain_id == chain_id)
            .map(|atom| (atom.residue_seq, atom.residue_name.clone()))
            .collect();
        residues.sort();
        residues.dedup();
        residues
    }

    /// Gets the amino acid sequence for a specific chain from SEQRES records.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - The chain identifier to get the sequence for
    ///
    /// # Returns
    ///
    /// Returns a `Vec<String>` containing the sequence of residue names
    /// for the specified chain.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let sequence = structure.get_sequence("A");
    /// println!("Chain A sequence: {:?}", sequence);
    /// ```
    pub fn get_sequence(&self, chain_id: &str) -> Vec<String> {
        self.seqres
            .iter()
            .filter(|sr| sr.chain_id == chain_id)
            .flat_map(|sr| sr.residues.clone())
            .collect()
    }

    /// Gets all remarks with a specific remark number.
    ///
    /// # Arguments
    ///
    /// * `number` - The remark number to filter by
    ///
    /// # Returns
    ///
    /// Returns a `Vec` of references to `Remark` structs matching the specified number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let resolution_remarks = structure.get_remarks_by_number(2);
    /// for remark in resolution_remarks {
    ///     println!("Resolution info: {}", remark.content);
    /// }
    /// ```
    pub fn get_remarks_by_number(&self, number: i32) -> Vec<&Remark> {
        self.remarks
            .iter()
            .filter(|r| r.number == number)
            .collect()
    }

    /// Gets all atoms that are connected to a specific atom.
    ///
    /// # Arguments
    ///
    /// * `atom_serial` - The serial number of the atom to find connections for
    ///
    /// # Returns
    ///
    /// Returns a `Vec` of references to `Atom` structs that are connected
    /// to the specified atom according to CONECT records.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let connected_atoms = structure.get_connected_atoms(1);
    /// for atom in connected_atoms {
    ///     println!("Connected to atom {} at ({}, {}, {})",
    ///              atom.serial, atom.x, atom.y, atom.z);
    /// }
    /// ```
    pub fn get_connected_atoms(&self, atom_serial: i32) -> Vec<&Atom> {
        // First find all atoms that are connected to this atom
        let connected_serials: Vec<i32> = self.connects
            .iter()
            .filter(|c| c.atom_serial == atom_serial)
            .flat_map(|c| c.bonded_atoms.clone())
            .collect();

        // Then get the actual atom references
        self.atoms
            .iter()
            .filter(|a| connected_serials.contains(&a.serial))
            .collect()
    }
    /// Translate all atoms in the structure by a given amount.     
    ///
    /// # Arguments
    ///
    /// * `dx` - The amount to translate along the x-axis
    /// * `dy` - The amount to translate along the y-axis
    /// * `dz` - The amount to translate along the z-axis
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mut structure = PdbStructure::new();
    /// structure.translate(1.0, 2.0, 3.0);
    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        for atom in &mut self.atoms {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
        }
    }

    /// Writes the structure to a PDB file.
    ///
    /// # Arguments
    ///
    /// * `path` - Path where the PDB file should be written
    ///
    /// # Returns
    ///
    /// Returns a `Result` indicating success or containing a `PdbError`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mut structure = PdbStructure::new();
    /// // Add atoms and other records...
    /// structure.to_file("output.pdb")?;
    /// ```
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), PdbError> {
        let mut file = File::create(path)?;

        // Write header if present
        if let Some(header) = &self.header {
            writeln!(file, "HEADER    {:<72}", header)?;
        }

        // Write title if present
        if let Some(title) = &self.title {
            writeln!(file, "TITLE     {:<70}", title)?;
        }

        // Write remarks
        for remark in &self.remarks {
            writeln!(file, "{}", Self::format_remark_record(remark))?;
        }

        // Write SEQRES records
        for seqres in &self.seqres {
            writeln!(file, "{}", Self::format_seqres_record(seqres))?;
        }

        // Write SSBOND records
        for ssbond in &self.ssbonds {
            writeln!(file, "{}", Self::format_ssbond_record(ssbond))?;
        }

        // If there are models, write them with MODEL/ENDMDL records
        if !self.models.is_empty() {
            for model in &self.models {
                writeln!(file, "MODEL     {:>4}", model.serial)?;
                for atom in &model.atoms {
                    writeln!(file, "{}", Self::format_atom_record(atom))?;
                }
                writeln!(file, "ENDMDL")?;
            }
        } else {
            // Otherwise write atoms directly
            for atom in &self.atoms {
                writeln!(file, "{}", Self::format_atom_record(atom))?;
            }
        }

        // Write CONECT records
        for conect in &self.connects {
            writeln!(file, "{}", Self::format_conect_record(conect))?;
        }

        // Write END record
        writeln!(file, "END")?;

        Ok(())
    }

    /// Formats an ATOM record according to the PDB format specification.
    fn format_atom_record(atom: &Atom) -> String {
        format!(
            "{:<6}{:>5} {:<4}{:3} {}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            if atom.name.starts_with("H") { "HETATM" } else { "ATOM" },
            atom.serial,
            atom.name,
            atom.residue_name,
            atom.chain_id,
            atom.residue_seq,
            atom.x,
            atom.y,
            atom.z,
            atom.occupancy,
            atom.temp_factor,
            atom.element
        )
    }

    /// Formats a SEQRES record according to the PDB format specification.
    fn format_seqres_record(seqres: &SeqRes) -> String {
        let mut residues_str = String::new();
        for residue in &seqres.residues {
            residues_str.push_str(&format!("{:<4}", residue));
        }
        format!(
            "SEQRES  {:>3} {} {:>4}  {}",
            seqres.serial,
            seqres.chain_id,
            seqres.num_residues,
            residues_str
        )
    }

    /// Formats a CONECT record according to the PDB format specification.
    fn format_conect_record(conect: &Conect) -> String {
        let mut result = format!("CONECT{:>5}", conect.atom_serial);
        for bonded in &conect.bonded_atoms {
            result.push_str(&format!("{:>5}", bonded));
        }
        result
    }

    /// Formats an SSBOND record according to the PDB format specification.
    fn format_ssbond_record(ssbond: &SSBond) -> String {
        format!(
            "SSBOND {:>3} {} {} {:>4}    {} {} {:>4}                       {:>6.2}",
            ssbond.serial,
            ssbond.residue1_name,
            ssbond.chain1_id,
            ssbond.residue1_seq,
            ssbond.residue2_name,
            ssbond.chain2_id,
            ssbond.residue2_seq,
            ssbond.distance.unwrap_or(0.0)
        )
    }

    /// Formats a REMARK record according to the PDB format specification.
    fn format_remark_record(remark: &Remark) -> String {
        format!("REMARK {:>3} {}", remark.number, remark.content)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_atom_record() {
        let atom_line = "ATOM      1  N   ASP A  30      31.904  -0.904  -0.904  1.00 20.00           N  ";
        let atom = PdbStructure::parse_atom_record(atom_line).unwrap();
        
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "N");
        assert_eq!(atom.residue_name, "ASP");
        assert_eq!(atom.chain_id, "A");
        assert_eq!(atom.residue_seq, 30);
        assert!((atom.x - 31.904).abs() < 1e-6);
        assert!((atom.y - (-0.904)).abs() < 1e-6);
        assert!((atom.z - (-0.904)).abs() < 1e-6);
        assert!((atom.occupancy - 1.0).abs() < 1e-6);
        assert!((atom.temp_factor - 20.0).abs() < 1e-6);
        assert_eq!(atom.element, "N");
    }

    #[test]
    fn test_parse_seqres_record() {
        let seqres_line = "SEQRES   1 A  20  MET ALA CYS PRO GLY VAL GLY PRO GLU CYS ARG GLU MET SER PRO";
        let seqres = PdbStructure::parse_seqres_record(seqres_line).unwrap();
        
        assert_eq!(seqres.serial, 1);
        assert_eq!(seqres.chain_id, "A");
        assert_eq!(seqres.num_residues, 20);
        let expected_residues = vec!["MET", "ALA", "CYS", "PRO", "GLY", "VAL", "GLY", "PRO", "GLU", "CYS", "ARG", "GLU", "MET", "SER", "PRO"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        assert_eq!(seqres.residues, expected_residues);
    }

    #[test]
    fn test_parse_conect_record() {
        let conect_line = "CONECT 1179 1174 1177 1180";
        let conect = PdbStructure::parse_conect_record(conect_line).unwrap();
        
        assert_eq!(conect.atom_serial, 1179);
        let expected_bonds = vec![1174, 1177, 1180];
        assert_eq!(conect.bonded_atoms, expected_bonds);
    }

    #[test]
    fn test_parse_ssbond_record() {
        let ssbond_line = "SSBOND   1 CYS A   85    CYS A  101                               2.03";
        let ssbond = PdbStructure::parse_ssbond_record(ssbond_line).unwrap();
        
        assert_eq!(ssbond.serial, 1);
        assert_eq!(ssbond.residue1_name, "CYS");
        assert_eq!(ssbond.chain1_id, "A");
        assert_eq!(ssbond.residue1_seq, 85);
        assert_eq!(ssbond.residue2_name, "CYS");
        assert_eq!(ssbond.chain2_id, "A");
        assert_eq!(ssbond.residue2_seq, 101);
        assert!(ssbond.distance.is_some());
        assert!((ssbond.distance.unwrap() - 2.03).abs() < 1e-6);
    }

    #[test]
    fn test_new_structure() {
        let structure = PdbStructure::new();
        assert!(structure.atoms.is_empty());
        assert!(structure.header.is_none());
        assert!(structure.title.is_none());
        assert!(structure.seqres.is_empty());
        assert!(structure.connects.is_empty());
        assert!(structure.ssbonds.is_empty());
        assert!(structure.remarks.is_empty());
        assert!(structure.models.is_empty());
        assert!(structure.current_model.is_none());
    }

    #[test]
    fn test_get_chain_ids() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "N".to_string(),
            ins_code: None,
        });
        structure.atoms.push(Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "B".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        });

        let chain_ids = structure.get_chain_ids();
        assert_eq!(chain_ids, vec!["A", "B"]);
    }

    #[test]
    fn test_get_residues_for_chain() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "N".to_string(),
            ins_code: None,
        });
        structure.atoms.push(Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        });

        let residues = structure.get_residues_for_chain("A");
        assert_eq!(residues, vec![(1, "ASP".to_string())]);
    }

    #[test]
    fn test_get_sequence() {
        let mut structure = PdbStructure::new();
        structure.seqres.push(SeqRes {
            serial: 1,
            chain_id: "A".to_string(),
            num_residues: 3,
            residues: vec!["ALA".to_string(), "GLY".to_string(), "VAL".to_string()],
        });

        let sequence = structure.get_sequence("A");
        assert_eq!(sequence, vec!["ALA", "GLY", "VAL"]);
    }

    #[test]
    fn test_get_connected_atoms() {
        let mut structure = PdbStructure::new();
        
        // Add some atoms
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "N".to_string(),
            ins_code: None,
        });
        structure.atoms.push(Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        });

        // Add connectivity
        structure.connects.push(Conect {
            atom_serial: 1,
            bonded_atoms: vec![2],
        });

        let connected = structure.get_connected_atoms(1);
        assert_eq!(connected.len(), 1);
        assert_eq!(connected[0].serial, 2);
        assert_eq!(connected[0].name, "CA");
    }

    #[test]
    fn test_write_and_read_pdb() -> Result<(), Box<dyn std::error::Error>> {
        // Create a test structure
        let mut structure = PdbStructure::new();
        structure.header = Some("TEST STRUCTURE".to_string());
        structure.title = Some("TEST TITLE".to_string());
        
        // Add an atom
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 2.0,
            z: 3.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "N".to_string(),
            ins_code: None,
        });

        // Add a SEQRES record
        structure.seqres.push(SeqRes {
            serial: 1,
            chain_id: "A".to_string(),
            num_residues: 1,
            residues: vec!["ASP".to_string()],
        });

        // Add a CONECT record
        structure.connects.push(Conect {
            atom_serial: 1,
            bonded_atoms: vec![2],
        });

        // Create a temporary file
        let temp_file = NamedTempFile::new()?;
        let temp_path = temp_file.path().to_owned();

        // Write the structure
        structure.to_file(&temp_path)?;

        // Read the file contents
        let mut contents = String::new();
        File::open(&temp_path)?.read_to_string(&mut contents)?;

        // Verify the contents
        assert!(contents.contains("HEADER    TEST STRUCTURE"));
        assert!(contents.contains("TITLE     TEST TITLE"));
        assert!(contents.contains("ATOM      1 N   ASP A   1       1.000   2.000   3.000  1.00 20.00           N"));
        assert!(contents.contains("SEQRES    1 A    1  ASP"));
        assert!(contents.contains("CONECT    1    2"));
        assert!(contents.contains("END"));

        // Read the structure back
        let read_structure = PdbStructure::from_file(&temp_path)?;

        // Verify the structure
        assert_eq!(read_structure.header, Some("TEST STRUCTURE".to_string()));
        assert_eq!(read_structure.title, Some("TEST TITLE".to_string()));
        assert_eq!(read_structure.atoms.len(), 1);
        assert_eq!(read_structure.atoms[0].serial, 1);
        assert_eq!(read_structure.atoms[0].name, "N");
        assert_eq!(read_structure.seqres.len(), 1);
        assert_eq!(read_structure.connects.len(), 1);

        Ok(())
    }

    #[test]
    fn test_format_atom_record() {
        let atom = Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 2.0,
            z: 3.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "N".to_string(),
            ins_code: None,
        };

        let formatted = PdbStructure::format_atom_record(&atom);
        assert_eq!(
            formatted,
            "ATOM      1 N   ASP A   1       1.000   2.000   3.000  1.00 20.00           N"
        );
    }

    #[test]
    fn test_format_hetatm_record() {
        let atom = Atom {
            serial: 1,
            name: "HOH".to_string(),
            alt_loc: None,
            residue_name: "WAT".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 2.0,
            z: 3.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "O".to_string(),
            ins_code: None,
        };

        let formatted = PdbStructure::format_atom_record(&atom);
        assert_eq!(
            formatted,
            "HETATM    1 HOH WAT A   1       1.000   2.000   3.000  1.00 20.00           O"
        );
    }
}
