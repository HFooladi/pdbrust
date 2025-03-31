/// Bond types between two atoms
#[derive(Debug, Clone)]
pub enum Bond {
    /// A covalent bond
    Covalent,
    /// A disulfide bond S-S
    Disulfide,
    /// A hydrogen bond H-H
    Hydrogen,
    /// ?
    MetalCoordination,
    /// ?
    MisMatchedBasePairs,
    /// ?
    SaltBridge,
    /// ?
    CovalentModificationResidue,
    /// ?
    CovalentModificationNucleotideBase,
    /// ?
    CovalentModificationNucleotideSugar,
    /// ?
    CovalentModificationNucleotidePhosphate,
}

impl Bond {
    pub fn is_covalent(&self) -> bool {
        matches!(self, Bond::Covalent)
    }
}