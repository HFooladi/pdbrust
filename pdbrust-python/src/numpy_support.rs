//! Numpy array support for efficient coordinate access

use numpy::{PyArray2, PyArrayMethods};
use pdbrust::PdbStructure;
use pyo3::prelude::*;

/// Get all atom coordinates as a numpy array (N x 3)
pub fn get_coords_array<'py>(
    py: Python<'py>,
    structure: &PdbStructure,
) -> Bound<'py, PyArray2<f64>> {
    let atoms = &structure.atoms;
    let n = atoms.len();

    // Create array and fill it
    let array = PyArray2::zeros_bound(py, [n, 3], false);

    unsafe {
        let ptr: *mut f64 = array.as_raw_array_mut().as_mut_ptr();
        for (i, atom) in atoms.iter().enumerate() {
            *ptr.add(i * 3) = atom.x;
            *ptr.add(i * 3 + 1) = atom.y;
            *ptr.add(i * 3 + 2) = atom.z;
        }
    }

    array
}

/// Get CA atom coordinates as a numpy array (N x 3)
#[cfg(feature = "filter")]
pub fn get_ca_coords_array<'py>(
    py: Python<'py>,
    structure: &PdbStructure,
    chain_id: Option<&str>,
) -> Bound<'py, PyArray2<f64>> {
    // get_ca_coords takes Option<&str> directly (inherent method on PdbStructure)
    let coords: Vec<(f64, f64, f64)> = structure.get_ca_coords(chain_id);

    let n = coords.len();
    let array = PyArray2::zeros_bound(py, [n, 3], false);

    unsafe {
        let ptr: *mut f64 = array.as_raw_array_mut().as_mut_ptr();
        for (i, (x, y, z)) in coords.iter().enumerate() {
            *ptr.add(i * 3) = *x;
            *ptr.add(i * 3 + 1) = *y;
            *ptr.add(i * 3 + 2) = *z;
        }
    }

    array
}

/// Get backbone atom coordinates as a numpy array (N x 3)
#[cfg(feature = "filter")]
pub fn get_backbone_coords_array<'py>(
    py: Python<'py>,
    structure: &PdbStructure,
) -> Bound<'py, PyArray2<f64>> {
    // keep_only_backbone is an inherent method on PdbStructure
    let backbone = structure.clone().keep_only_backbone();
    get_coords_array(py, &backbone)
}
