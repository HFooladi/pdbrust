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

/// Convert a 2D Vec<Vec<f64>> to a numpy array (N x N)
///
/// Used for distance matrices.
#[cfg(feature = "descriptors")]
pub fn vec2d_to_array2<'py>(py: Python<'py>, data: &[Vec<f64>]) -> Bound<'py, PyArray2<f64>> {
    let n = data.len();
    if n == 0 {
        return PyArray2::zeros_bound(py, [0, 0], false);
    }

    let m = data[0].len();
    let array = PyArray2::zeros_bound(py, [n, m], false);

    unsafe {
        let ptr: *mut f64 = array.as_raw_array_mut().as_mut_ptr();
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                *ptr.add(i * m + j) = val;
            }
        }
    }

    array
}

/// Convert a 2D Vec<Vec<bool>> to a numpy boolean array (N x N)
///
/// Used for contact maps.
#[cfg(feature = "descriptors")]
pub fn vec2d_bool_to_array2<'py>(py: Python<'py>, data: &[Vec<bool>]) -> Bound<'py, PyArray2<bool>> {
    let n = data.len();
    if n == 0 {
        return PyArray2::zeros_bound(py, [0, 0], false);
    }

    let m = data[0].len();
    let array = PyArray2::zeros_bound(py, [n, m], false);

    unsafe {
        let ptr: *mut bool = array.as_raw_array_mut().as_mut_ptr();
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                *ptr.add(i * m + j) = val;
            }
        }
    }

    array
}
