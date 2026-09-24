//! R bindings for doublet-rs.
//!
//! These `rs_*` functions are internal: the exported, documented R API in
//! `R/doubletrs.R` validates arguments and calls them. Inputs are borrowed
//! straight from R memory (dgCMatrix slots, numeric matrices) without
//! copying.

use doublet_rs::knn::{HnswParams, KnnMethod};
use extendr_api::prelude::*;
use extendr_api::{Error, Result};
use ndarray::{ArrayView2, ShapeBuilder};
use sprs::CsMatViewI;

fn fail(msg: impl std::fmt::Display) -> Error {
    Error::Other(msg.to_string())
}

fn to_usize(value: i32, name: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| fail(format!("`{name}` must be non-negative, got {value}")))
}

fn knn_method(
    method: &str,
    max_connections: i32,
    ef_construction: i32,
    ef_search: i32,
) -> Result<KnnMethod> {
    match method {
        "exact" => Ok(KnnMethod::Exact),
        "hnsw" => Ok(KnnMethod::Hnsw(HnswParams {
            max_connections: to_usize(max_connections, "max_connections")?,
            ef_construction: to_usize(ef_construction, "ef_construction")?,
            ef_search: to_usize(ef_search, "ef_search")?,
        })),
        other => Err(fail(format!(
            "unknown method `{other}`, use \"exact\" or \"hnsw\""
        ))),
    }
}

fn embedding(pcs: &RMatrix<f64>) -> Result<ArrayView2<'_, f64>> {
    // R matrices are column-major.
    ArrayView2::from_shape((pcs.nrows(), pcs.ncols()).f(), pcs.data()).map_err(fail)
}

/// Simulate doublets from dgCMatrix slots.
/// Returns the doublet matrix slots and 1-based parent cell indices.
/// @noRd
#[extendr]
fn rs_simulate_doublets(
    i: Robj,
    p: Robj,
    x: Robj,
    n_genes: i32,
    n_cells: i32,
    n_doublets: i32,
    seed: i32,
) -> Result<List> {
    let i = i
        .as_integer_slice()
        .ok_or_else(|| fail("`i` slot must be integer"))?;
    let p = p
        .as_integer_slice()
        .ok_or_else(|| fail("`p` slot must be integer"))?;
    let x = x
        .as_real_slice()
        .ok_or_else(|| fail("`x` slot must be double"))?;
    let shape = (to_usize(n_genes, "n_genes")?, to_usize(n_cells, "n_cells")?);
    let expr = CsMatViewI::try_new_csc(shape, p, i, x)
        .map_err(|(.., e)| fail(format!("invalid dgCMatrix: {e}")))?;

    let doublets = doublet_rs::simulate_doublets(
        expr,
        to_usize(n_doublets, "n_doublets")?,
        seed as i64 as u64,
    )
    .map_err(fail)?;

    let (cell_a, cell_b): (Vec<i32>, Vec<i32>) = doublets
        .pairs
        .iter()
        .map(|&(a, b)| (a as i32 + 1, b as i32 + 1))
        .unzip();
    let (indptr, indices, data) = doublets.matrix.into_raw_storage();
    Ok(list!(
        i = indices,
        p = indptr,
        x = data,
        cell_a = cell_a,
        cell_b = cell_b
    ))
}

/// k nearest neighbours of the first `n_real` rows, 1-based.
/// @noRd
#[extendr]
fn rs_find_neighbors(
    pcs: RMatrix<f64>,
    n_real: i32,
    k: i32,
    method: &str,
    max_connections: i32,
    ef_construction: i32,
    ef_search: i32,
) -> Result<RMatrix<i32>> {
    let method = knn_method(method, max_connections, ef_construction, ef_search)?;
    let (n_real, k) = (to_usize(n_real, "n_real")?, to_usize(k, "k")?);
    let neighbors =
        doublet_rs::find_neighbors(embedding(&pcs)?, n_real, k, method).map_err(fail)?;

    // Approximate search can in principle return fewer than k; pad with NA.
    Ok(RMatrix::new_matrix(n_real, k, |r, c| {
        neighbors[r].get(c).map_or(i32::MIN, |&j| j as i32 + 1)
    }))
}

/// pANN for the first `n_real` rows.
/// @noRd
#[extendr]
fn rs_compute_pann(
    pcs: RMatrix<f64>,
    n_real: i32,
    k: i32,
    method: &str,
    max_connections: i32,
    ef_construction: i32,
    ef_search: i32,
) -> Result<Vec<f64>> {
    let method = knn_method(method, max_connections, ef_construction, ef_search)?;
    doublet_rs::pann_from_embedding(
        embedding(&pcs)?,
        to_usize(n_real, "n_real")?,
        to_usize(k, "k")?,
        method,
    )
    .map_err(fail)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod doubletrs;
    fn rs_simulate_doublets;
    fn rs_find_neighbors;
    fn rs_compute_pann;
}
