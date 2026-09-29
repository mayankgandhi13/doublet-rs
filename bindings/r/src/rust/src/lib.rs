//! R bindings for doublet-rs.
//!
//! These `rs_*` functions are internal: the exported, documented R API in
//! `R/` validates arguments and calls them. Count matrices are borrowed
//! straight from R memory (dgCMatrix slots) without copying; PC matrices are
//! transposed to row-major once, which is small next to the counts.
//!
//! Computation runs on a thread pool of `threads` workers, and R objects
//! are only created afterwards, on R's own thread.

use doublet_rs::matrix::col_major_to_row_major;
use doublet_rs::{with_threads, CscView, Embedding, KnnMethod};
use extendr_api::prelude::*;
use extendr_api::{Error, Result};

fn fail(msg: impl std::fmt::Display) -> Error {
    Error::Other(msg.to_string())
}

fn to_usize(value: i32, name: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| fail(format!("`{name}` must be non-negative, got {value}")))
}

/// Row-major copy of an R numeric matrix, as (data, n_points, dims).
fn row_major(pcs: &RMatrix<f64>) -> (Vec<f64>, usize, usize) {
    let (n, d) = (pcs.nrows(), pcs.ncols());
    (col_major_to_row_major(pcs.data(), n, d), n, d)
}

/// Build doublets from dgCMatrix slots and 1-based parent cell indices
/// (drawn in R, so `set.seed()` controls them). Returns the doublet matrix
/// slots.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_build_doublets(
    i: Robj,
    p: Robj,
    x: Robj,
    n_genes: i32,
    n_cells: i32,
    cell_a: Vec<i32>,
    cell_b: Vec<i32>,
    threads: i32,
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
    let expr = CscView::new(
        to_usize(n_genes, "n_genes")?,
        to_usize(n_cells, "n_cells")?,
        p,
        i,
        x,
    )
    .map_err(fail)?;
    if cell_a.len() != cell_b.len() {
        return Err(fail("`cell_a` and `cell_b` must have the same length"));
    }
    let pairs: Vec<(usize, usize)> = cell_a
        .iter()
        .zip(&cell_b)
        .map(|(&a, &b)| Ok((to_usize(a - 1, "cell_a")?, to_usize(b - 1, "cell_b")?)))
        .collect::<Result<_>>()?;

    let matrix = with_threads(to_usize(threads, "threads")?, || {
        doublet_rs::simulate::build_doublets(expr, &pairs)
    })
    .and_then(|r| r)
    .map_err(fail)?;

    let (indptr, indices, data) = matrix.into_parts();
    Ok(list!(i = indices, p = indptr, x = data))
}

/// k nearest neighbours of the first `n_real` rows, 1-based.
/// @noRd
#[extendr]
fn rs_find_neighbors(pcs: RMatrix<f64>, n_real: i32, k: i32, threads: i32) -> Result<RMatrix<i32>> {
    let (data, n_points, dims) = row_major(&pcs);
    let (n_real, k) = (to_usize(n_real, "n_real")?, to_usize(k, "k")?);
    let neighbors = with_threads(to_usize(threads, "threads")?, || {
        let emb = Embedding::new(&data, n_points, dims)?;
        doublet_rs::find_neighbors(emb, n_real, k, KnnMethod::Exact)
    })
    .and_then(|r| r)
    .map_err(fail)?;

    Ok(RMatrix::new_matrix(n_real, k, |r, c| {
        neighbors[r][c] as i32 + 1
    }))
}

/// pANN for the first `n_real` rows.
/// @noRd
#[extendr]
fn rs_compute_pann(pcs: RMatrix<f64>, n_real: i32, k: i32, threads: i32) -> Result<Vec<f64>> {
    let (data, n_points, dims) = row_major(&pcs);
    let (n_real, k) = (to_usize(n_real, "n_real")?, to_usize(k, "k")?);
    with_threads(to_usize(threads, "threads")?, || {
        let emb = Embedding::new(&data, n_points, dims)?;
        doublet_rs::pann_from_embedding(emb, n_real, k, KnnMethod::Exact)
    })
    .and_then(|r| r)
    .map_err(fail)
}

/// pANN for the first `n_real` rows at each k in `ks`: an n_real x length(ks)
/// matrix.
/// @noRd
#[extendr]
fn rs_pann_sweep(
    pcs: RMatrix<f64>,
    n_real: i32,
    ks: Vec<i32>,
    threads: i32,
) -> Result<RMatrix<f64>> {
    let (data, n_points, dims) = row_major(&pcs);
    let n_real = to_usize(n_real, "n_real")?;
    let ks: Vec<usize> = ks
        .iter()
        .map(|&k| to_usize(k, "k"))
        .collect::<Result<_>>()?;
    let sweep = with_threads(to_usize(threads, "threads")?, || {
        let emb = Embedding::new(&data, n_points, dims)?;
        doublet_rs::pann_sweep(emb, n_real, &ks)
    })
    .and_then(|r| r)
    .map_err(fail)?;

    Ok(RMatrix::new_matrix(n_real, ks.len(), |r, c| sweep[c][r]))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod doubletrs;
    fn rs_build_doublets;
    fn rs_find_neighbors;
    fn rs_compute_pann;
    fn rs_pann_sweep;
}
