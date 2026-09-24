//! Memory-efficient reimplementation of the DoubletFinder core engine.
//!
//! The pipeline is split so the host language (R / Python) keeps
//! normalization and PCA, while Rust handles the memory-heavy steps:
//!
//! 1. [`simulate::simulate_doublets`] builds artificial doublets from the
//!    raw gene × cell count matrix and returns them as a sparse matrix.
//! 2. The host merges real cells + doublets, normalizes, and runs PCA.
//! 3. [`knn`] finds the `k` nearest neighbours of every real cell in the
//!    merged PC embedding.
//! 4. [`score`] turns neighbourhoods into pANN scores and doublet calls.
//!
//! Conventions: embeddings are `points × dims`, with the `n_real` real
//! cells first and artificial doublets after them, as in DoubletFinder.

pub mod io;
pub mod knn;
pub mod score;
pub mod simulate;

use ndarray::ArrayView2;

pub use knn::{KnnMethod, find_neighbors};
pub use score::{call_doublets, compute_pann, n_expected_doublets};
pub use simulate::{n_doublets_for_pn, simulate_doublets};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error: {0}")]
    Parse(String),
    #[error("invalid argument: {0}")]
    InvalidArgument(String),
}

pub type Result<T> = std::result::Result<T, Error>;

/// Neighbourhood size used by DoubletFinder: `round(n_points * pK)`, where
/// `n_points` counts real cells plus artificial doublets.
pub fn k_from_pk(n_points: usize, pk: f64) -> usize {
    ((n_points as f64) * pk).round_ties_even() as usize
}

/// Runs steps 3 and 4 on a merged embedding: returns pANN for each real cell.
pub fn pann_from_embedding(
    embedding: ArrayView2<'_, f64>,
    n_real: usize,
    k: usize,
    method: KnnMethod,
) -> Result<Vec<f64>> {
    let neighbors = find_neighbors(embedding, n_real, k, method)?;
    Ok(compute_pann(&neighbors, n_real, k))
}
