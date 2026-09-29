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
pub mod matrix;
pub mod score;
pub mod simulate;

use std::fmt;

pub use knn::{find_neighbors, pann_sweep, KnnMethod};
pub use matrix::{CscMatrix, CscView, Embedding, SparseIndex};
pub use score::{call_doublets, compute_pann, n_expected_doublets};
pub use simulate::{n_doublets_for_pn, simulate_doublets};

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    Parse(String),
    InvalidArgument(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Parse(msg) => write!(f, "parse error: {msg}"),
            Error::InvalidArgument(msg) => write!(f, "invalid argument: {msg}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

pub type Result<T> = std::result::Result<T, Error>;

/// Neighbourhood size used by DoubletFinder: `round(n_points * pK)`, where
/// `n_points` counts real cells plus artificial doublets.
pub fn k_from_pk(n_points: usize, pk: f64) -> usize {
    round_half_even((n_points as f64) * pk) as usize
}

/// Runs steps 3 and 4 on a merged embedding: returns pANN for each real cell.
pub fn pann_from_embedding(
    embedding: Embedding<'_>,
    n_real: usize,
    k: usize,
    method: KnnMethod,
) -> Result<Vec<f64>> {
    let neighbors = find_neighbors(embedding, n_real, k, method)?;
    Ok(compute_pann(&neighbors, n_real, k))
}

/// Runs `f` on a thread pool of `threads` workers, so callers (and CRAN
/// checks) control how many cores the parallel steps use.
pub fn with_threads<T: Send>(threads: usize, f: impl FnOnce() -> T + Send) -> Result<T> {
    if threads == 0 {
        return Err(Error::InvalidArgument("threads must be at least 1".into()));
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| Error::InvalidArgument(format!("could not start thread pool: {e}")))?;
    Ok(pool.install(f))
}

/// Rounds halves to even, like R's `round()`: `round(0.5) == 0`,
/// `round(1.5) == 2`. Rust's `f64::round` rounds halves away from zero.
pub(crate) fn round_half_even(x: f64) -> f64 {
    if (x - x.trunc()).abs() == 0.5 {
        2.0 * (x / 2.0).round()
    } else {
        x.round()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounds_like_r() {
        let cases = [
            (0.5, 0.0),
            (1.5, 2.0),
            (2.5, 2.0),
            (-0.5, 0.0),
            (2.4, 2.0),
            (2.6, 3.0),
        ];
        for (x, want) in cases {
            assert_eq!(round_half_even(x), want, "round({x})");
        }
    }

    #[test]
    fn thread_pool_runs_and_rejects_zero() {
        assert_eq!(with_threads(2, rayon::current_num_threads).unwrap(), 2);
        assert!(with_threads(0, || ()).is_err());
    }
}
