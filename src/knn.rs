//! k-nearest-neighbour search over the merged PC embedding.
//!
//! Two backends:
//! - [`KnnMethod::Exact`]: brute force, O(n_real × n_points). The reference
//!   for validating against DoubletFinder, which also searches exactly.
//! - [`KnnMethod::Hnsw`]: approximate search with an HNSW graph built once
//!   and queried in parallel. The one to use on large datasets.
//!
//! Both return, for each real cell, its `k` nearest neighbours among all
//! points (real + artificial), excluding the cell itself, ordered by
//! increasing distance.

use hnsw_rs::prelude::{DistL2, Hnsw};
use ndarray::ArrayView2;
use rayon::prelude::*;

use crate::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KnnMethod {
    Exact,
    Hnsw(HnswParams),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HnswParams {
    /// Links per node per layer (HNSW `M`). Must be ≤ 256.
    pub max_connections: usize,
    /// Candidate list size while building the graph.
    pub ef_construction: usize,
    /// Candidate list size while searching; raised to `k + 1` if smaller.
    pub ef_search: usize,
}

impl Default for HnswParams {
    fn default() -> Self {
        Self {
            max_connections: 24,
            ef_construction: 200,
            ef_search: 128,
        }
    }
}

const HNSW_MAX_LAYERS: usize = 16;

/// Finds the `k` nearest neighbours of each of the first `n_real` points.
///
/// `embedding` is `points × dims` with real cells in rows `0..n_real`.
pub fn find_neighbors(
    embedding: ArrayView2<'_, f64>,
    n_real: usize,
    k: usize,
    method: KnnMethod,
) -> Result<Vec<Vec<usize>>> {
    validate(embedding, n_real, k)?;
    match method {
        KnnMethod::Exact => Ok(exact(embedding, n_real, k)),
        KnnMethod::Hnsw(params) => hnsw(embedding, n_real, k, params),
    }
}

fn validate(embedding: ArrayView2<'_, f64>, n_real: usize, k: usize) -> Result<()> {
    let (n_points, dims) = embedding.dim();
    if dims == 0 {
        return Err(Error::InvalidArgument("embedding has no dimensions".into()));
    }
    if n_real > n_points {
        return Err(Error::InvalidArgument(format!(
            "n_real ({n_real}) exceeds number of points ({n_points})"
        )));
    }
    if k == 0 || k >= n_points {
        return Err(Error::InvalidArgument(format!(
            "k must be in 1..{n_points} (number of points), got {k}"
        )));
    }
    if embedding.iter().any(|x| !x.is_finite()) {
        return Err(Error::InvalidArgument(
            "embedding contains NaN or infinite values".into(),
        ));
    }
    Ok(())
}

fn exact(embedding: ArrayView2<'_, f64>, n_real: usize, k: usize) -> Vec<Vec<usize>> {
    let n_points = embedding.nrows();
    (0..n_real)
        .into_par_iter()
        .map(|i| {
            let query = embedding.row(i);
            let mut dists: Vec<(f64, usize)> = (0..n_points)
                .filter(|&j| j != i)
                .map(|j| {
                    let d = query
                        .iter()
                        .zip(embedding.row(j).iter())
                        .map(|(a, b)| (a - b) * (a - b))
                        .sum::<f64>();
                    (d, j)
                })
                .collect();
            // Ties break on index so results are deterministic.
            let cmp = |a: &(f64, usize), b: &(f64, usize)| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1));
            dists.select_nth_unstable_by(k - 1, cmp);
            dists.truncate(k);
            dists.sort_unstable_by(cmp);
            dists.into_iter().map(|(_, j)| j).collect()
        })
        .collect()
}

fn hnsw(
    embedding: ArrayView2<'_, f64>,
    n_real: usize,
    k: usize,
    params: HnswParams,
) -> Result<Vec<Vec<usize>>> {
    if params.max_connections == 0 || params.max_connections > 256 {
        // hnsw_rs exits the process on > 256, so check before calling it.
        return Err(Error::InvalidArgument(format!(
            "max_connections must be in 1..=256, got {}",
            params.max_connections
        )));
    }
    let n_points = embedding.nrows();

    // f32 halves the index's memory and uses hnsw_rs's SIMD L2 kernel.
    let points: Vec<Vec<f32>> = embedding
        .rows()
        .into_iter()
        .map(|row| row.iter().map(|&x| x as f32).collect())
        .collect();

    let index = Hnsw::<f32, DistL2>::new(
        params.max_connections,
        n_points,
        HNSW_MAX_LAYERS,
        params.ef_construction,
        DistL2,
    );
    let with_ids: Vec<(&Vec<f32>, usize)> = points.iter().zip(0..).collect();
    index.parallel_insert(&with_ids);

    let ef = params.ef_search.max(k + 1);
    Ok((0..n_real)
        .into_par_iter()
        .map(|i| {
            // Ask for one extra so we can drop the query point itself.
            index
                .search(&points[i], k + 1, ef)
                .into_iter()
                .map(|n| n.d_id)
                .filter(|&j| j != i)
                .take(k)
                .collect()
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array2, array};

    #[test]
    fn exact_finds_nearest_on_a_line() {
        // Points at 0, 1, 3, 7, 15 on a line.
        let emb = array![[0.0], [1.0], [3.0], [7.0], [15.0]];
        let nn = find_neighbors(emb.view(), 3, 2, KnnMethod::Exact).unwrap();
        assert_eq!(nn, vec![vec![1, 2], vec![0, 2], vec![1, 0]]);
    }

    #[test]
    fn exact_breaks_ties_by_index() {
        let emb = array![[0.0], [1.0], [-1.0]];
        let nn = find_neighbors(emb.view(), 1, 1, KnnMethod::Exact).unwrap();
        assert_eq!(nn, vec![vec![1]]);
    }

    #[test]
    fn rejects_bad_arguments() {
        let emb = array![[0.0], [1.0], [2.0]];
        assert!(find_neighbors(emb.view(), 4, 1, KnnMethod::Exact).is_err());
        assert!(find_neighbors(emb.view(), 3, 0, KnnMethod::Exact).is_err());
        assert!(find_neighbors(emb.view(), 3, 3, KnnMethod::Exact).is_err());
        let nan = array![[0.0], [f64::NAN]];
        assert!(find_neighbors(nan.view(), 1, 1, KnnMethod::Exact).is_err());
        let bad = HnswParams {
            max_connections: 300,
            ..Default::default()
        };
        assert!(find_neighbors(emb.view(), 3, 1, KnnMethod::Hnsw(bad)).is_err());
    }

    #[test]
    fn hnsw_agrees_with_exact_on_random_points() {
        use rand::{RngExt, SeedableRng};
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(42);
        let (n, dims, k) = (2000, 10, 15);
        let emb = Array2::from_shape_fn((n, dims), |_| rng.random_range(-1.0..1.0));

        let exact = find_neighbors(emb.view(), 500, k, KnnMethod::Exact).unwrap();
        let approx =
            find_neighbors(emb.view(), 500, k, KnnMethod::Hnsw(HnswParams::default())).unwrap();

        let hits: usize = exact
            .iter()
            .zip(&approx)
            .map(|(e, a)| a.iter().filter(|j| e.contains(j)).count())
            .sum();
        let recall = hits as f64 / (500 * k) as f64;
        assert!(recall > 0.95, "HNSW recall too low: {recall}");
        assert!(
            approx
                .iter()
                .enumerate()
                .all(|(i, nn)| nn.len() == k && !nn.contains(&i))
        );
    }
}
