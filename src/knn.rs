//! k-nearest-neighbour search over the merged PC embedding.
//!
//! Backends:
//! - [`KnnMethod::Exact`]: parallel brute force, O(n_real × n_points). Gives
//!   the same neighbours as DoubletFinder, which also searches exactly.
//! - `KnnMethod::Hnsw` (cargo feature `hnsw`): approximate search with an
//!   HNSW graph. On the 16 benchmark datasets (up to 26k cells, 10 PCs) it
//!   was slower than exact search, so it is off by default.
//!
//! Both return, for each real cell, its `k` nearest neighbours among all
//! points (real + artificial), excluding the cell itself, ordered by
//! increasing distance.

use rayon::prelude::*;

use crate::matrix::Embedding;
use crate::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KnnMethod {
    Exact,
    #[cfg(feature = "hnsw")]
    Hnsw(HnswParams),
}

/// Finds the `k` nearest neighbours of each of the first `n_real` points.
pub fn find_neighbors(
    embedding: Embedding<'_>,
    n_real: usize,
    k: usize,
    method: KnnMethod,
) -> Result<Vec<Vec<usize>>> {
    validate(embedding, n_real, k)?;
    match method {
        KnnMethod::Exact => Ok((0..n_real)
            .into_par_iter()
            .map(|i| nearest(embedding, i, k))
            .collect()),
        #[cfg(feature = "hnsw")]
        KnnMethod::Hnsw(params) => hnsw::search(embedding, n_real, k, params),
    }
}

/// pANN of every real cell at several neighbourhood sizes at once, for pK
/// parameter sweeps.
///
/// Each cell's neighbours are found once, up to the largest `k`, and the
/// artificial-neighbour count is read off at each `k`. Returns one pANN
/// vector per entry of `ks`, in the same order; each equals
/// [`compute_pann`](crate::compute_pann) on exact neighbours at that `k`.
pub fn pann_sweep(embedding: Embedding<'_>, n_real: usize, ks: &[usize]) -> Result<Vec<Vec<f64>>> {
    let k_max = ks.iter().copied().max().unwrap_or(0);
    validate(embedding, n_real, k_max)?;
    if ks.contains(&0) {
        return Err(Error::InvalidArgument("every k must be at least 1".into()));
    }

    let per_cell: Vec<Vec<f64>> = (0..n_real)
        .into_par_iter()
        .map(|i| {
            let nn = nearest(embedding, i, k_max);
            // artificial[m] = number of artificial points among the first m.
            let mut artificial = Vec::with_capacity(k_max + 1);
            artificial.push(0usize);
            for &j in &nn {
                artificial.push(artificial.last().unwrap() + usize::from(j >= n_real));
            }
            ks.iter()
                .map(|&k| artificial[k] as f64 / k as f64)
                .collect()
        })
        .collect();

    Ok((0..ks.len())
        .map(|s| per_cell.iter().map(|cell| cell[s]).collect())
        .collect())
}

fn validate(embedding: Embedding<'_>, n_real: usize, k: usize) -> Result<()> {
    let n_points = embedding.n_points();
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
    Ok(())
}

/// Exact `k` nearest neighbours of point `i`, excluding itself. Ties break
/// on index so results are deterministic.
fn nearest(embedding: Embedding<'_>, i: usize, k: usize) -> Vec<usize> {
    let query = embedding.row(i);
    let mut dists: Vec<(f64, usize)> = (0..embedding.n_points())
        .filter(|&j| j != i)
        .map(|j| {
            let d = query
                .iter()
                .zip(embedding.row(j))
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>();
            (d, j)
        })
        .collect();
    let cmp = |a: &(f64, usize), b: &(f64, usize)| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1));
    dists.select_nth_unstable_by(k - 1, cmp);
    dists.truncate(k);
    dists.sort_unstable_by(cmp);
    dists.into_iter().map(|(_, j)| j).collect()
}

#[cfg(feature = "hnsw")]
pub use hnsw::HnswParams;

#[cfg(feature = "hnsw")]
mod hnsw {
    use hnsw_rs::prelude::{DistL2, Hnsw};
    use rayon::prelude::*;

    use crate::matrix::Embedding;
    use crate::{Error, Result};

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

    const MAX_LAYERS: usize = 16;

    pub(super) fn search(
        embedding: Embedding<'_>,
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
        let n_points = embedding.n_points();
        let points: Vec<Vec<f32>> = (0..n_points)
            .map(|i| embedding.row(i).iter().map(|&x| x as f32).collect())
            .collect();

        let index = Hnsw::<f32, DistL2>::new(
            params.max_connections,
            n_points,
            MAX_LAYERS,
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute_pann;
    use rand::{Rng, SeedableRng};

    fn random_points(n: usize, dims: usize, seed: u64) -> Vec<f64> {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        (0..n * dims).map(|_| rng.random_range(-1.0..1.0)).collect()
    }

    #[test]
    fn exact_finds_nearest_on_a_line() {
        // Points at 0, 1, 3, 7, 15 on a line.
        let data = [0.0, 1.0, 3.0, 7.0, 15.0];
        let emb = Embedding::new(&data, 5, 1).unwrap();
        let nn = find_neighbors(emb, 3, 2, KnnMethod::Exact).unwrap();
        assert_eq!(nn, vec![vec![1, 2], vec![0, 2], vec![1, 0]]);
    }

    #[test]
    fn exact_breaks_ties_by_index() {
        let data = [0.0, 1.0, -1.0];
        let emb = Embedding::new(&data, 3, 1).unwrap();
        assert_eq!(
            find_neighbors(emb, 1, 1, KnnMethod::Exact).unwrap(),
            vec![vec![1]]
        );
    }

    #[test]
    fn rejects_bad_arguments() {
        let data = [0.0, 1.0, 2.0];
        let emb = Embedding::new(&data, 3, 1).unwrap();
        assert!(find_neighbors(emb, 4, 1, KnnMethod::Exact).is_err());
        assert!(find_neighbors(emb, 3, 0, KnnMethod::Exact).is_err());
        assert!(find_neighbors(emb, 3, 3, KnnMethod::Exact).is_err());
        assert!(pann_sweep(emb, 3, &[1, 0]).is_err());
        assert!(pann_sweep(emb, 3, &[3]).is_err());
    }

    #[test]
    fn sweep_matches_single_k_pann() {
        let (n, dims, n_real) = (400, 5, 300);
        let data = random_points(n, dims, 1);
        let emb = Embedding::new(&data, n, dims).unwrap();
        let ks = [1, 5, 17, 40, 399];
        let sweep = pann_sweep(emb, n_real, &ks).unwrap();
        for (s, &k) in ks.iter().enumerate() {
            let nn = find_neighbors(emb, n_real, k, KnnMethod::Exact).unwrap();
            assert_eq!(sweep[s], compute_pann(&nn, n_real, k), "k = {k}");
        }
    }

    #[cfg(feature = "hnsw")]
    #[test]
    fn hnsw_agrees_with_exact_on_random_points() {
        let (n, dims, k) = (2000, 10, 15);
        let data = random_points(n, dims, 42);
        let emb = Embedding::new(&data, n, dims).unwrap();
        let exact = find_neighbors(emb, 500, k, KnnMethod::Exact).unwrap();
        let approx = find_neighbors(emb, 500, k, KnnMethod::Hnsw(HnswParams::default())).unwrap();
        let hits: usize = exact
            .iter()
            .zip(&approx)
            .map(|(e, a)| a.iter().filter(|j| e.contains(j)).count())
            .sum();
        let recall = hits as f64 / (500 * k) as f64;
        assert!(recall > 0.95, "HNSW recall too low: {recall}");
        let bad = HnswParams {
            max_connections: 300,
            ..Default::default()
        };
        assert!(find_neighbors(emb, 1, 1, KnnMethod::Hnsw(bad)).is_err());
    }
}
