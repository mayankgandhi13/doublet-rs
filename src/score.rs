//! pANN scoring and doublet calling.

use rayon::prelude::*;

/// pANN for each real cell: the fraction of its `k` neighbours that are
/// artificial doublets (point index ≥ `n_real`).
///
/// Divides by `k` rather than the neighbour list length, as DoubletFinder
/// does, so a short list from approximate search can only lower the score.
pub fn compute_pann(neighbors: &[Vec<usize>], n_real: usize, k: usize) -> Vec<f64> {
    neighbors
        .par_iter()
        .map(|nn| nn.iter().filter(|&&j| j >= n_real).count() as f64 / k as f64)
        .collect()
}

/// Expected number of doublets: `round(rate × n_cells)`, rounding halves to
/// even like R.
pub fn n_expected_doublets(n_cells: usize, doublet_rate: f64) -> usize {
    ((n_cells as f64) * doublet_rate).round_ties_even() as usize
}

/// Flags the `n_doublets` cells with the highest pANN as doublets.
///
/// Ties keep the lower cell index first, matching R's stable
/// `order(pANN, decreasing = TRUE)`.
pub fn call_doublets(pann: &[f64], n_doublets: usize) -> Vec<bool> {
    let mut order: Vec<usize> = (0..pann.len()).collect();
    order.sort_by(|&a, &b| pann[b].total_cmp(&pann[a]));
    let mut calls = vec![false; pann.len()];
    for &i in order.iter().take(n_doublets) {
        calls[i] = true;
    }
    calls
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pann_counts_artificial_neighbours() {
        // 3 real cells (0..3), artificial points are 3 and 4.
        let nn = vec![vec![1, 3, 4], vec![0, 2, 3], vec![0, 1, 2]];
        assert_eq!(compute_pann(&nn, 3, 3), vec![2.0 / 3.0, 1.0 / 3.0, 0.0]);
    }

    #[test]
    fn calls_top_n_with_stable_ties() {
        let pann = [0.1, 0.5, 0.5, 0.9, 0.0];
        assert_eq!(
            call_doublets(&pann, 2),
            vec![false, true, false, true, false]
        );
        assert_eq!(call_doublets(&pann, 0), vec![false; 5]);
        assert_eq!(call_doublets(&pann, 10), vec![true; 5]);
    }

    #[test]
    fn expected_doublets_rounds() {
        assert_eq!(n_expected_doublets(1000, 0.075), 75);
        assert_eq!(n_expected_doublets(10, 0.05), 0); // R rounds 0.5 to even
        assert_eq!(n_expected_doublets(30, 0.05), 2); // 1.5 -> 2
    }
}
