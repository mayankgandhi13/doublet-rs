//! Artificial doublet simulation.
//!
//! Like DoubletFinder, each doublet is the average of two cells sampled
//! independently and with replacement. Only the `(cell_a, cell_b)` pairs are
//! drawn up front (two integers per doublet); expression columns are built
//! in parallel, either all at once or in fixed-size batches so the full
//! doublet matrix never has to exist in memory.

use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use sprs::{CsMat, CsMatView};

use crate::{Error, Result};

/// Number of doublets DoubletFinder simulates for a given `pN`, the
/// fraction of the merged (real + artificial) dataset that is artificial:
/// `round(n_real / (1 - pN) - n_real)`.
pub fn n_doublets_for_pn(n_real: usize, pn: f64) -> Result<usize> {
    if !(0.0..1.0).contains(&pn) {
        return Err(Error::InvalidArgument(format!(
            "pN must be in [0, 1), got {pn}"
        )));
    }
    let n = n_real as f64;
    Ok((n / (1.0 - pn) - n).round_ties_even() as usize)
}

/// Draws `n_doublets` cell pairs from `0..n_cells`, reproducibly for a seed.
pub fn sample_pairs(n_cells: usize, n_doublets: usize, seed: u64) -> Result<Vec<(usize, usize)>> {
    if n_cells == 0 {
        return Err(Error::InvalidArgument(
            "expression matrix has no cells".into(),
        ));
    }
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    Ok((0..n_doublets)
        .map(|_| (rng.random_range(0..n_cells), rng.random_range(0..n_cells)))
        .collect())
}

/// Simulated doublets plus the parent cells each one was built from.
#[derive(Debug, Clone)]
pub struct Doublets {
    /// genes × doublets, CSC.
    pub matrix: CsMat<f64>,
    pub pairs: Vec<(usize, usize)>,
}

/// Simulates `n_doublets` doublets from a genes × cells matrix.
///
/// `expr` should be CSC (columns = cells); a CSR matrix is converted first,
/// which costs one copy.
pub fn simulate_doublets(
    expr: CsMatView<'_, f64>,
    n_doublets: usize,
    seed: u64,
) -> Result<Doublets> {
    let pairs = sample_pairs(expr.cols(), n_doublets, seed)?;
    let matrix = with_csc(expr, |csc| build_doublets(csc, &pairs))?;
    Ok(Doublets { matrix, pairs })
}

/// Builds the genes × `pairs.len()` doublet matrix for the given pairs.
pub fn build_doublets(expr: CsMatView<'_, f64>, pairs: &[(usize, usize)]) -> Result<CsMat<f64>> {
    if !expr.is_csc() {
        return Err(Error::InvalidArgument(
            "expression matrix must be CSC".into(),
        ));
    }
    if let Some(&(a, b)) = pairs
        .iter()
        .find(|&&(a, b)| a >= expr.cols() || b >= expr.cols())
    {
        return Err(Error::InvalidArgument(format!(
            "pair ({a}, {b}) out of range for {} cells",
            expr.cols()
        )));
    }

    let columns: Vec<(Vec<usize>, Vec<f64>)> = pairs
        .par_iter()
        .map(|&(a, b)| average_columns(expr, a, b))
        .collect();

    let mut indptr = Vec::with_capacity(columns.len() + 1);
    indptr.push(0);
    let nnz = columns.iter().map(|(idx, _)| idx.len()).sum();
    let mut indices = Vec::with_capacity(nnz);
    let mut data = Vec::with_capacity(nnz);
    for (idx, vals) in columns {
        indices.extend(idx);
        data.extend(vals);
        indptr.push(indices.len());
    }
    Ok(CsMat::new_csc(
        (expr.rows(), pairs.len()),
        indptr,
        indices,
        data,
    ))
}

/// Iterator over doublets in batches of at most `batch_size` columns.
///
/// Each batch is built in parallel; peak memory is one batch rather than
/// the whole doublet matrix.
pub struct DoubletBatches<'a> {
    expr: CsMatView<'a, f64>,
    pairs: Vec<(usize, usize)>,
    batch_size: usize,
    next: usize,
}

impl<'a> DoubletBatches<'a> {
    pub fn new(
        expr: CsMatView<'a, f64>,
        n_doublets: usize,
        batch_size: usize,
        seed: u64,
    ) -> Result<Self> {
        if !expr.is_csc() {
            return Err(Error::InvalidArgument(
                "expression matrix must be CSC".into(),
            ));
        }
        if batch_size == 0 {
            return Err(Error::InvalidArgument("batch_size must be > 0".into()));
        }
        let pairs = sample_pairs(expr.cols(), n_doublets, seed)?;
        Ok(Self {
            expr,
            pairs,
            batch_size,
            next: 0,
        })
    }

    pub fn pairs(&self) -> &[(usize, usize)] {
        &self.pairs
    }
}

impl Iterator for DoubletBatches<'_> {
    type Item = CsMat<f64>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.pairs.len() {
            return None;
        }
        let end = (self.next + self.batch_size).min(self.pairs.len());
        let batch = build_doublets(self.expr, &self.pairs[self.next..end])
            .expect("pairs were validated when sampled");
        self.next = end;
        Some(batch)
    }
}

fn with_csc<T>(expr: CsMatView<'_, f64>, f: impl FnOnce(CsMatView<'_, f64>) -> T) -> T {
    if expr.is_csc() {
        f(expr)
    } else {
        let csc = expr.to_csc();
        f(csc.view())
    }
}

/// Merges two sorted sparse columns into their element-wise average.
fn average_columns(expr: CsMatView<'_, f64>, a: usize, b: usize) -> (Vec<usize>, Vec<f64>) {
    let col_a = expr.outer_view(a).expect("column index validated");
    let col_b = expr.outer_view(b).expect("column index validated");
    let (ia, va) = (col_a.indices(), col_a.data());
    let (ib, vb) = (col_b.indices(), col_b.data());

    let mut indices = Vec::with_capacity(ia.len() + ib.len());
    let mut values = Vec::with_capacity(ia.len() + ib.len());
    let (mut i, mut j) = (0, 0);
    while i < ia.len() || j < ib.len() {
        let (row, sum) = if j == ib.len() || (i < ia.len() && ia[i] < ib[j]) {
            i += 1;
            (ia[i - 1], va[i - 1])
        } else if i == ia.len() || ib[j] < ia[i] {
            j += 1;
            (ib[j - 1], vb[j - 1])
        } else {
            i += 1;
            j += 1;
            (ia[i - 1], va[i - 1] + vb[j - 1])
        };
        indices.push(row);
        values.push(sum / 2.0);
    }
    (indices, values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use sprs::TriMat;

    /// 3 genes × 3 cells:
    /// cell0 = [2, 0, 4], cell1 = [0, 6, 2], cell2 = [0, 0, 0]
    fn toy() -> CsMat<f64> {
        let mut t = TriMat::new((3, 3));
        t.add_triplet(0, 0, 2.0);
        t.add_triplet(2, 0, 4.0);
        t.add_triplet(1, 1, 6.0);
        t.add_triplet(2, 1, 2.0);
        t.to_csc()
    }

    fn dense_col(m: &CsMat<f64>, c: usize) -> Vec<f64> {
        (0..m.rows())
            .map(|r| *m.get(r, c).unwrap_or(&0.0))
            .collect()
    }

    #[test]
    fn n_doublets_matches_doubletfinder() {
        // DoubletFinder: round(1000 / (1 - 0.25) - 1000) = 333
        assert_eq!(n_doublets_for_pn(1000, 0.25).unwrap(), 333);
        assert_eq!(n_doublets_for_pn(1000, 0.0).unwrap(), 0);
        assert!(n_doublets_for_pn(1000, 1.0).is_err());
    }

    #[test]
    fn doublet_is_average_of_parents() {
        let m = build_doublets(toy().view(), &[(0, 1), (0, 2), (1, 1)]).unwrap();
        assert_eq!(m.shape(), (3, 3));
        assert_eq!(dense_col(&m, 0), vec![1.0, 3.0, 3.0]);
        assert_eq!(dense_col(&m, 1), vec![1.0, 0.0, 2.0]);
        assert_eq!(dense_col(&m, 2), vec![0.0, 6.0, 2.0]);
    }

    #[test]
    fn same_seed_same_doublets() {
        let expr = toy();
        let a = simulate_doublets(expr.view(), 50, 7).unwrap();
        let b = simulate_doublets(expr.view(), 50, 7).unwrap();
        let c = simulate_doublets(expr.view(), 50, 8).unwrap();
        assert_eq!(a.pairs, b.pairs);
        assert_eq!(a.matrix, b.matrix);
        assert_ne!(a.pairs, c.pairs);
    }

    #[test]
    fn csr_input_gives_same_result() {
        let expr = toy();
        let csr = expr.to_csr();
        let a = simulate_doublets(expr.view(), 20, 1).unwrap();
        let b = simulate_doublets(csr.view(), 20, 1).unwrap();
        assert_eq!(a.matrix, b.matrix);
    }

    #[test]
    fn batches_concatenate_to_full_matrix() {
        let expr = toy();
        let full = simulate_doublets(expr.view(), 10, 3).unwrap();
        let batches: Vec<_> = DoubletBatches::new(expr.view(), 10, 4, 3)
            .unwrap()
            .collect();
        assert_eq!(
            batches.iter().map(|b| b.cols()).collect::<Vec<_>>(),
            vec![4, 4, 2]
        );

        let mut col = 0;
        for batch in &batches {
            for c in 0..batch.cols() {
                for (x, y) in dense_col(batch, c).iter().zip(dense_col(&full.matrix, col)) {
                    assert_relative_eq!(*x, y);
                }
                col += 1;
            }
        }
    }

    #[test]
    fn rejects_out_of_range_pairs() {
        assert!(build_doublets(toy().view(), &[(0, 3)]).is_err());
    }
}
