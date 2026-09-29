//! Artificial doublet simulation.
//!
//! Like DoubletFinder, each doublet is the average of two cells sampled
//! independently and with replacement. Only the `(cell_a, cell_b)` pairs are
//! drawn up front (two integers per doublet); expression columns are built
//! in parallel, either all at once or in fixed-size batches so the full
//! doublet matrix never has to exist in memory.
//!
//! Functions are generic over the sparse index type so callers can pass
//! borrowed matrices without converting indices, e.g. R's `dgCMatrix`,
//! which stores 32-bit `i` and `p` slots.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::matrix::{CscMatrix, CscView, SparseIndex};
use crate::{round_half_even, Error, Result};

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
    Ok(round_half_even(n / (1.0 - pn) - n) as usize)
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
pub struct Doublets<I> {
    /// genes × doublets.
    pub matrix: CscMatrix<I>,
    pub pairs: Vec<(usize, usize)>,
}

/// Simulates `n_doublets` doublets from a genes × cells matrix.
pub fn simulate_doublets<I: SparseIndex>(
    expr: CscView<'_, I>,
    n_doublets: usize,
    seed: u64,
) -> Result<Doublets<I>> {
    let pairs = sample_pairs(expr.n_cols(), n_doublets, seed)?;
    let matrix = build_doublets(expr, &pairs)?;
    Ok(Doublets { matrix, pairs })
}

/// Builds the genes × `pairs.len()` doublet matrix for the given pairs.
pub fn build_doublets<I: SparseIndex>(
    expr: CscView<'_, I>,
    pairs: &[(usize, usize)],
) -> Result<CscMatrix<I>> {
    if let Some(&(a, b)) = pairs
        .iter()
        .find(|&&(a, b)| a >= expr.n_cols() || b >= expr.n_cols())
    {
        return Err(Error::InvalidArgument(format!(
            "pair ({a}, {b}) out of range for {} cells",
            expr.n_cols()
        )));
    }

    let columns: Vec<(Vec<I>, Vec<f64>)> = pairs
        .par_iter()
        .map(|&(a, b)| average_columns(expr, a, b))
        .collect();

    let nnz: usize = columns.iter().map(|(idx, _)| idx.len()).sum();
    let to_index = |n: usize| {
        I::from_usize(n).ok_or_else(|| {
            Error::InvalidArgument(format!("{nnz} non-zeros do not fit the matrix index type"))
        })
    };
    let mut indptr = Vec::with_capacity(columns.len() + 1);
    indptr.push(to_index(0)?);
    let mut indices = Vec::with_capacity(nnz);
    let mut data = Vec::with_capacity(nnz);
    for (idx, vals) in columns {
        indices.extend(idx);
        data.extend(vals);
        indptr.push(to_index(indices.len())?);
    }
    CscMatrix::new(expr.n_rows(), pairs.len(), indptr, indices, data)
}

/// Iterator over doublets in batches of at most `batch_size` columns.
///
/// Each batch is built in parallel; peak memory is one batch rather than
/// the whole doublet matrix.
pub struct DoubletBatches<'a, I> {
    expr: CscView<'a, I>,
    pairs: Vec<(usize, usize)>,
    batch_size: usize,
    next: usize,
}

impl<'a, I: SparseIndex> DoubletBatches<'a, I> {
    pub fn new(
        expr: CscView<'a, I>,
        n_doublets: usize,
        batch_size: usize,
        seed: u64,
    ) -> Result<Self> {
        if batch_size == 0 {
            return Err(Error::InvalidArgument("batch_size must be > 0".into()));
        }
        let pairs = sample_pairs(expr.n_cols(), n_doublets, seed)?;
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

impl<I: SparseIndex> Iterator for DoubletBatches<'_, I> {
    type Item = Result<CscMatrix<I>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.pairs.len() {
            return None;
        }
        let end = (self.next + self.batch_size).min(self.pairs.len());
        let batch = build_doublets(self.expr, &self.pairs[self.next..end]);
        self.next = end;
        Some(batch)
    }
}

/// Merges two sorted sparse columns into their element-wise average.
fn average_columns<I: SparseIndex>(expr: CscView<'_, I>, a: usize, b: usize) -> (Vec<I>, Vec<f64>) {
    let (ia, va) = expr.column(a);
    let (ib, vb) = expr.column(b);

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

    /// 3 genes × 3 cells, laid out like R's dgCMatrix slots:
    /// cell0 = [2, 0, 4], cell1 = [0, 6, 2], cell2 = [0, 0, 0]
    const P: [i32; 4] = [0, 2, 4, 4];
    const I: [i32; 4] = [0, 2, 1, 2];
    const X: [f64; 4] = [2.0, 4.0, 6.0, 2.0];

    fn toy() -> CscView<'static, i32> {
        CscView::new(3, 3, &P, &I, &X).unwrap()
    }

    fn dense_col<T: SparseIndex>(m: &CscMatrix<T>, c: usize) -> Vec<f64> {
        (0..3).map(|r| m.view().get(r, c)).collect()
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
        let m = build_doublets(toy(), &[(0, 1), (0, 2), (1, 1)]).unwrap();
        assert_eq!(dense_col(&m, 0), vec![1.0, 3.0, 3.0]);
        assert_eq!(dense_col(&m, 1), vec![1.0, 0.0, 2.0]);
        assert_eq!(dense_col(&m, 2), vec![0.0, 6.0, 2.0]);
    }

    #[test]
    fn same_seed_same_doublets() {
        let a = simulate_doublets(toy(), 50, 7).unwrap();
        let b = simulate_doublets(toy(), 50, 7).unwrap();
        let c = simulate_doublets(toy(), 50, 8).unwrap();
        assert_eq!(a.pairs, b.pairs);
        assert_eq!(a.matrix, b.matrix);
        assert_ne!(a.pairs, c.pairs);
    }

    #[test]
    fn index_type_does_not_change_results() {
        let p: Vec<usize> = P.iter().map(|&v| v as usize).collect();
        let i: Vec<usize> = I.iter().map(|&v| v as usize).collect();
        let wide = CscView::new(3, 3, &p, &i, &X).unwrap();
        let a = simulate_doublets(toy(), 20, 5).unwrap();
        let b = simulate_doublets(wide, 20, 5).unwrap();
        assert_eq!(a.pairs, b.pairs);
        for c in 0..20 {
            assert_eq!(dense_col(&a.matrix, c), dense_col(&b.matrix, c));
        }
    }

    #[test]
    fn batches_concatenate_to_full_matrix() {
        let full = simulate_doublets(toy(), 10, 3).unwrap();
        let batches: Vec<_> = DoubletBatches::new(toy(), 10, 4, 3)
            .unwrap()
            .map(|b| b.unwrap())
            .collect();
        assert_eq!(
            batches
                .iter()
                .map(|b| b.view().n_cols())
                .collect::<Vec<_>>(),
            vec![4, 4, 2]
        );

        let mut col = 0;
        for batch in &batches {
            for c in 0..batch.view().n_cols() {
                assert_eq!(dense_col(batch, c), dense_col(&full.matrix, col));
                col += 1;
            }
        }
    }

    #[test]
    fn rejects_out_of_range_pairs() {
        assert!(build_doublets(toy(), &[(0, 3)]).is_err());
    }
}
