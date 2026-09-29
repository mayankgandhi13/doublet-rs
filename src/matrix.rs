//! Minimal matrix types: a compressed sparse column (CSC) matrix for counts
//! and a dense row-major embedding for PCs.
//!
//! CSC types are generic over the index type so a caller can lend its own
//! buffers without converting, e.g. R's `dgCMatrix`, which stores 32-bit
//! `i` and `p` slots.

use std::fmt::Debug;

use crate::{Error, Result};

/// Integer types usable as sparse row indices and column pointers.
pub trait SparseIndex: Copy + Ord + Debug + Send + Sync {
    /// `None` if negative or too large for `usize`.
    fn to_usize(self) -> Option<usize>;
    /// `None` if `value` does not fit.
    fn from_usize(value: usize) -> Option<Self>;
}

macro_rules! impl_sparse_index {
    ($($t:ty),*) => {$(
        impl SparseIndex for $t {
            fn to_usize(self) -> Option<usize> {
                usize::try_from(self).ok()
            }
            fn from_usize(value: usize) -> Option<Self> {
                <$t>::try_from(value).ok()
            }
        }
    )*};
}
impl_sparse_index!(usize, u32, u64, i32, i64);

/// Borrowed CSC matrix. Columns are cells for gene × cell count matrices.
#[derive(Debug, Clone, Copy)]
pub struct CscView<'a, I> {
    n_rows: usize,
    n_cols: usize,
    indptr: &'a [I],
    indices: &'a [I],
    data: &'a [f64],
}

impl<'a, I: SparseIndex> CscView<'a, I> {
    /// Checks the CSC invariants: `indptr` has `n_cols + 1` non-decreasing
    /// entries from 0 to `nnz`, and each column's row indices are strictly
    /// increasing and below `n_rows`.
    pub fn new(
        n_rows: usize,
        n_cols: usize,
        indptr: &'a [I],
        indices: &'a [I],
        data: &'a [f64],
    ) -> Result<Self> {
        let invalid =
            |msg: String| Err(Error::InvalidArgument(format!("invalid CSC matrix: {msg}")));
        if indptr.len() != n_cols + 1 {
            return invalid(format!(
                "indptr has {} entries, expected {}",
                indptr.len(),
                n_cols + 1
            ));
        }
        if indices.len() != data.len() {
            return invalid(format!(
                "{} indices but {} values",
                indices.len(),
                data.len()
            ));
        }
        let ptr: Vec<usize> = match indptr.iter().map(|p| p.to_usize()).collect() {
            Some(p) => p,
            None => return invalid("negative column pointer".into()),
        };
        if ptr[0] != 0 || ptr[n_cols] != indices.len() {
            return invalid("column pointers must run from 0 to the number of non-zeros".into());
        }
        for col in 0..n_cols {
            if ptr[col] > ptr[col + 1] {
                return invalid(format!("column pointers decrease at column {col}"));
            }
            let mut prev: Option<usize> = None;
            for idx in &indices[ptr[col]..ptr[col + 1]] {
                match idx.to_usize() {
                    Some(row) if row < n_rows && prev.map_or(true, |p| row > p) => prev = Some(row),
                    _ => {
                        return invalid(format!(
                            "row indices in column {col} must be increasing and below {n_rows}"
                        ))
                    }
                }
            }
        }
        Ok(Self {
            n_rows,
            n_cols,
            indptr,
            indices,
            data,
        })
    }

    pub fn n_rows(&self) -> usize {
        self.n_rows
    }

    pub fn n_cols(&self) -> usize {
        self.n_cols
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// Row indices and values of column `col`.
    ///
    /// # Panics
    /// If `col >= n_cols`.
    pub fn column(&self, col: usize) -> (&'a [I], &'a [f64]) {
        // Validated in `new`, so these conversions cannot fail.
        let start = self.indptr[col].to_usize().unwrap();
        let end = self.indptr[col + 1].to_usize().unwrap();
        (&self.indices[start..end], &self.data[start..end])
    }

    /// Value at (`row`, `col`), zero if not stored.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        let (idx, vals) = self.column(col);
        idx.iter()
            .position(|r| r.to_usize() == Some(row))
            .map_or(0.0, |p| vals[p])
    }
}

/// Owned CSC matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct CscMatrix<I> {
    n_rows: usize,
    n_cols: usize,
    indptr: Vec<I>,
    indices: Vec<I>,
    data: Vec<f64>,
}

impl<I: SparseIndex> CscMatrix<I> {
    pub fn new(
        n_rows: usize,
        n_cols: usize,
        indptr: Vec<I>,
        indices: Vec<I>,
        data: Vec<f64>,
    ) -> Result<Self> {
        CscView::new(n_rows, n_cols, &indptr, &indices, &data)?;
        Ok(Self {
            n_rows,
            n_cols,
            indptr,
            indices,
            data,
        })
    }

    pub fn view(&self) -> CscView<'_, I> {
        CscView {
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            indptr: &self.indptr,
            indices: &self.indices,
            data: &self.data,
        }
    }

    /// `(indptr, indices, data)`.
    pub fn into_parts(self) -> (Vec<I>, Vec<I>, Vec<f64>) {
        (self.indptr, self.indices, self.data)
    }
}

impl CscMatrix<usize> {
    /// Builds a CSC matrix from `(row, col, value)` triplets in any order;
    /// duplicate positions are summed.
    pub fn from_triplets(
        n_rows: usize,
        n_cols: usize,
        mut triplets: Vec<(usize, usize, f64)>,
    ) -> Result<Self> {
        if let Some(&(r, c, _)) = triplets
            .iter()
            .find(|&&(r, c, _)| r >= n_rows || c >= n_cols)
        {
            return Err(Error::InvalidArgument(format!(
                "entry ({r}, {c}) outside a {n_rows} x {n_cols} matrix"
            )));
        }
        triplets.sort_unstable_by_key(|&(r, c, _)| (c, r));

        let mut indptr = vec![0; n_cols + 1];
        let mut indices = Vec::with_capacity(triplets.len());
        let mut data: Vec<f64> = Vec::with_capacity(triplets.len());
        let mut last: Option<(usize, usize)> = None;
        for (r, c, v) in triplets {
            if last == Some((r, c)) {
                *data.last_mut().unwrap() += v;
                continue;
            }
            indices.push(r);
            data.push(v);
            indptr[c + 1] += 1;
            last = Some((r, c));
        }
        for c in 0..n_cols {
            indptr[c + 1] += indptr[c];
        }
        Ok(Self {
            n_rows,
            n_cols,
            indptr,
            indices,
            data,
        })
    }
}

/// Dense points × dims embedding in row-major order, e.g. merged PCs with
/// real cells in the first rows.
#[derive(Debug, Clone, Copy)]
pub struct Embedding<'a> {
    data: &'a [f64],
    n_points: usize,
    dims: usize,
}

impl<'a> Embedding<'a> {
    pub fn new(data: &'a [f64], n_points: usize, dims: usize) -> Result<Self> {
        if dims == 0 {
            return Err(Error::InvalidArgument("embedding has no dimensions".into()));
        }
        if data.len() != n_points * dims {
            return Err(Error::InvalidArgument(format!(
                "embedding has {} values, expected {n_points} x {dims}",
                data.len()
            )));
        }
        if data.iter().any(|x| !x.is_finite()) {
            return Err(Error::InvalidArgument(
                "embedding contains NaN or infinite values".into(),
            ));
        }
        Ok(Self {
            data,
            n_points,
            dims,
        })
    }

    pub fn n_points(&self) -> usize {
        self.n_points
    }

    pub fn dims(&self) -> usize {
        self.dims
    }

    pub fn row(&self, i: usize) -> &'a [f64] {
        &self.data[i * self.dims..(i + 1) * self.dims]
    }
}

/// Converts a column-major (R / Fortran order) matrix to row-major.
pub fn col_major_to_row_major(data: &[f64], n_rows: usize, n_cols: usize) -> Vec<f64> {
    let mut out = vec![0.0; data.len()];
    for c in 0..n_cols {
        for r in 0..n_rows {
            out[r * n_cols + c] = data[c * n_rows + r];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_csc_view() {
        // [[1, 0], [0, 2], [3, 0]]
        let v = CscView::new(3, 2, &[0i32, 2, 3], &[0, 2, 1], &[1.0, 3.0, 2.0]).unwrap();
        assert_eq!(v.nnz(), 3);
        assert_eq!(v.column(0), (&[0i32, 2][..], &[1.0, 3.0][..]));
        assert_eq!(v.get(1, 1), 2.0);
        assert_eq!(v.get(0, 1), 0.0);
    }

    #[test]
    fn rejects_invalid_csc() {
        let d = [1.0, 2.0];
        assert!(CscView::new(3, 2, &[0i32, 1], &[0, 1], &d).is_err()); // short indptr
        assert!(CscView::new(3, 2, &[0i32, 1, 3], &[0, 1], &d).is_err()); // bad nnz
        assert!(CscView::new(3, 1, &[0i32, 2], &[1, 1], &d).is_err()); // repeated row
        assert!(CscView::new(3, 1, &[0i32, 2], &[2, 1], &d).is_err()); // unsorted
        assert!(CscView::new(3, 1, &[0i32, 2], &[0, 3], &d).is_err()); // row out of range
        assert!(CscView::new(3, 2, &[0i32, -1, 2], &[0, 1], &d).is_err()); // negative
    }

    #[test]
    fn triplets_are_sorted_and_summed() {
        let m = CscMatrix::from_triplets(
            2,
            2,
            vec![(1, 1, 4.0), (0, 0, 1.0), (1, 1, 1.0), (1, 0, 2.0)],
        )
        .unwrap();
        let v = m.view();
        assert_eq!(v.get(0, 0), 1.0);
        assert_eq!(v.get(1, 0), 2.0);
        assert_eq!(v.get(1, 1), 5.0);
        assert_eq!(v.nnz(), 3);
        assert!(CscMatrix::from_triplets(2, 2, vec![(2, 0, 1.0)]).is_err());
    }

    #[test]
    fn embedding_checks_shape_and_values() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let e = Embedding::new(&data, 3, 2).unwrap();
        assert_eq!(e.row(1), &[3.0, 4.0]);
        assert!(Embedding::new(&data, 4, 2).is_err());
        assert!(Embedding::new(&data, 6, 0).is_err());
        assert!(Embedding::new(&[1.0, f64::NAN], 1, 2).is_err());
    }

    #[test]
    fn transposes_column_major() {
        // R matrix(1:6, 2, 3) is column-major [1,2 | 3,4 | 5,6].
        let cm = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(
            col_major_to_row_major(&cm, 2, 3),
            vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]
        );
    }
}
