//! Matrix Market (`.mtx`) reading, e.g. 10x Genomics `matrix.mtx(.gz)` output.
//!
//! `sprs::io::read_matrix_market` rejects `integer` files when asked for
//! `f64` values, and 10x count matrices are always `integer`, so we parse
//! the coordinate format ourselves and stream entries straight into a
//! triplet matrix.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use sprs::{CsMat, TriMat};

use crate::{Error, Result};

/// Reads a `coordinate` Matrix Market file into a CSC matrix.
///
/// 10x files are gene × cell, so columns are cells. Supports `integer`,
/// `real`, and `pattern` fields with `general` symmetry.
pub fn read_mtx<P: AsRef<Path>>(path: P) -> Result<CsMat<f64>> {
    let file = File::open(path)?;
    read_mtx_from(BufReader::new(file))
}

pub fn read_mtx_from<R: BufRead>(reader: R) -> Result<CsMat<f64>> {
    let mut lines = reader.lines();

    let header = lines
        .next()
        .ok_or_else(|| Error::Parse("empty file".into()))??;
    let fields: Vec<String> = header
        .split_whitespace()
        .map(|s| s.to_ascii_lowercase())
        .collect();
    if fields.len() != 5 || fields[0] != "%%matrixmarket" || fields[1] != "matrix" {
        return Err(Error::Parse(format!(
            "not a Matrix Market header: {header}"
        )));
    }
    if fields[2] != "coordinate" {
        return Err(Error::Parse("only coordinate format is supported".into()));
    }
    let pattern = match fields[3].as_str() {
        "integer" | "real" => false,
        "pattern" => true,
        other => return Err(Error::Parse(format!("unsupported field type: {other}"))),
    };
    if fields[4] != "general" {
        return Err(Error::Parse(format!("unsupported symmetry: {}", fields[4])));
    }

    // Skip comments to reach the size line.
    let size_line = loop {
        let line = lines
            .next()
            .ok_or_else(|| Error::Parse("missing size line".into()))??;
        if !line.starts_with('%') && !line.trim().is_empty() {
            break line;
        }
    };
    let dims: Vec<usize> = size_line
        .split_whitespace()
        .map(|s| {
            s.parse()
                .map_err(|_| Error::Parse(format!("bad size line: {size_line}")))
        })
        .collect::<Result<_>>()?;
    let [n_rows, n_cols, nnz] = dims[..] else {
        return Err(Error::Parse(format!("bad size line: {size_line}")));
    };

    let mut tri = TriMat::with_capacity((n_rows, n_cols), nnz);
    for line in lines {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('%') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let mut next_index = |bound: usize| -> Result<usize> {
            let i: usize = parts
                .next()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| Error::Parse(format!("bad entry: {line}")))?;
            if i == 0 || i > bound {
                return Err(Error::Parse(format!("index out of range: {line}")));
            }
            Ok(i - 1) // Matrix Market is 1-based
        };
        let row = next_index(n_rows)?;
        let col = next_index(n_cols)?;
        let value = if pattern {
            1.0
        } else {
            parts
                .next()
                .and_then(|s| s.parse::<f64>().ok())
                .ok_or_else(|| Error::Parse(format!("bad value: {line}")))?
        };
        tri.add_triplet(row, col, value);
    }
    if tri.nnz() != nnz {
        return Err(Error::Parse(format!(
            "expected {nnz} entries, found {}",
            tri.nnz()
        )));
    }

    Ok(tri.to_csc())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_integer_coordinate_file() {
        let mtx = "%%MatrixMarket matrix coordinate integer general\n\
                   % comment\n\
                   3 2 3\n\
                   1 1 5\n\
                   3 1 2\n\
                   2 2 7\n";
        let m = read_mtx_from(mtx.as_bytes()).unwrap();
        assert_eq!(m.shape(), (3, 2));
        assert!(m.is_csc());
        assert_eq!(m.get(0, 0), Some(&5.0));
        assert_eq!(m.get(2, 0), Some(&2.0));
        assert_eq!(m.get(1, 1), Some(&7.0));
        assert_eq!(m.get(0, 1), None);
    }

    #[test]
    fn rejects_wrong_entry_count() {
        let mtx = "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1.5\n";
        assert!(read_mtx_from(mtx.as_bytes()).is_err());
    }

    #[test]
    fn rejects_out_of_range_index() {
        let mtx = "%%MatrixMarket matrix coordinate real general\n2 2 1\n3 1 1.5\n";
        assert!(read_mtx_from(mtx.as_bytes()).is_err());
    }
}
