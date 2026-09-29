//! Matrix Market (`.mtx`) reading, e.g. 10x Genomics `matrix.mtx(.gz)` output.
//!
//! Parses the coordinate format directly, accepting the `integer` fields
//! that 10x count matrices always use.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::matrix::CscMatrix;
use crate::{Error, Result};

/// Reads a `coordinate` Matrix Market file into a CSC matrix.
///
/// 10x files are gene × cell, so columns are cells. Supports `integer`,
/// `real`, and `pattern` fields with `general` symmetry.
pub fn read_mtx<P: AsRef<Path>>(path: P) -> Result<CscMatrix<usize>> {
    let file = File::open(path)?;
    read_mtx_from(BufReader::new(file))
}

pub fn read_mtx_from<R: BufRead>(reader: R) -> Result<CscMatrix<usize>> {
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

    let mut triplets = Vec::with_capacity(nnz);
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
        triplets.push((row, col, value));
    }
    if triplets.len() != nnz {
        return Err(Error::Parse(format!(
            "expected {nnz} entries, found {}",
            triplets.len()
        )));
    }

    CscMatrix::from_triplets(n_rows, n_cols, triplets)
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
        let v = m.view();
        assert_eq!((v.n_rows(), v.n_cols(), v.nnz()), (3, 2, 3));
        assert_eq!(v.get(0, 0), 5.0);
        assert_eq!(v.get(2, 0), 2.0);
        assert_eq!(v.get(1, 1), 7.0);
        assert_eq!(v.get(0, 1), 0.0);
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
