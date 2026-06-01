use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A labeled symmetric distance matrix in redundant (square) form.
pub struct DistanceMatrix {
    pub ids: Vec<String>,
    /// Row-major `n*n` distances.
    pub data: Vec<f64>,
}

impl DistanceMatrix {
    pub fn n(&self) -> usize {
        self.ids.len()
    }

    #[inline]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.ids.len() + j]
    }

    /// Parse the scikit-bio / lsmat layout: an empty top-left cell, column IDs
    /// across the header, each data row prefixed by its ID.
    ///
    /// # Errors
    /// Errors on a ragged matrix, a row/column label mismatch, or a non-numeric cell.
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<DistanceMatrix> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty distance matrix".into())),
            }
        };
        let col_ids: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        if col_ids.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "distance-matrix header has no column IDs".into(),
            ));
        }
        let n = col_ids.len();

        let mut ids = Vec::with_capacity(n);
        let mut data = Vec::with_capacity(n * n);
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split(delim);
            let row_id = fields
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("empty matrix row".into()))?
                .trim();
            ids.push(row_id.to_string());
            let mut cells = 0usize;
            for field in fields {
                let v: f64 = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row '{row_id}': '{}' is not a number",
                        field.trim()
                    ))
                })?;
                data.push(v);
                cells += 1;
            }
            if cells != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "row '{row_id}' has {cells} cells, header has {n} columns"
                )));
            }
        }
        if ids.len() != n {
            return Err(RsomicsError::InvalidInput(format!(
                "matrix has {} rows but {n} columns",
                ids.len()
            )));
        }
        if ids != col_ids {
            return Err(RsomicsError::InvalidInput(
                "row IDs do not match column IDs in the same order".into(),
            ));
        }
        Ok(DistanceMatrix { ids, data })
    }
}
