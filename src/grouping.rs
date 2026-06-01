use std::collections::HashMap;
use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// Grouping factor aligned to a distance matrix's ID order.
pub struct Grouping {
    /// Integer factor per sample, indexing the sorted-unique label set.
    pub codes: Vec<usize>,
    pub group_sizes: Vec<usize>,
    /// Sorted-unique labels, parallel to factor codes.
    pub labels: Vec<String>,
}

impl Grouping {
    pub fn num_groups(&self) -> usize {
        self.labels.len()
    }
}

/// Read an `id<delim>group` table (an optional header is detected and skipped),
/// align it to `ids`, then factor-encode the labels the way
/// `numpy.unique(return_inverse=True)` does: lexicographically-sorted unique
/// labels mapped to 0..g.
///
/// # Errors
/// Errors when an ID is missing a group, when all labels are equal, or when
/// every sample is its own group.
pub fn parse<R: BufRead>(reader: R, ids: &[String], delim: char) -> Result<Grouping> {
    let mut map: HashMap<String, String> = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split(delim);
        let id = it.next().unwrap_or("").trim();
        let group = it
            .next()
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "grouping line {} has no group column",
                    i + 1
                ))
            })?
            .trim();
        if i == 0 && !ids.iter().any(|x| x == id) {
            // a header row whose first field is not a sample ID
            continue;
        }
        map.insert(id.to_string(), group.to_string());
    }

    let raw: Vec<String> = ids
        .iter()
        .map(|id| {
            map.get(id)
                .cloned()
                .ok_or_else(|| RsomicsError::InvalidInput(format!("no group for sample '{id}'")))
        })
        .collect::<Result<_>>()?;

    let mut labels: Vec<String> = raw.clone();
    labels.sort();
    labels.dedup();
    let index: HashMap<&str, usize> = labels
        .iter()
        .enumerate()
        .map(|(i, l)| (l.as_str(), i))
        .collect();
    let codes: Vec<usize> = raw.iter().map(|l| index[l.as_str()]).collect();

    if labels.len() == 1 {
        return Err(RsomicsError::InvalidInput(
            "all samples are in a single group — PERMANOVA needs ≥2 groups".into(),
        ));
    }
    if labels.len() == codes.len() {
        return Err(RsomicsError::InvalidInput(
            "every sample is its own group — no within-group distances".into(),
        ));
    }

    let mut group_sizes = vec![0usize; labels.len()];
    for &c in &codes {
        group_sizes[c] += 1;
    }
    Ok(Grouping {
        codes,
        group_sizes,
        labels,
    })
}
