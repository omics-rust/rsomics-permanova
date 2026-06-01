use std::io::{BufRead, Write};

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

mod dm;
mod grouping;
mod rng;

pub use dm::DistanceMatrix;
pub use grouping::{Grouping, parse as parse_grouping};

use rng::Pcg64;

pub struct PermanovaResult {
    pub sample_size: usize,
    pub num_groups: usize,
    pub statistic: f64,
    pub permutations: usize,
    pub p_value: Option<f64>,
}

/// Squared distances kept once: the permutation loop only ever needs `d²`, so
/// pay the square up front and let the inner loop be a masked multiply-add.
struct SquaredMatrix {
    n: usize,
    sq: Vec<f64>,
    total: f64,
}

impl SquaredMatrix {
    fn build(dm: &DistanceMatrix) -> SquaredMatrix {
        let n = dm.n();
        let sq: Vec<f64> = dm.data.iter().map(|&d| d * d).collect();
        let total = (sq.iter().sum::<f64>() / n as f64) / 2.0;
        SquaredMatrix { n, sq, total }
    }

    /// `SS_within` for one assignment. The membership test is folded into a
    /// branchless `0.0/1.0` multiplier so the row scan vectorises; the `1/size`
    /// divisor is constant within a row, applied once per row.
    fn s_within(&self, codes: &[u32], inv_size: &[f64]) -> f64 {
        let n = self.n;
        let mut total = 0.0;
        for i in 0..n - 1 {
            let gi = codes[i];
            let row = &self.sq[i * n + i + 1..i * n + n];
            let other = &codes[i + 1..];
            let mut local = 0.0;
            for (&gj, &d2) in other.iter().zip(row) {
                local += d2 * f64::from(u8::from(gj == gi));
            }
            total += local * inv_size[gi as usize];
        }
        total
    }
}

fn pseudo_f(s_t: f64, s_w: f64, n: usize, g: usize) -> f64 {
    let s_a = s_t - s_w;
    (s_a / (g - 1) as f64) / (s_w / (n - g) as f64)
}

/// Run PERMANOVA. The statistic is exact; the p-value (when `permutations > 0`)
/// is a seeded permutation estimate.
///
/// # Errors
/// Errors only if `permutations` overflows the worker arithmetic; inputs are
/// validated upstream.
pub fn permanova(
    dm: &DistanceMatrix,
    grouping: &Grouping,
    permutations: usize,
    seed: u64,
    threads: usize,
) -> PermanovaResult {
    let n = dm.n();
    let g = grouping.num_groups();
    let sq = SquaredMatrix::build(dm);
    let s_t = sq.total;

    let codes: Vec<u32> = grouping.codes.iter().map(|&c| c as u32).collect();
    let inv_size: Vec<f64> = grouping
        .group_sizes
        .iter()
        .map(|&s| 1.0 / s as f64)
        .collect();

    let s_w = sq.s_within(&codes, &inv_size);
    let stat = pseudo_f(s_t, s_w, n, g);

    let p_value = if permutations == 0 {
        None
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .expect("rayon pool");
        let ge = pool.install(|| {
            (0..permutations)
                .into_par_iter()
                .map(|i| {
                    let mut rng =
                        Pcg64::seed(seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
                    let mut perm = codes.clone();
                    rng.shuffle(&mut perm);
                    let perm_sw = sq.s_within(&perm, &inv_size);
                    let f = pseudo_f(s_t, perm_sw, n, g);
                    usize::from(f >= stat)
                })
                .sum::<usize>()
        });
        Some((ge + 1) as f64 / (permutations + 1) as f64)
    };

    PermanovaResult {
        sample_size: n,
        num_groups: g,
        statistic: stat,
        permutations,
        p_value,
    }
}

pub struct Config {
    pub permutations: usize,
    pub seed: u64,
    pub threads: usize,
    pub delim: char,
    pub precision: usize,
}

/// Parse a distance matrix + grouping, run PERMANOVA, and write the result table.
///
/// # Errors
/// Propagates parse errors.
pub fn run<R: BufRead, G: BufRead, W: Write>(
    dm_reader: R,
    grouping_reader: G,
    mut out: W,
    cfg: &Config,
) -> Result<()> {
    let dm = DistanceMatrix::parse(dm_reader, cfg.delim)?;
    if dm.n() < 2 {
        return Err(RsomicsError::InvalidInput(
            "distance matrix needs at least 2 samples".into(),
        ));
    }
    let grouping = parse_grouping(grouping_reader, &dm.ids, cfg.delim)?;
    let res = permanova(&dm, &grouping, cfg.permutations, cfg.seed, cfg.threads);
    write_result(&mut out, &res, cfg.precision)
}

fn write_result<W: Write>(out: &mut W, res: &PermanovaResult, precision: usize) -> Result<()> {
    writeln!(out, "method name\tPERMANOVA").map_err(RsomicsError::Io)?;
    writeln!(out, "test statistic name\tpseudo-F").map_err(RsomicsError::Io)?;
    writeln!(out, "sample size\t{}", res.sample_size).map_err(RsomicsError::Io)?;
    writeln!(out, "number of groups\t{}", res.num_groups).map_err(RsomicsError::Io)?;
    writeln!(out, "test statistic\t{:.*}", precision, res.statistic).map_err(RsomicsError::Io)?;
    writeln!(out, "number of permutations\t{}", res.permutations).map_err(RsomicsError::Io)?;
    match res.p_value {
        Some(p) => writeln!(out, "p-value\t{:.*}", precision, p),
        None => writeln!(out, "p-value\tNA"),
    }
    .map_err(RsomicsError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dm() -> DistanceMatrix {
        let text = "\tA\tB\tC\tD\nA\t0\t1\t4\t5\nB\t1\t0\t3\t6\nC\t4\t3\t0\t2\nD\t5\t6\t2\t0\n";
        DistanceMatrix::parse(text.as_bytes(), '\t').unwrap()
    }

    #[test]
    fn parses_square_dm() {
        let m = dm();
        assert_eq!(m.ids, ["A", "B", "C", "D"]);
        assert_eq!(m.at(0, 2), 4.0);
        assert_eq!(m.at(3, 1), 6.0);
    }

    #[test]
    fn factor_encoding_sorted() {
        let m = dm();
        let g = "A\tg2\nB\tg2\nC\tg1\nD\tg1\n";
        let grp = parse_grouping(g.as_bytes(), &m.ids, '\t').unwrap();
        // labels sorted: g1->0, g2->1
        assert_eq!(grp.labels, ["g1", "g2"]);
        assert_eq!(grp.codes, [1, 1, 0, 0]);
        assert_eq!(grp.group_sizes, [2, 2]);
    }

    #[test]
    fn statistic_is_deterministic() {
        let m = dm();
        let g = "A\tx\nB\tx\nC\ty\nD\ty\n";
        let grp = parse_grouping(g.as_bytes(), &m.ids, '\t').unwrap();
        let r1 = permanova(&m, &grp, 0, 42, 1);
        let r2 = permanova(&m, &grp, 0, 42, 1);
        assert_eq!(r1.statistic.to_bits(), r2.statistic.to_bits());
    }

    #[test]
    fn p_value_reproducible_with_seed() {
        let m = dm();
        let g = "A\tx\nB\tx\nC\ty\nD\ty\n";
        let grp = parse_grouping(g.as_bytes(), &m.ids, '\t').unwrap();
        let r1 = permanova(&m, &grp, 99, 7, 4);
        let r2 = permanova(&m, &grp, 99, 7, 1);
        assert_eq!(
            r1.p_value, r2.p_value,
            "p-value must not depend on thread count"
        );
    }

    #[test]
    fn single_group_errors() {
        let m = dm();
        let g = "A\tx\nB\tx\nC\tx\nD\tx\n";
        assert!(parse_grouping(g.as_bytes(), &m.ids, '\t').is_err());
    }
}
