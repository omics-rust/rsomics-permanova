# rsomics-permanova

PERMANOVA — Permutational Multivariate Analysis of Variance (Anderson 2001).
Tests whether two or more groups of samples differ, given a symmetric distance
matrix and a categorical grouping. Reports the pseudo-F statistic and a
permutation p-value.

```
rsomics-permanova dm.tsv -g groups.tsv --seed 42
rsomics-permanova dm.tsv -g groups.tsv -p 0     # statistic only
```

## Input

- **Distance matrix** (`dm.tsv`): square, symmetric, lsmat layout — a blank
  top-left cell, sample IDs across the header, each data row prefixed by its
  ID. The same format `rsomics-pcoa` / `rsomics-beta-diversity` emit.
- **Grouping** (`-g groups.tsv`): `id<tab>group` per line; an optional header
  is detected and skipped. Every matrix ID must be present; extra IDs are
  ignored.

## Output

A `key<tab>value` table: method name, test statistic name (pseudo-F), sample
size, number of groups, the pseudo-F statistic, number of permutations, and the
p-value.

## Method

The pseudo-F statistic follows Anderson (2001):

```
F = (SS_between / (g-1)) / (SS_within / (n-g))
```

with `SS_within` computed from the within-group squared-distance sums divided
by group size, and `SS_total = Σ d² / n` (each pair once). `SS_between =
SS_total − SS_within`. The p-value is `(1 + #{F_perm ≥ F}) / (1 + permutations)`,
permuting the grouping labels.

The **pseudo-F statistic is value-exact** versus scikit-bio (to 1e-9). The
**p-value is a seeded permutation estimate**: scikit-bio permutes with NumPy's
`Generator.permutation` (PCG64); this crate uses its own seeded PCG64
Fisher-Yates, so for a given seed the two p-values are independent Monte-Carlo
estimates of the same true p (they agree within Monte-Carlo tolerance, not
bit-for-bit). Increasing `-p` tightens both. The estimate is reproducible for a
fixed `--seed` and independent of thread count.

## Origin

This crate is an independent Rust reimplementation of the PERMANOVA operation
provided by `scikit-bio` (`skbio.stats.distance.permanova`), based on:

- The published method: Anderson, M. J. "A new method for non-parametric
  multivariate analysis of variance." *Austral Ecology* 26.1 (2001): 32–46.
  doi:10.1111/j.1442-9993.2001.01070.pp.x
- The scikit-bio source (BSD-3-Clause), read and cited for the exact
  sum-of-squares decomposition (within-group squared-distance sums normalised
  by group size; total SS halved over the redundant matrix).

scikit-bio is BSD-3-Clause and was read and cited. Test fixtures are
independently generated distance matrices with planted group structure.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio <https://scikit-bio.org> (BSD-3-Clause).

## Compatibility & performance

`tests/compat.rs` compares this binary against committed scikit-bio golden
outputs (always runs) and, when scikit-bio is importable, against the live
oracle (`tests/oracle_skbio.py`). The pseudo-F statistic is asserted exact to
1e-9; the permutation p-value is asserted within Monte-Carlo tolerance.
