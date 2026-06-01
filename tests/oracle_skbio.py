#!/usr/bin/env python3
"""scikit-bio PERMANOVA oracle for rsomics-permanova compat tests.

argv: dm.tsv grouping.tsv [permutations] [seed] [precision]

Emits the same key<tab>value result table rsomics-permanova writes. The
pseudo-F statistic is deterministic; the p-value is the scikit-bio
permutation estimate at the given seed (a different RNG from ours).
"""

import sys

from skbio import DistanceMatrix
from skbio.stats.distance import permanova


def main():
    dm_path = sys.argv[1]
    grp_path = sys.argv[2]
    permutations = int(sys.argv[3]) if len(sys.argv) > 3 else 999
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    precision = int(sys.argv[5]) if len(sys.argv) > 5 else 6

    dm = DistanceMatrix.read(dm_path)

    groups = {}
    with open(grp_path) as fh:
        rows = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    ids = set(dm.ids)
    start = 0
    if rows and rows[0].split("\t")[0] not in ids:
        start = 1
    for ln in rows[start:]:
        a, b = ln.split("\t")[:2]
        groups[a] = b
    grouping = [groups[i] for i in dm.ids]

    res = permanova(dm, grouping, permutations=permutations, seed=seed)
    stat = res["test statistic"]
    p = res["p-value"]

    print("method name\tPERMANOVA")
    print("test statistic name\tpseudo-F")
    print(f"sample size\t{int(res['sample size'])}")
    print(f"number of groups\t{int(res['number of groups'])}")
    print(f"test statistic\t{stat:.{precision}f}")
    print(f"number of permutations\t{int(res['number of permutations'])}")
    if permutations == 0:
        print("p-value\tNA")
    else:
        print(f"p-value\t{p:.{precision}f}")


if __name__ == "__main__":
    main()
