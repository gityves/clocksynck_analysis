"""A162 scoring: V162.1-2 and the descriptive tables from `compute_a154.py --kinds P162 --out ../output/a162`.

    python score_a162.py
Writes (../output/a162): target_coverage.csv (design x c_anchor x c_target: coverage with Wilson interval, median bias, boundary share),
anchor_coverage.csv (leave-one-out host, K = 8 designs, first 50 replicates), bias_vs_anchor.csv (median bias by c_anchor, with the A154
cells c_anchor = c_target as reference), summary.csv (V162.1, V162.2).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "output"
A = OUT / "a162"


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    raw = pd.concat([pd.read_csv(f) for f in sorted((A / "raw" / "P162").glob("*.csv"))], ignore_index=True)
    raw = raw[raw.lo.notna()]
    raw["covered"] = (raw.lo <= raw.true_c) & (raw.true_c <= raw.hi)
    raw["bias"] = raw.coupling - raw.true_c
    rows = {}
    for role, name in (("target", "target_coverage.csv"), ("anchor", "anchor_coverage.csv")):
        g = raw[raw.role == role].groupby(["design", "c_anchor", "c_target"])
        t = g.agg(reps=("rep", "nunique"), intervals=("covered", "size"), coverage=("covered", "mean"), median_bias=("bias", "median"),
                  median_estimate=("coupling", "median"), boundary_share=("boundary", "mean")).reset_index()
        t[["wilson_lo", "wilson_hi"]] = [wilson(c * n, n) for c, n in zip(t.coverage, t.intervals)]
        if role == "anchor":
            t = t.drop(columns="c_target").groupby(["design", "c_anchor"]).mean(numeric_only=True).reset_index()
        t.to_csv(A / name, index=False)
        rows[role] = t
    t = rows["target"]
    t["cell_pass"] = t.coverage >= 0.90
    t.to_csv(A / "target_coverage.csv", index=False)
    summ = []
    for ca, name in ((0.8, "V162.1"), (0.6, "V162.2")):
        s = t[t.c_anchor == ca]
        summ.append(dict(criterion=name, c_anchor=ca, cells=len(s), cells_passing=int(s.cell_pass.sum()), min_coverage=s.coverage.min(),
                         median_coverage=s.coverage.median(), result=bool(len(s) == 32 and s.cell_pass.all())))
    pd.DataFrame(summ).to_csv(A / "summary.csv", index=False)
    ref = pd.read_csv(OUT / "v154_3_coverage.csv").rename(columns={"true_c": "c_target"})
    ref = ref.assign(c_anchor=ref.c_target, median_bias=ref.median_estimate - ref.c_target, source="A154")[
        ["design", "c_anchor", "c_target", "coverage", "median_bias", "boundary_share", "source"]]
    cur = t.assign(source="A162")[ref.columns]
    pd.concat([cur, ref]).sort_values(["design", "c_target", "c_anchor"]).to_csv(A / "bias_vs_anchor.csv", index=False)
    print(pd.DataFrame(summ).to_string())


if __name__ == "__main__":
    main()
