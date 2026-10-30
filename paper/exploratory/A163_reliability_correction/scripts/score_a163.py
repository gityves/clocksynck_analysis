"""A163 scoring: coverage, median bias and share > 1 of the raw and corrected coupling, per design x cell x kind, from
`compute_a163.py` raw files; V163.1-4 in validate mode (registered in paper/docs/registered_experiments.md, A163).

    python score_a163.py --mode dev|validate
Writes output/<mode>/: cells.csv (targets; anchors in anchors.csv), summary.csv (per kind; validate: V163.1-4).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "output"


def table(raw: pd.DataFrame, kinds: list) -> pd.DataFrame:
    rows = []
    for (dsg, cell, ca, ct), g in raw.groupby(["design", "cell", "c_anchor", "c_target"]):
        for k in kinds:
            est, lo, hi = g[f"{k}_est"], g[f"{k}_lo"], g[f"{k}_hi"]
            ok = est.notna() & lo.notna() & hi.notna()
            rows.append(dict(design=dsg, cell=cell, c_anchor=ca, c_target=ct, kind=k, n=int(ok.sum()),
                             coverage=float(((lo <= g.true_c) & (g.true_c <= hi))[ok].mean()),
                             median_bias=float((est - g.true_c)[ok].median()),
                             median_ratio=float((est / g.true_c)[ok].median()),
                             share_above_1=float((est[ok] > 1).mean()),
                             median_R=float(g[f"{k}_R"].median()) if f"{k}_R" in g and k != "raw" else np.nan))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", choices=("dev", "validate"), required=True)
    args = ap.parse_args()
    d = OUT / args.mode
    raw = pd.concat([pd.read_csv(f) for f in sorted((d / "raw").glob("*.csv"))], ignore_index=True)
    kinds = [c[:-4] for c in raw.columns if c.endswith("_est")]
    tg = table(raw[raw.role == "target"], kinds)
    tg.to_csv(d / "cells.csv", index=False)
    table(raw[raw.role == "anchor"], kinds).to_csv(d / "anchors.csv", index=False)
    main_cells = tg[~tg.cell.str.startswith("eq")]
    summ = []
    for k in kinds:
        x = main_cells[main_cells.kind == k]
        r = dict(kind=k, cells=len(x), reps_min=int(x.n.min()), cov_min=x.coverage.min(), cov_median=x.coverage.median(),
                 cov_min_ca08=x[x.c_anchor == 0.8].coverage.min(), cov_min_ca06=x[x.c_anchor == 0.6].coverage.min(),
                 absbias_max=x.median_bias.abs().max(), cov_min_ct01=x[x.c_target == 0.1].coverage.min(),
                 cells_cov_ge_090=int((x.coverage >= 0.90).sum()), share_above_1_max=x.share_above_1.max())
        if args.mode == "validate" and k != "raw":
            r.update({"V163.1": bool((x[x.c_anchor == 0.8].coverage >= 0.90).all()),
                      "V163.2": bool((x[x.c_anchor == 0.6].coverage >= 0.90).all()),
                      "V163.3": bool((x.median_bias.abs() <= 0.05).all()),
                      "V163.4": bool((x[x.c_target == 0.1].coverage >= 0.90).all())})
        summ.append(r)
    s = pd.DataFrame(summ)
    s.to_csv(d / "summary.csv", index=False)
    print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
