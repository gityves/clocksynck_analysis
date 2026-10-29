"""A154 scoring: V154.1-3 and the descriptive tables from `compute_a154.py` raw files.

    python score_a154.py [--out DIR]
Writes (DIR): v154_1_n0.csv, v154_2_null.csv, v154_3_coverage.csv, lower_bound_coverage.csv, n1_map.csv, p_sensitivity.csv,
summary.csv (one line per criterion).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent.parent / "output"


def read(out: Path, kind: str) -> pd.DataFrame:
    fs = sorted((out / "raw" / kind).glob("*.csv"))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    out = ap.parse_args().out
    summ = []
    # V154.1
    n0 = read(out, "N0")
    if len(n0):
        r = n0.groupby("design").passes.agg(["sum", "count"]).rename(columns={"sum": "pass", "count": "reps"})
        r["rate"] = r["pass"] / r.reps
        r["binom_p_greater_than_0.05"] = [stats.binomtest(int(k), int(n), 0.05, alternative="greater").pvalue for k, n in zip(r["pass"], r.reps)]
        r["rate_le_0.075"] = r.rate <= 0.075
        r["binom_not_rejected_p>=0.01"] = r["binom_p_greater_than_0.05"] >= 0.01
        r["pass"] = r["pass"].astype(int)
        r["V154.1_design_pass"] = r["rate_le_0.075"] & r["binom_not_rejected_p>=0.01"]
        r.to_csv(out / "v154_1_n0.csv")
        summ.append(dict(criterion="V154.1", value=f"max N0 pass rate {r.rate.max():.3f}", designs=len(r),
                         designs_passing=int(r["V154.1_design_pass"].sum()), result=bool(r["V154.1_design_pass"].all() and len(r) == 8)))
    # V154.2
    nu = read(out, "GNULL")
    if len(nu):
        nu = nu[nu.p.notna()]
        rows = []
        for dsg, g in nu.groupby("design"):
            k = int((g.p < 0.05).sum())
            lo, hi = wilson(k, len(g))
            # one target per replicate (the first target) removes the within-replicate dependence
            first = g[g.compartment == g.compartment.min()]
            rows.append(dict(design=dsg, reps=g.rep.nunique(), tests=len(g), **{"share_p_lt_0.05": k / len(g)}, wilson_lo=lo, wilson_hi=hi,
                             ks_p=stats.kstest(g.p, "uniform").pvalue, share_first_target=(first.p < 0.05).mean(),
                             ks_p_first_target=stats.kstest(first.p, "uniform").pvalue))
        r = pd.DataFrame(rows).set_index("design")
        r["V154.2_design_pass"] = r["share_p_lt_0.05"] <= 0.075
        r.to_csv(out / "v154_2_null.csv")
        summ.append(dict(criterion="V154.2", value=f"max share p<0.05 {r['share_p_lt_0.05'].max():.3f}", designs=len(r),
                         designs_passing=int(r["V154.2_design_pass"].sum()), result=bool(r["V154.2_design_pass"].all() and len(r) == 8)))
    # V154.3 and lower bound
    cv = read(out, "COVER")
    if len(cv):
        cv["covered"] = (cv.lo <= cv.true_c) & (cv.true_c <= cv.hi)
        rows = []
        for (dsg, c), g in cv.groupby(["design", "true_c"]):
            lo, hi = wilson(int(g.covered.sum()), len(g))
            nb = g[~g.boundary]
            rows.append(dict(design=dsg, true_c=c, reps=g.rep.nunique(), intervals=len(g), coverage=g.covered.mean(), wilson_lo=lo,
                             wilson_hi=hi, boundary_share=g.boundary.mean(), coverage_off_boundary=nb.covered.mean() if len(nb) else np.nan,
                             median_estimate=g.coupling.median()))
        r = pd.DataFrame(rows)
        r["V154.3_cell_pass"] = r.coverage >= 0.90
        r.to_csv(out / "v154_3_coverage.csv", index=False)
        summ.append(dict(criterion="V154.3", value=f"min coverage {r.coverage.min():.3f}", designs=len(r),
                         designs_passing=int(r["V154.3_cell_pass"].sum()), result=bool(r["V154.3_cell_pass"].all() and len(r) == 24)))
        b = cv[cv.boundary & cv.c_lo.notna()].copy()
        if len(b):
            b["bound_holds"] = b.true_c >= b.c_lo - 1e-9
            lb = b.groupby(["design", "true_c"]).bound_holds.agg(boundary_fits="count", coverage="mean").reset_index()
            lb.to_csv(out / "lower_bound_coverage.csv", index=False)
            pooled = b.groupby("true_c").bound_holds.agg(["count", "mean"])
            pooled.to_csv(out / "lower_bound_coverage_by_c.csv")
    # descriptive: N1 map and P sensitivity
    n1 = read(out, "N1")
    if len(n1):
        r = n1.groupby(["design", "cell"]).passes.agg(["sum", "count"]).reset_index()
        r["rate"] = r["sum"] / r["count"]
        r[["wilson_lo", "wilson_hi"]] = [wilson(k, n) for k, n in zip(r["sum"], r["count"])]
        r = r.rename(columns={"sum": "pass", "count": "reps"})
        r.to_csv(out / "n1_map.csv", index=False)
    pp = read(out, "P")
    if len(pp):
        r = pp.groupby(["design", "cell"]).passes.agg(["sum", "count"]).reset_index()
        r["rate"] = r["sum"] / r["count"]
        r = r.rename(columns={"sum": "pass", "count": "reps"})
        r.to_csv(out / "p_sensitivity.csv", index=False)
    pd.DataFrame(summ).to_csv(out / "summary.csv", index=False)
    print(pd.DataFrame(summ).to_string())


if __name__ == "__main__":
    main()
