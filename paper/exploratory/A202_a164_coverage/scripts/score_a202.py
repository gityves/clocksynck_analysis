"""A202 scoring: V202.1 (coverage of the registered A164 interval, arm B) and V202.2 (equivalence size at +/-0.15,
arm P), with the descriptive tables, from `compute_a202.py`.

    python score_a202.py [--subsample]
Writes (../output): reread.csv (registered re-reading of A164 when a criterion fails), arm_p.csv (design x cell: realised truth, bias, SD_sim against the A164 SE and null SD, coverage
with each half-width, TOST rate), arm_b.csv (design x cell: coverage of the replicate's own interval, SE_boot
distribution), per_cohort.csv, attenuation.csv (A164 intervals divided by the A162 range), summary.csv (V202.1-2).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
A164 = HERE.parents[1] / "A164_gtex_metabolic_matched/output"
Z90, BOUND = 1.645, 0.15
COVER_MIN, TOST_MAX = 0.85, 0.075
ATTENUATION = (0.80, 0.90)  # A162 range of the attenuation of c at realistic anchor coupling


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def a164():
    """Pooled module weights, bootstrap SE and null SD per design (a164_power.csv), and the observed intervals."""
    pw = pd.read_csv(A164 / "a164_power.csv")
    pw = pw[(pw.statistic == "module") & (pw.set == "pooled")].set_index("design")
    w = {d: {k: float(v) for k, v in (x.split(":") for x in pw.loc[d, "weights"].split(";"))} for d in pw.index}
    obs = pd.read_csv(A164 / "a164_tests.csv")
    obs = obs[(obs.statistic == "module") & (obs.set == "pooled")].set_index("design")
    return w, pw.boot_se.to_dict(), pw.null_sd.to_dict(), obs


def pooled(raw: pd.DataFrame, w: dict) -> pd.DataFrame:
    """Per design x cell x replicate: the fixed-weight pooled module estimate, truth and (arm B) SE."""
    rows = []
    for (dsg, cell, delta, rep), g in raw.groupby(["design", "cell", "delta", "rep"]):
        g = g.set_index("cohort")
        ww = w[dsg]
        if not all(c in g.index and np.isfinite(g.loc[c, "D"]) for c in ww):
            continue
        r = dict(design=dsg, cell=cell, delta=delta, rep=rep,
                 D=sum(ww[c] * g.loc[c, "D"] for c in ww), truth=sum(ww[c] * g.loc[c, "truth"] for c in ww))
        if "se_boot" in g:
            r["se"] = float(np.sqrt(sum(ww[c] ** 2 * g.loc[c, "se_boot"] ** 2 for c in ww)))
        rows.append(r)
    return pd.DataFrame(rows)


def load(base: Path, arm: str) -> pd.DataFrame:
    files = sorted((base / "raw" / arm).glob("*.csv"))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--subsample", action="store_true")
    args = ap.parse_args()
    base = OUT / "subsample" if args.subsample else OUT
    w, se164, null164, obs = a164()
    summ = []

    rawp = load(base, "P")
    if len(rawp):
        P = pooled(rawp, w)
        rows = []
        for (dsg, delta), g in P.groupby(["design", "delta"]):
            err = g.D - g.truth
            sd_sim = float(g.D.std(ddof=1))
            r = dict(design=dsg, delta=delta, n=len(g), truth_mean=g.truth.mean(), D_median=g.D.median(),
                     bias_median=float(err.median()), sd_sim=sd_sim, se_a164=se164[dsg], null_sd_a164=null164[dsg],
                     sd_sim_over_se=sd_sim / se164[dsg], sd_sim_over_null=sd_sim / null164[dsg])
            for name, h in (("se_a164", se164[dsg]), ("sd_sim", sd_sim), ("null_sd", null164[dsg])):
                k = int((err.abs() <= Z90 * h).sum())
                r[f"cover_{name}"] = k / len(g)
            inside = ((g.D - Z90 * se164[dsg]) > -BOUND) & ((g.D + Z90 * se164[dsg]) < BOUND)
            r["tost_rate"] = float(inside.mean())
            r["tost_lo"], r["tost_hi"] = wilson(int(inside.sum()), len(g))
            rows.append(r)
        AP = pd.DataFrame(rows)
        AP.to_csv(base / "arm_p.csv", index=False)
        rawp.groupby(["design", "delta", "cohort"]).agg(
            n=("D", "size"), D_median=("D", "median"), truth_mean=("truth", "mean"), sd=("D", "std"),
            n_tissues=("n_tissues", "mean"), thin_failed=("n_thin_failed", "sum")).reset_index().to_csv(
            base / "per_cohort.csv", index=False)
        for dsg in w:
            b = AP[(AP.design == dsg) & (AP.delta.abs() == BOUND)]
            if len(b):
                summ.append(dict(criterion="V202.2", design=dsg, value=float(b.tost_rate.max()), threshold=TOST_MAX,
                                 passes=bool(len(b) == 2 and (b.tost_rate <= TOST_MAX).all()),
                                 note="max TOST rate over delta -0.15 and +0.15; realised truths "
                                      + ", ".join(f"{x:+.3f}" for x in b.truth_mean)))

    rawb = load(base, "B")
    if len(rawb):
        B = pooled(rawb, w)
        rows = []
        for (dsg, delta), g in B.groupby(["design", "delta"]):
            cov = (g.D - g.truth).abs() <= Z90 * g.se
            lo, hi = wilson(int(cov.sum()), len(g))
            rows.append(dict(design=dsg, delta=delta, n=len(g), truth_mean=g.truth.mean(), coverage=float(cov.mean()),
                             cover_lo=lo, cover_hi=hi, se_median=g.se.median(), se_q05=g.se.quantile(.05),
                             se_q95=g.se.quantile(.95), sd_D=g.D.std(ddof=1), se_a164=se164[dsg]))
        AB = pd.DataFrame(rows)
        AB.to_csv(base / "arm_b.csv", index=False)
        for dsg in w:
            b = AB[AB.design == dsg]
            if len(b):
                summ.append(dict(criterion="V202.1", design=dsg, value=float(b.coverage.min()), threshold=COVER_MIN,
                                 passes=bool(len(b) == 2 and (b.coverage >= COVER_MIN).all()),
                                 note="min coverage over delta 0 and -0.15"))

    if len(rawp) and len(summ) and not all(r["passes"] for r in summ):
        # registered re-reading after a failed V202.1 or V202.2: half-width 1.645 x max(SD_sim over the arm P cells,
        # A164 null SD); bias correction only if |median bias| > 0.02 in some cell
        rr = []
        for dsg, r in obs.iterrows():
            g = AP[AP.design == dsg]
            sd = max(float(g.sd_sim.max()), null164[dsg])
            bias = float(g.bias_median.abs().max())
            D = r.D - (float(g.bias_median.median()) if bias > 0.02 else 0.0)
            rr.append(dict(design=dsg, D=D, half_width=Z90 * sd, lo90=D - Z90 * sd, hi90=D + Z90 * sd,
                           sd_used=sd, max_abs_bias=bias, bias_corrected=bias > 0.02,
                           inside=bool(D - Z90 * sd > -BOUND and D + Z90 * sd < BOUND),
                           a164_lo90=r.lo90, a164_hi90=r.hi90))
        pd.DataFrame(rr).to_csv(base / "reread.csv", index=False)

    att = []
    for dsg, r in obs.iterrows():
        for a in ATTENUATION:
            lo, hi = r.lo90 / a, r.hi90 / a
            att.append(dict(design=dsg, attenuation=a, lo90=lo, hi90=hi, inside=bool(lo > -BOUND and hi < BOUND)))
    pd.DataFrame(att).to_csv(base / "attenuation.csv", index=False)
    S = pd.DataFrame(summ)
    S.to_csv(base / "summary.csv", index=False)
    pd.set_option("display.width", 220)
    for f in ("arm_p.csv", "arm_b.csv", "attenuation.csv", "summary.csv", "reread.csv"):
        if (base / f).exists():
            print(f, "\n", pd.read_csv(base / f).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
