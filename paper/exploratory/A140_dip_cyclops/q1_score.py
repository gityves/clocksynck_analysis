"""A140 Q1: simulator gate, H140.1-H140.3 with Holm, and the descriptive grid summary (reads output/q1/*/*.csv).

    python q1_score.py gate        simulator gate (arm D, seed 1,500 without core clock genes; with core genes alongside)
    python q1_score.py primary     H140.1-3, exact binomial thresholds, Holm over the three -> output/q1/primary_scores.csv
    python q1_score.py grid        per-cell medians and shares -> output/q1/grid_summary.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent / "output" / "q1"


def read(step: str) -> pd.DataFrame:
    fs = sorted((OUT / step).glob("*.csv"))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def gate():
    d = read("gate")
    rows = []
    for core, g in d.groupby("core"):
        g = g[g.status == "ok"]
        rows.append(dict(core=core, n_ok=len(g), n_total=int((d.core == core).sum()), median_truth_r=g.truth_r_median.median(),
                         share_ge_0_7=(g.truth_r_median >= 0.7).mean(), median_dip_truth_r=g.dip_truth_r.median(),
                         min_truth_r=g.truth_r_median.min()))
    t = pd.DataFrame(rows)
    t["gate_passes"] = t.median_truth_r >= 0.7
    t.to_csv(OUT / "simulator_gate.csv", index=False)
    print(t.to_string())


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for rank, i in enumerate(o):
        run = max(run, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = run
    return adj


def primary():
    d = read("primary")
    b = d[(d.arm == "B") & (d.status == "ok")]
    n_b = len(b)
    x1 = int((b.cc_median >= 0.5).sum())
    x2 = int((b.n_rhy_min >= 500).sum())
    p1 = stats.binom.sf(x1 - 1, 50, 0.5)
    p2 = stats.binom.sf(x2 - 1, 50, 0.5)
    gates = d.drop_duplicates(["arm", "rep"]).groupby("arm").gate_pass.agg(["sum", "count"])
    xB, xC, xA = (int(gates.loc[a, "sum"]) for a in "BCA")
    pB, pC = stats.binom.cdf(xB, 50, 0.10), stats.binom.cdf(xC, 50, 0.10)  # P(X <= x | 0.10): small = specific
    pA = stats.binom.cdf(xA, 50, 0.80)  # P(X <= x | 0.80): >= 0.05 = not rejected downward
    p3 = max(pB, pC)
    adj = holm([p1, p2, p3])
    rows = [dict(h="H140.1", count=x1, n=50, n_fits_ok=n_b, threshold=">= 32 of 50", raw_p=p1, holm_p=adj[0],
                 passes=bool(x1 >= 32 and adj[0] < 0.05)),
            dict(h="H140.2", count=x2, n=50, n_fits_ok=n_b, threshold=">= 32 of 50", raw_p=p2, holm_p=adj[1],
                 passes=bool(x2 >= 32 and adj[1] < 0.05)),
            dict(h="H140.3", count=f"B {xB}, C {xC}, A {xA}", n=50, n_fits_ok=n_b,
                 threshold="B <= 1, C <= 1, A >= 35", raw_p=p3, holm_p=adj[2],
                 passes=bool(xB <= 1 and xC <= 1 and xA >= 35 and adj[2] < 0.05))]
    t = pd.DataFrame(rows)
    t["extra"] = ["", "", f"p_B {pB:.4g}, p_C {pC:.4g}, A: P(X<=x|0.8) {pA:.4g}"]
    t.to_csv(OUT / "primary_scores.csv", index=False)
    desc = b[["cc_median", "cc_cd4_cd8", "cc_cd4_myel", "cc_cd8_myel", "n_rhy_cd4", "n_rhy_cd8", "n_rhy_myel", "n_rhy_min",
              "dip_R", "dip_cover", "acro_R_mean", "clock_order_R", "n_bins"]].describe().T
    desc.to_csv(OUT / "primary_arm_B_descriptive.csv")
    print(t.to_string())
    print(desc.to_string())


def grid():
    fs = [f for sub in ("grid", "primary", "gate") for f in (OUT / sub).glob("*.csv")]
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.n_seed > 0]
    ok = d[d.status == "ok"]
    key = ["arm", "K", "n_seed", "core"]
    s = ok.groupby(key).agg(n=("rep", "size"), cc_median=("cc_median", "median"), share_cc_ge_0_5=("cc_median", lambda x: (x >= 0.5).mean()),
                            n_rhy_min_median=("n_rhy_min", "median"), share_rhy_ge_500=("n_rhy_min", lambda x: (x >= 500).mean()),
                            dip_cover_median=("dip_cover", "median"), dip_R_median=("dip_R", "median"),
                            acro_R_median=("acro_R_mean", "median"), clock_order_R_median=("clock_order_R", "median"),
                            truth_r_median=("truth_r_median", "median"), gate_pass_share=("gate_pass", "mean")).reset_index()
    s.to_csv(OUT / "grid_summary.csv", index=False)
    print(s.to_string())


if __name__ == "__main__":
    {"gate": gate, "primary": primary, "grid": grid}[sys.argv[1]]()
