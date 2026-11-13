"""A258: bias of two-group contrasts at unequal host spread, pooled against arm-specific host prior.

Registered in paper/docs/registered_experiments.md (A258). Two arms of simulated donors (`clocksync.simulate`), three anchors
and one target (compartments of the review's feasibility look), identical loadings and noise in both arms (truth d* = 0);
reference window 6-18 h, case window 6-18 / 8-16 / 10-14 / 11-13 h; anchor noise psi 0.002 / 0.05 / 0.15; donors per arm
40 / 80 / 200. A power arm lowers the target's coupling by 0.15 in the case arm (80 donors). Each data set is read by
`Dataset.contrast` with host_prior "pooled" and "arm".

    python simulate_a258.py --part bias|size|power [--seeds 100] [--jobs 8] [--smoke]

Outputs (output/): bias_runs.csv (one row per seed, cell, estimator, compartment), size_runs.csv (permutation p of d*,
200 permutations), power_runs.csv, and the summaries summary_bias.csv, summary_size.csv, summary_power.csv.

Choices not fixed by the registration
- The target's own noise is psi 0.05 (as the review's look); in the power arm its psi in the case arm is raised so that its
  coupling at the reference spread is 0.15 lower (`coupling = |z|^2 V / (|z|^2 V + 2 psi)` with V of the 6-18 h window).
- The size study uses the 12 cells at 80 donors, 30 seeds each; p is two-sided on d*, blocks "none".
- Anchors are also read (their leave-one-out host): rows for "fibroblast" are kept as the anchor case.
"""

from __future__ import annotations

import argparse
import itertools
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402

SEED = 258
WF = cs.gtex()
ANCH = ["fibroblast", "endothelial", "pericyte"]
Z = {"fibroblast": 1.0, "endothelial": 0.8 * np.exp(1j * 2 / HOURS), "pericyte": 1.25 * np.exp(-1j * 1.5 / HOURS),
     "T cell": 0.6 * np.exp(-1j * 3 / HOURS)}
PSI_T = 0.05
WINDOWS = [(6, 18), (8, 16), (10, 14), (11, 13)]
PSIS = [0.002, 0.05, 0.15]
NS = [40, 80, 200]


def cohort(n, seed, window, psi_a, prefix, psi_t=PSI_T):
    comps = {k: cs.Compartment(Z[k], psi=psi_a, mesor_sd=0.0) for k in ANCH}
    comps["T cell"] = cs.Compartment(Z["T cell"], psi=psi_t, mesor_sd=0.0)
    df, _ = cs.simulate(WF, comps, n=n, seed=seed, window_h=window, jitter_h=1.0, tau=0.0)
    df["resection"] = prefix + df.resection.str[1:]
    return df.rename(columns={"resection": "donor"})


def psi_for_drop(drop=0.15):
    """Target psi in the case arm that lowers its coupling by `drop` at the 6-18 h spread: z2 V of the target from a fit of
    one large reference cohort (1,500 donors), then psi = z2 V (1 - c1) / (2 c1) with c1 = c0 - drop."""
    t = cs.Dataset(cohort(1500, [SEED, 7], (6, 18), 0.002, "r"), anchors=ANCH).fit(n_boot=0).table.loc["T cell"]
    s = float(t.z2_debiased * t.V)
    c0 = s / (s + 2 * PSI_T)
    c1 = c0 - drop
    return s * (1 - c1) / (2 * c1), c0


def one(task):
    part, psi_a, wb, n, s, psi_case = task
    limit_threads(1)
    warnings.simplefilter("ignore")
    a = cohort(n, [SEED, 1, s, n, int(psi_a * 1000)], (6, 18), psi_a, "a")
    b = cohort(n, [SEED, 2, s, n, int(psi_a * 1000), wb[0]], wb, psi_a, "b", psi_t=psi_case)
    g = pd.Series(["control"] * n + ["case"] * n, index=sorted(a.donor.unique()) + sorted(b.donor.unique()))
    ds = cs.Dataset(pd.concat([a, b], ignore_index=True), anchors=ANCH)
    rows = []
    for prior in ("pooled", "arm"):
        n_perm = 200 if part == "size" or part == "power" else 0
        c = ds.contrast(g, n_perm=n_perm, n_boot=0, blocks="none", host_prior=prior, seed=s)
        t = c.table
        for k in ("T cell", "fibroblast"):
            r = t.loc[k]
            row = dict(part=part, psi_anchor=psi_a, case_window=f"{wb[0]}-{wb[1]}", n=n, seed=s, prior=prior, comp=k,
                       d=r.d, d_star=r.d_star, log_V=r.log_V_ratio, log_z2=r.log_z2_ratio, c_control=r.c_control, c_case=r.c_case)
            if n_perm:
                row["p_star"] = r.get("p_d_star", r.get("p_star", np.nan))
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["bias", "size", "power"], required=True)
    ap.add_argument("--seeds", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / ("smoke" if a.smoke else "")
    out.mkdir(parents=True, exist_ok=True)
    if a.part == "bias":
        tasks = [("bias", p, w, n, s, PSI_T) for p, w, n in itertools.product(PSIS, WINDOWS, NS) for s in range(a.seeds)]
    elif a.part == "size":
        tasks = [("size", p, w, 80, s, PSI_T) for p, w in itertools.product(PSIS, WINDOWS) for s in range(a.seeds)]
    else:
        psi_case, c0 = psi_for_drop()
        print("power arm: target coupling at 6-18 h", round(c0, 3), "case psi", round(psi_case, 4), flush=True)
        tasks = [("power", p, w, 80, s, psi_case) for p, w in itertools.product(PSIS, WINDOWS) for s in range(a.seeds)]
    if a.smoke:
        tasks = tasks[:: max(1, len(tasks) // 12)][:12]
    t0 = time.time()
    rows = []
    with Pool(a.jobs) as pool:
        for i, r in enumerate(pool.imap_unordered(one, tasks, chunksize=1)):
            rows += r
            if i % 50 == 0:
                print(time.strftime("%F %T"), f"{i + 1}/{len(tasks)}", round(time.time() - t0), "s", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(out / f"{a.part}_runs.csv", index=False)
    key = ["psi_anchor", "case_window", "n", "prior", "comp"]
    if a.part == "bias":
        S = R.groupby(key).agg(d_star=("d_star", "mean"), d_star_sem=("d_star", "sem"), d=("d", "mean"), log_z2=("log_z2", "mean"),
                               log_V=("log_V", "mean"), rmse=("d_star", lambda x: float(np.sqrt(np.mean(x ** 2)))), runs=("d_star", "size")).reset_index()
    else:
        S = R.groupby(key).agg(reject=("p_star", lambda x: float(np.mean(x < 0.05))), d_star=("d_star", "mean"), log_V=("log_V", "mean"),
                               runs=("p_star", "size")).reset_index()
    S.to_csv(out / f"summary_{a.part}.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
