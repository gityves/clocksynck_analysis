"""Roadmap item A1 (b), simulation (not a registered hypothesis): are offset contrasts between two groups biased when
their host spreads differ, and does label permutation hold its size? A258 answered this for d*, log z2 and log V only.

Design and code of A258 (`A258_spread_bias/scripts/simulate_a258.py`): two arms of simulated donors (`clocksync.simulate`),
three anchors and one target (T cell, offset -3 h, psi 0.05); reference window 6-18 h, case window 6-18 / 8-16 / 10-14 /
11-13 h; anchor noise psi 0.002 / 0.05 / 0.15; donors per arm 40 / 80 / 200. Readout: `offset_diff_h` of
`Dataset.contrast` (arg z_case / z_control in h; case minus control) and its permutation p (`p_offset`), with host_prior
"pooled" (used by A208, A231, A233, A235) and "arm".
  bias   truth 0 h, 100 seeds per cell, no permutations
  size   truth 0 h, 80 donors per arm, 30 seeds per cell, 200 permutations
  shift  truth +0.5 h in the case arm's target, 80 donors per arm, 30 seeds per cell: recovery of a known shift

Reading fixed before the run: an equivalence within 1 h stands if, at the windows of the scored contrasts, |mean bias| of
the pooled-prior offset is below 0.1 h (a tenth of the bound) and the permutation rejection rate at 0.05 is at most 0.10
(A258's tolerance for size); a known +0.5 h shift is recovered within 0.1 h. Otherwise the cell is flagged.

    python simulate_offsets.py --part bias|size|shift [--seeds N] [--jobs 6] [--smoke]
Outputs (output/): offsets_<part>_runs.csv, offsets_<part>_summary.csv.
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
sys.path.insert(0, str(HERE.parents[1] / "A258_spread_bias" / "scripts"))
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from simulate_a258 import ANCH, NS, PSI_T, PSIS, WF, WINDOWS, Z  # noqa: E402

SEED = 1001


def cohort(n, seed, window, psi_a, prefix, shift_h=0.0):
    comps = {k: cs.Compartment(Z[k], psi=psi_a, mesor_sd=0.0) for k in ANCH}
    comps["T cell"] = cs.Compartment(Z["T cell"] * np.exp(1j * shift_h / HOURS), psi=PSI_T, mesor_sd=0.0)
    df, _ = cs.simulate(WF, comps, n=n, seed=seed, window_h=window, jitter_h=1.0, tau=0.0)
    df["resection"] = prefix + df.resection.str[1:]
    return df.rename(columns={"resection": "donor"})


def one(task):
    part, psi_a, wb, n, s, shift = task
    limit_threads(1)
    warnings.simplefilter("ignore")
    a = cohort(n, [SEED, 1, s, n, int(psi_a * 1000)], (6, 18), psi_a, "a")
    b = cohort(n, [SEED, 2, s, n, int(psi_a * 1000), wb[0]], wb, psi_a, "b", shift)
    g = pd.Series(["control"] * n + ["case"] * n, index=sorted(a.donor.unique()) + sorted(b.donor.unique()))
    ds = cs.Dataset(pd.concat([a, b], ignore_index=True), anchors=ANCH)
    rows = []
    for prior in ("pooled", "arm"):
        n_perm = 200 if part == "size" else 0
        t = ds.contrast(g, n_perm=n_perm, n_boot=0, blocks="none", host_prior=prior, seed=s).table
        for k in ("T cell", "fibroblast"):
            r = t.loc[k]
            rows.append(dict(part=part, psi_anchor=psi_a, case_window=f"{wb[0]}-{wb[1]}", n=n, seed=s, prior=prior,
                             comp=k, truth_h=shift if k == "T cell" else 0.0, offset_h=r.offset_diff_h,
                             p_offset=r.p_offset if n_perm else np.nan, log_V=r.log_V_ratio, d_star=r.d_star))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["bias", "size", "shift"], required=True)
    ap.add_argument("--seeds", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    seeds = a.seeds or (100 if a.part == "bias" else 30)
    out = OUT / ("smoke" if a.smoke else "")
    out.mkdir(parents=True, exist_ok=True)
    if a.part == "bias":
        tasks = [("bias", p, w, n, s, 0.0) for p, w, n in itertools.product(PSIS, WINDOWS, NS) for s in range(seeds)]
    else:
        shift = 0.5 if a.part == "shift" else 0.0
        tasks = [(a.part, p, w, 80, s, shift) for p, w in itertools.product(PSIS, WINDOWS) for s in range(seeds)]
    if a.smoke:
        tasks = tasks[:: max(1, len(tasks) // 12)][:12]
    t0, rows = time.time(), []
    with Pool(a.jobs) as pool:
        for i, r in enumerate(pool.imap_unordered(one, tasks, chunksize=1)):
            rows += r
            if i % 100 == 0:
                print(time.strftime("%F %T"), f"{i + 1}/{len(tasks)}", round(time.time() - t0), "s", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(out / f"offsets_{a.part}_runs.csv", index=False)
    R["err"] = R.offset_h - R.truth_h
    key = ["psi_anchor", "case_window", "n", "prior", "comp"]
    agg = dict(bias_h=("err", "mean"), bias_sem=("err", "sem"), median_err_h=("err", "median"),
               rmse_h=("err", lambda x: float(np.sqrt(np.nanmean(x ** 2)))), log_V=("log_V", "mean"),
               runs=("err", "size"))
    if a.part == "size":
        agg["reject"] = ("p_offset", lambda x: float(np.mean(x < 0.05)))
    S = R.groupby(key).agg(**agg).reset_index()
    S.to_csv(out / f"offsets_{a.part}_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
