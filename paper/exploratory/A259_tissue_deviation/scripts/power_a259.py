"""A259 power note (simulation only): bias, interval coverage and detection limit of `clocksync.deviation.incoherence`
at the designs of A259, before any real data are read.

Simulated donors (`clocksync.simulate`, deaths uniform over 24 h): ten host compartments (relative amplitudes 0.8-1.25,
own noise psi 0.01) and one target with relative amplitude 0.7 and a person-specific deviation of SD sigma hours
(isotropic in the clock plane: psi = kappa |z|^2 A^2 / 2 with kappa = (sigma pi / 12)^2). Cells: clock amplitude
A 0.5 / 1.0, donors 100 / 200 / 400, sigma 0 / 0.25 / 0.5 / 0.75 / 1.0 / 1.5 h; a shared non-clock program across a
donor's compartments (shared_sd 0.1) at A 1.0, 200 donors; the pair statistic (one host compartment, relative amplitude
0.95, psi 0.01) at 680 donors, A 0.5 / 1.0. 40 data sets per cell, 10 gene splits, 10 host splits, 200 bootstraps.

    python power_a259.py [--reps 40] [--jobs 4]
Outputs (output/power/): runs.csv, summary.csv (bias in kappa and hours, rate of lower 90 % bound > 0, rate of upper
90 % bound below the truth).
"""

from __future__ import annotations

import argparse
import itertools
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync import deviation as D
from clocksync.config import limit_threads

OUT = Path(__file__).resolve().parents[1] / "output" / "power"
SEED = 259
WF = cs.gtex()
SIGMAS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5]
ZT = 0.7 * np.exp(-1j)


def one(task):
    mode, amp, n, sigma, shared, rep = task
    limit_threads(1)
    warnings.simplefilter("ignore")
    kappa = (sigma * np.pi / 12) ** 2
    psi_t = kappa * abs(ZT) ** 2 * amp ** 2 / 2
    if mode == "target":
        comps = {f"h{i}": cs.Compartment((0.8 + 0.05 * i) * np.exp(1j * 0.3 * i), psi=0.01) for i in range(10)}
        hosts = list(comps)
    else:
        comps = {"h": cs.Compartment(0.95 * np.exp(0.5j), psi=0.01)}
        hosts = ["h"]
    comps["T"] = cs.Compartment(ZT, psi=psi_t)
    df, _ = cs.simulate(WF, comps, n=n, seed=[SEED, rep, n, int(amp * 10), int(sigma * 100), int(shared * 10)],
                        window_h=(0, 24), jitter_h=0.0, amplitude=amp, shared_sd=shared)
    r = D.incoherence(df.rename(columns={"resection": "donor"}), "T", hosts, n_splits=10, n_host_splits=10,
                      n_boot=200, seed=rep)
    k_h = 2 * 0.01 / (0.95 * amp) ** 2  # the host compartment's own kappa (pair mode)
    truth = kappa if mode == "target" else (1 + kappa) * (1 + k_h) - 1
    return dict(mode=mode, amp=amp, n=n, sigma=sigma, shared=shared, rep=rep, kappa_true=truth, kappa=r.kappa,
                lo=r.lo, hi=r.hi, hours=r.hours, hours_hi=r.hours_hi, hours_true=float(D.hours(truth)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=40)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [("target", A, n, s, 0.0, r) for A, n, s in itertools.product((0.5, 1.0), (100, 200, 400), SIGMAS)
             for r in range(a.reps)]
    tasks += [("target", 1.0, 200, s, 0.1, r) for s in SIGMAS for r in range(a.reps)]
    tasks += [("pair", A, 680, s, 0.0, r) for A in (0.5, 1.0) for s in SIGMAS for r in range(a.reps)]
    t0, rows = time.time(), []
    with Pool(a.jobs) as pool:
        for i, r in enumerate(pool.imap_unordered(one, tasks, chunksize=2)):
            rows.append(r)
            if i % 100 == 0:
                print(time.strftime("%F %T"), f"{i + 1}/{len(tasks)}", round(time.time() - t0), "s", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "runs.csv", index=False)
    S = R.groupby(["mode", "amp", "n", "shared", "sigma"]).apply(lambda d: pd.Series(dict(
        kappa_true=d.kappa_true.iloc[0], kappa_mean=d.kappa.mean(), kappa_bias=(d.kappa - d.kappa_true).mean(),
        hours_true=d.hours_true.iloc[0], hours_median=d.hours.median(), hours_hi_median=d.hours_hi.median(),
        detect=(d.lo > 0).mean(), hi_below_truth=(d.hi < d.kappa_true).mean(), lo_above_truth=(d.lo > d.kappa_true).mean(),
        runs=len(d))), include_groups=False).reset_index()
    S.to_csv(OUT / "summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
