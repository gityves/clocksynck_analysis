"""A205: gate and offset when the genes of one compartment shift by different amounts. Fits, then the registered summary.

Registered in paper/docs/registered_experiments.md (A205). One pickle per (scenario, window, rep) under <out>/runs/.

    RUN_CAP_GB=10 paper/common/run.sh paper/exploratory/A205_partial_shift_gate/scripts/compute_a205.py --jobs 8

Outputs: per_rep.csv, cells.csv, tests.csv.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import multiprocessing as mp
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent))
sys.path.insert(0, str(PAPER))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402

OUT = HERE.parent / "output"
HOURS = 24.0 / (2 * np.pi)
WF = cs.gtex()
GENES = list(WF.genes)
SEED, N_DONORS, C = 205, 80, 0.9
ANCHORS = ["A1", "A2", "A3", "A4"]
LOOP = ("ARNTL", "CLOCK", "NPAS2", "CRY1")
SCENARIOS = {
    "P0": None,
    "P12": {g: 12.0 for g in GENES},
    "P_partial": {g: 6.0 if g in LOOP else 12.0 for g in GENES},
    "P_partial3": {g: 3.0 if g in LOOP else 12.0 for g in GENES},
}
WINDOWS = {"F": ((0.0, 24.0), 0.0), "W": ((8.0, 18.0), 1.0)}


def host_variance(window, jitter) -> float:
    r = np.random.default_rng(0)
    t = r.uniform(*window, 400_000) + (r.normal(0, jitter, 400_000) if jitter > 0 else 0.0)
    return float(np.trace(np.cov(np.column_stack([np.cos(t / HOURS), np.sin(t / HOURS)]).T)))


def task(scen: str, window: str, rep: int, runs: str) -> dict:
    from clocksync.config import limit_threads

    limit_threads(1)
    f = Path(runs) / f"{scen}_{window}_{rep:03d}.pkl"
    if f.exists():
        return pd.read_pickle(f)
    win, jit = WINDOWS[window]
    V = host_variance(win, jit)
    psi = V * (1 - C) / (2 * C)  # c = 0.9 at this window (A145 psi_of, z = 1)
    comps = {k: Compartment(z=1.0, psi=psi) for k in ANCHORS}
    comps["T"] = Compartment(z=1.0, psi=psi, gene_shift_h=SCENARIOS[scen])
    tab, _ = cs.simulate(WF, comps, n=N_DONORS, window_h=win, jitter_h=jit,
                         seed=np.random.SeedSequence([SEED, list(SCENARIOS).index(scen), list(WINDOWS).index(window), rep]))
    tab = tab.rename(columns={"resection": "donor"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = cs.Dataset(tab, anchors=ANCHORS, name="a205").fit(n_boot=0, seed=rep, force=True)
    t, g = fit.table, fit.gate
    ref = S.circmean(t.loc[ANCHORS, "offset_h"].to_numpy())
    out = dict(scenario=scen, window=window, rep=rep, gate_pass=bool(g.passes), clock_share=g.clock_share,
               p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
               offset_h=float(S.wrap(t.loc["T", "offset_h"] - ref)) % 24,
               coupling=float(t.loc["T", "coupling"]),
               log_amp_rel=float(t.loc["T", "log_amp_rel"] - t.loc[ANCHORS, "log_amp_rel"].mean()))
    pd.to_pickle(out, f)
    return out


def summarise(df: pd.DataFrame, out: Path):
    from scipy import stats

    df.to_csv(out / "per_rep.csv", index=False)
    cells = df.groupby(["scenario", "window"]).agg(
        gate_pass=("gate_pass", "mean"), clock_share=("clock_share", "median"),
        offset_h=("offset_h", lambda x: S.circmean(x) % 24), offset_in_6_12=("offset_h", lambda x: float(((x >= 6) & (x <= 12)).mean())),
        coupling=("coupling", "median"), log_amp_rel=("log_amp_rel", "median"), reps=("rep", "size")).reset_index()
    cells.to_csv(out / "cells.csv", index=False)
    w = df.query("scenario == 'P_partial' and window == 'W'")
    n = len(w)
    k1, k2 = int(w.gate_pass.sum()), int(((w.offset_h >= 6) & (w.offset_h <= 12)).sum())
    tests = pd.DataFrame([
        dict(v="V205.1", readout="W P_partial gate pass share", k=k1, n=n, share=k1 / n, rule=">= 0.80",
             binom_lo95=stats.binomtest(k1, n).proportion_ci().low, verdict="PASS" if k1 / n >= 0.8 else "FAIL"),
        dict(v="V205.2", readout="W P_partial offset in [6, 12] h", k=k2, n=n, share=k2 / n, rule=">= 0.90",
             binom_lo95=stats.binomtest(k2, n).proportion_ci().low, verdict="PASS" if k2 / n >= 0.9 else "FAIL")])
    tests.to_csv(out / "tests.csv", index=False)
    print(tests.to_string(), flush=True)
    print(cells.round(3).to_string(), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    out = Path(args.out)
    runs = out / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    todo = [(s, w, r) for r in range(args.reps) for s in SCENARIOS for w in WINDOWS]
    print(f"workers = {args.jobs}; {len(todo)} tasks", flush=True)
    t0, rows, last = time.time(), [], time.time()
    with cf.ProcessPoolExecutor(args.jobs, mp_context=mp.get_context("spawn")) as ex:
        for i, r in enumerate(cf.as_completed([ex.submit(task, s, w, r, str(runs)) for s, w, r in todo])):
            rows.append(r.result())
            if time.time() - last > 300 or i + 1 == len(todo):
                print(f"[{time.time() - t0:.0f} s] {i + 1}/{len(todo)}", flush=True)
                last = time.time()
    summarise(pd.DataFrame(rows).sort_values(["scenario", "window", "rep"]), out)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
