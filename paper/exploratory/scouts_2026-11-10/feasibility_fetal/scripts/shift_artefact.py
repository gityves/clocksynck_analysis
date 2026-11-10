"""A shift common to the clock genes of a sample is read as clock-gene coherence by the unsigned statistic.

Simulation (`clocksync.simulate`, one compartment, no clock or an adult clock, per-sample shift of SD 0, 0.1 or 0.2
natural-log units on all clock genes): the unsigned coherence R of `decompose.components` (the L1 statistic up to
A189) and the signed coherence C of `decompose.signed_coherence` (offset fitted per sample, sign kept; A237 onwards),
each with its per-gene permutation p. Found while testing `common.development` (2026-11-10).

    python shift_artefact.py [--n-seed 40] [--n-perm 500] [--jobs 4]

Output (`../output/`): shift_artefact.csv (one row per scenario: mean statistic and share of p < 0.05).
"""
import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, logfrac, signed_coherence  # noqa: E402

OUT = HERE.parents[1] / "output"
WF = cs.gtex()
N_PERM = 500
STATS = {"R_unsigned": lambda y: components(y, np.zeros(y.shape[1]), 1, WF)["coherence"],
         "C_signed": lambda y: signed_coherence(y, WF)}


def one(args):
    amplitude, shift, n, seed = args
    tab, _ = cs.simulate(WF, {"liver": cs.Compartment(1.0, mesor_sd=0.0)}, n=n, seed=seed, window_h=(0, 24),
                         amplitude=amplitude, shift_sd=shift)
    y = logfrac(tab, list(WF.genes))[0]
    y = y - y.mean(0)
    rng, row = np.random.default_rng([1110, seed]), dict(amplitude=amplitude, shift_sd=shift, n=n, seed=seed)
    perms = [np.column_stack([rng.permutation(c) for c in y.T]) for _ in range(N_PERM)]
    for name, f in STATS.items():
        obs = f(y)
        row[name], row[f"p_{name}"] = obs, (1 + sum(f(p) >= obs for p in perms)) / (1 + N_PERM)
    return row


def main():
    global N_PERM
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-seed", type=int, default=40)
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    N_PERM = a.n_perm
    limit_threads(1)
    tasks = [(amp, s, n, seed) for amp in (0.0, 1.0) for s in (0.0, 0.1, 0.2) for n in (10, 40) for seed in range(a.n_seed)]
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        d = pd.DataFrame(pool.map(one, tasks))
    out = d.groupby(["amplitude", "shift_sd", "n"]).agg(
        n_sim=("seed", "size"), R_unsigned=("R_unsigned", "mean"), C_signed=("C_signed", "mean"),
        pass_R_unsigned=("p_R_unsigned", lambda p: (p < .05).mean()),
        pass_C_signed=("p_C_signed", lambda p: (p < .05).mean())).reset_index()
    out.to_csv(OUT / "shift_artefact.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
