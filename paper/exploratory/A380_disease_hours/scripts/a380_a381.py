"""A380 (hour bounds of organ movement per chronic condition, relative to matched controls) and A381 (pooled metabolic
arrangement test), registered 2026-11-25.

Cases against all their matched controls (A379 matching; A381: T2D or BMI >= 30 against no T2D and BMI 18.5-25, up to
3 controls), the 12 A109 body tissues, pairs with >= 20 donors in both arms. A374 model: cases = control map + organ
scatter N(0, x^2) per organ + case bootstrap noise; control map observed with its own bootstrap noise.

Stages:
  power  A380: per condition, the upper bound of the 90 % compatibility set when rho_obs equals the median of
         rho_sim(0) (the detection limit in hours); A381: the same for the pooled metabolic arms, plus the SD and 5th
         percentile of r under 2000 within-set relabellings. No observed case-against-control value.
  score  (written after the power note)

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A380_disease_hours/scripts/a380_a381.py --stage power
Outputs (`../output/`): power.csv. Controlled access: labels and matched sets stay in memory; aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A379_disease_arrangement" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a379 as D  # noqa: E402
from a379 import B, H, ML, raw  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
GRID = np.round(np.arange(0.0, 6.001, 0.1), 2)
N_BOOT, N_DRAW, N_PERM, SEED = 500, 500, 2000, 380
PAIRS = np.array(H.PAIRS)


def metabolic_sets(q, T):
    pos = {d: i for i, d in enumerate(T.donors)}
    sets = []
    for coh in ("vent", "fast"):
        qc = q[q.cohort == coh]
        cas = (qc.t2d == 1) | (qc.bmi >= 30)
        con = (qc.t2d == 0) & (qc.bmi >= 18.5) & (qc.bmi < 25)
        sets += [[pos[d] for d in s] for s in ML.match_ratio(qc, cas, con, ratio=3)]
    return sets


def arms(sets):
    return np.array([s[0] for s in sets]), np.array([j for s in sets for j in s[1:]])


def phases(Z, M, idx, mask):
    return B.phase(Z, M, idx, mask)


def boot_eps(Z, M, idx, mask, rng):
    base = phases(Z, M, idx, mask)
    return np.array([np.angle(np.exp(1j * (phases(Z, M, rng.choice(idx, len(idx)), mask) - base)))
                     for _ in range(N_BOOT)])


def simulate(tc, ec, ek, mask, rng):
    """rho_sim quantiles per grid value x (h)."""
    pa, pb = PAIRS[mask, 0], PAIRS[mask, 1]
    out = []
    for x in GRID:
        sims = []
        for _ in range(N_DRAW):
            e = rng.normal(0, x * W, len(H.BODY))
            sims.append(float(ST.fisher_lee(tc + ec[rng.integers(N_BOOT)], tc + (e[pa] - e[pb]) + ek[rng.integers(N_BOOT)])))
        out.append(np.quantile(sims, [0.05, 0.5, 0.95]))
    return np.array(out)


def compat(q, rho_obs):
    ok = GRID[(q[:, 0] <= rho_obs) & (rho_obs <= q[:, 2])]
    return (float(ok.min()), float(ok.max())) if len(ok) else (np.nan, np.nan)


def power():
    limit_threads(4)
    t0 = time.time()
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    Z, M = raw(T, "clock", 0, T.R)
    q = D.covariates(T)
    print(f"states and covariates ({time.time() - t0:.0f}s)", flush=True)
    designs = {c: D.matched_sets(q, T, c) for c in D.CONDITIONS}
    designs["metabolic (A381)"] = metabolic_sets(q, T)
    rows = []
    for k, (name, sets) in enumerate(designs.items()):
        rng = np.random.default_rng([SEED, k])
        cases, ctl = arms(sets)
        mask = B.pair_mask(M, cases, ctl)
        tc = phases(Z, M, ctl, mask)
        sq = simulate(tc, boot_eps(Z, M, ctl, mask, rng), boot_eps(Z, M, cases, mask, rng), mask, rng)
        lo, hi = compat(sq, sq[0, 1])
        row = dict(design=name, sets=len(sets), cases=len(cases), controls=len(ctl), pairs=int(mask.sum()),
                   rho_sim0_median=float(sq[0, 1]), detection_upper_h=hi, upper_at_grid_edge=bool(hi >= GRID[-1]))
        if name.startswith("metabolic"):
            null = []
            for _ in range(N_PERM):
                a, b = [], []
                for s in sets:
                    p = rng.permutation(s)
                    a.append(p[0])
                    b += list(p[1:])
                null.append(float(ST.fisher_lee(phases(Z, M, np.array(b), mask), phases(Z, M, np.array(a), mask))))
            row.update(null_sd=float(np.std(null, ddof=1)), null_q05=float(np.quantile(null, 0.05)),
                       null_median=float(np.median(null)))
        rows.append(row)
        print(name, f"done ({time.time() - t0:.0f}s)", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power"], required=True)
    ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    power()


if __name__ == "__main__":
    main()
