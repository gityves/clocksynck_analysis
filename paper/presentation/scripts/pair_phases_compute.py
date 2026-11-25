"""Slide data (results, drives, 4.3): organ-pair phases in fast deaths, ventilated donors with a living brain and brain
death, the quantity behind the arrangement rho of A293 / A317.

A317 inputs as registered (`score_a317.py`): A109 `Tissues(200, seed 0)` states (one clock-plane projection per tissue,
no host), 12 body tissues, the 65 pairs with enough donors in every group; GTEx fast deaths (Hardy 1-2), ventilated
post-mortem donors (79) and brain-dead organ donors drawn to 79 (20 draws, seed 317). Per pair: phase of
E[z_a conj z_b] in hours. Fast and ventilated: all donors of the group; brain death: circular mean over the 20 draws.
Intervals: 95 % donor bootstrap (200 draws) of a group-sized estimate (brain death: 79 donors resampled from all organ
donors). Descriptive display of the registered statistic. Aggregates only (one row per tissue pair).

    python scripts/pair_phases_compute.py
Output: presentation/output/pair_phases/pairs.csv
"""
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER / "exploratory" / "A317_braindeath_confirm" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a317 as A  # noqa: E402
from score_a317 import B, H, gtex, raw  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = PAPER / "presentation" / "output" / "pair_phases"
N_BOOT, SEED = 200, 4202
HOURS = 12 / np.pi


def circ_mean(x, axis=0):
    return np.angle(np.exp(1j * x / HOURS).mean(axis)) * HOURS


def interval(boot, point):
    """95 % interval of bootstrap estimates around the point estimate (deviations wrapped to +-12 h)."""
    dev = np.angle(np.exp(1j * (boot - point) / HOURS)) * HOURS
    return point + np.quantile(dev, 0.025, axis=0), point + np.quantile(dev, 0.975, axis=0)


def main():
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(A.SEED)  # the A317 draws
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    Z, M = raw(T, "clock", 0, T.R)
    mask = B.pair_mask(M, fast, pm, *od_draws[:1]) & B.pair_mask(M, od)
    ph = lambda g: -B.phase(Z, M, g, mask) * HOURS  # noqa: E731  (clocksync sign negated: + = a peaks after b)
    pairs = [(H.BODY[a], H.BODY[b]) for (a, b), m in zip(H.PAIRS, mask) if m]
    tab = pd.DataFrame(pairs, columns=["tissue_a", "tissue_b"])
    rng = np.random.default_rng(SEED)
    for g, donors, pool in (("fast", fast, fast), ("ventilated", pm, pm), ("brain_death", None, od)):
        point = circ_mean(np.array([ph(d) for d in od_draws])) if donors is None else ph(donors)
        n = len(pm) if donors is None else len(donors)
        boot = np.array([ph(rng.choice(pool, n)) for _ in range(N_BOOT)])
        lo, hi = interval(boot, point)
        tab[g], tab[f"{g}_lo"], tab[f"{g}_hi"], tab[f"{g}_n"] = point, lo, hi, n
    tab["brain_death_all"] = ph(od)
    tab.to_csv(OUT / "pairs.csv", index=False)
    print(f"{len(tab)} pairs; fast {len(fast)}, ventilated {len(pm)}, brain death {len(od)} (drawn to {len(pm)})")
    print(tab.drop(columns=["tissue_a", "tissue_b"]).describe().round(2).to_string())


if __name__ == "__main__":
    main()
