"""A374 (registered): how far can the organs have moved, given each group's observed arrangement correlation?

Truth = the full fast-death map of pair phases (A109 tissues, 66 pairs). A group with organ scatter x (SD, h) has pair
phases theta_f + (e_a - e_b) + eps_g, e ~ N(0, x^2) per organ; the fast map is observed as theta_f + eps_f; eps from
donor bootstraps (joint over pairs). 90 % compatibility set of x for the observed Fisher-Lee rho per group.

    python a374_interval.py [--check]
Outputs (`../output/`): intervals.csv, sim.csv. Aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A349_icu_limits" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a349_common as C  # noqa: E402
from a349_common import B, H, S  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
GRID = np.round(np.arange(0.0, 6.001, 0.1), 2)
N_BOOT, N_DRAW = 500, 500
PAIRS = np.array(H.PAIRS)


def wrap(a):
    return np.angle(np.exp(1j * a))


def phases(Z, M, idx):
    return H.pair_phase(N._H_of(Z, M, idx))


def boot_eps(Z, M, idx, rng, mask):
    base = phases(Z, M, idx)[mask]
    return np.array([wrap(phases(Z, M, rng.choice(idx, len(idx)))[mask] - base) for _ in range(N_BOOT)])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    global N_BOOT, N_DRAW
    if a.check:
        N_BOOT, N_DRAW = 10, 10
    OUT.mkdir(parents=True, exist_ok=True)
    T, fast, pm = C.tissues()
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    h0 = T.group_donors("vent")
    groups = {"slow deaths (Hardy 4)": T.group_donors("H4"), "ventilated, living brain": pm,
              "brain-dead organ donors": h0[coh[h0] == B.OPO], "Hardy 3 (descriptive)": T.group_donors("H3")}
    print(f"states ({time.time() - t0:.0f}s): fast {len(fast)}, " + ", ".join(f"{k} {len(v)}" for k, v in groups.items()), flush=True)
    rows, srows = [], []
    rng = np.random.default_rng([374, 1])
    for gs, kind in (("14 genes", "clock"), ("core", "core")):
        Z, M = S.raw(T, kind, 0, T.R)
        for g, idx in groups.items():
            mask = B.pair_mask(M, fast, idx)
            tf, tg = phases(Z, M, fast)[mask], phases(Z, M, idx)[mask]
            rho_obs = float(ST.fisher_lee(tf, tg))
            ef, eg = boot_eps(Z, M, fast, rng, mask), boot_eps(Z, M, idx, rng, mask)
            pa, pb = PAIRS[mask, 0], PAIRS[mask, 1]
            ok = []
            for x in GRID:
                sims = []
                for _ in range(N_DRAW):
                    e = rng.normal(0, x * W, len(H.BODY))
                    f_obs = tf + ef[rng.integers(N_BOOT)]
                    g_obs = tf + (e[pa] - e[pb]) + eg[rng.integers(N_BOOT)]
                    sims.append(float(ST.fisher_lee(f_obs, g_obs)))
                q05, q50, q95 = np.quantile(sims, [0.05, 0.5, 0.95])
                srows.append(dict(gene_set=gs, group=g, x_h=x, rho_q05=q05, rho_median=q50, rho_q95=q95))
                if q05 <= rho_obs <= q95:
                    ok.append(x)
            rows.append(dict(gene_set=gs, group=g, donors=len(idx), pairs=int(mask.sum()), rho_obs=rho_obs,
                             lower_h=float(min(ok)) if ok else np.nan, upper_h=float(max(ok)) if ok else np.nan,
                             upper_at_grid_edge=bool(ok and max(ok) >= GRID[-1]),
                             median_abs_pair_offset_fast_h=float(np.median(np.abs(tf)) / W),
                             median_folded_pair_offset_fast_h=float(np.median(np.minimum(np.abs(tf), np.pi - np.abs(tf))) / W)))
            print(f"{gs} / {g}: rho {rho_obs:.3f}, interval {rows[-1]['lower_h']} - {rows[-1]['upper_h']} h ({time.time() - t0:.0f}s)", flush=True)
    if not a.check:
        pd.DataFrame(rows).to_csv(OUT / "intervals.csv", index=False)
        pd.DataFrame(srows).to_csv(OUT / "sim.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
