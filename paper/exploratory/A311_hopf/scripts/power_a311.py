"""A311 power stage (registered): detection rate of the bootstrap part of H311.1 at the real GTEx design.

A104 rows and fit (fast deaths, 27 peripheral tissues, every tissue an anchor, leave-one-out host). Per tissue,
`clocksync.orientation.simulator`: 200 null simulations (isotropic own variation) give the null mean; each power
replicate draws one data set per tissue with a true tilt of 0.3 x E_o (E_o: own variance over the matched floor,
A294b) or 0 (type I), forms share_o = (T - null mean) / 2 psi_o, the median over tissues and its tissue-bootstrap
(2000) 90 % lower bound. Only simulated statistics are computed and written; the observed tilt is not.

    python power_a311.py [--reps 200]
Output: ../output/power.csv (one row per arm), ../output/power_tissues.csv (per tissue: 2 psi, E, null SD).
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync import orientation  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
FLOOR = HERE.parents[2] / "A294_floor" / "output" / "tissues.csv"
TILT = 0.3


def lower_bound(shares, rng, n_boot=2000):
    idx = rng.integers(0, len(shares), (n_boot, len(shares)))
    return float(np.quantile(np.median(shares[idx], axis=1), 0.05))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=200)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    limit_threads(4)
    t0 = time.time()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, name="gtex_periphery")
    floor = pd.read_csv(FLOOR, index_col=0)
    sims, rows = {}, []
    for k, t in enumerate(comps):
        obs, sim = orientation.simulator(ds, t)
        own = 2 * obs["psi"]
        E = float(floor.loc[t, "own"] - floor.loc[t, "own_matched_median"])
        rng = np.random.default_rng([311, k])
        null = np.array([sim(0.0, rng) for _ in range(200)])
        sims[t] = (sim, own, E, null.mean())
        rows.append(dict(tissue=t, own=own, E=E, tilt_true=TILT * E, null_sd=null.std(ddof=1),
                         share_sd=null.std(ddof=1) / own))
        del obs
    pt = pd.DataFrame(rows)
    print(f"{len(comps)} tissues, null simulations done ({time.time() - t0:.0f}s)", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pt.to_csv(OUT / "power_tissues.csv", index=False)
    rng = np.random.default_rng([311, 999])
    res = []
    for arm, scale in (("null", 0.0), ("tilt_0.3E", TILT)):
        med, lb = [], []
        for r in range(a.reps):
            sh = np.array([(sim(scale * E, rng) - m0) / own for sim, own, E, m0 in sims.values()])
            med.append(np.median(sh))
            lb.append(lower_bound(sh, rng))
        lb = np.array(lb)
        res.append(dict(arm=arm, reps=a.reps, median_share_mean=float(np.mean(med)), rate_lower_bound_gt_0=float(np.mean(lb > 0)),
                        expected_share_median=float(np.median(scale * pt.E / pt.own))))
        print(res[-1], f"({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(res).to_csv(OUT / "power.csv", index=False)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
