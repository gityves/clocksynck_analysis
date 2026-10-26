"""PI Q10 (2026-10-26, post hoc diagnostic, no hypothesis): why does the A19 age model fail its count-level null?

At the six acceptance compartments (`calibration/compute_age_acceptance.py` picks), on the real count tables with ages
permuted within stratum (no age effect by construction), compares two null distributions of each A19 statistic:
  perm  the statistic over `--n-perm` within-stratum age permutations of the real data (the true null at this design);
  sim   the package's simulated null (`clocksync._core.trend.simulated_null`, `--n-sim` data sets), drawn at one permuted
        age vector, the reference `Dataset.covariate` uses for its `_p`.
Readouts per compartment and estimand: SD and IQR of each, their ratio perm / sim (> 1: the simulated null is too narrow,
p anti-conservative), the share of perm statistics outside the central 95 % of sim (0.05 if calibrated); and, per
permutation, the fit's coupling, log own noise and whether it sits at the psi floor (coupling-1 fits). Real ages are not
read (re-reading A19 on them waits for the PI's answer to Q10 (b)).

    python check_null.py [--n-perm 100] [--n-sim 200] [--seed 0]
Outputs: paper/exploratory/A19_acceptance_check/output/{null_compare,perm_fits}.csv
"""

from __future__ import annotations

import argparse
import sys
import zlib
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/extended_data/scripts/calibration"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core import covariate as C  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync._core import trend as T  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

import compute_age_acceptance as AA  # noqa: E402

STATS = [T._STAT[e] for e in T.ESTIMANDS]


def fit_stats(d, k, h, x):
    r = C.on_host(d, k, h, x)
    return {s: r.get(s, np.nan) for s in STATS} | {e: r.get(e, np.nan) for e in T.ESTIMANDS} | \
        dict(coupling=r.get("coupling", np.nan), x0=r.get("x0", np.nan), n=r.get("n"))


def compartment(name, k, ds, age, n_perm, n_sim, seed):
    rng = np.random.default_rng([seed, zlib.crc32(k.encode())])
    d = ds.data
    h = A._hosts(d, ds.anchors, [k])[k]
    perm = []
    for b in range(n_perm):
        a = AA.permuted_ages(ds, age, rng) / 10.0
        x = a - a.median()
        perm.append(fit_stats(d, k, h, x) | dict(perm=b))
    perm = pd.DataFrame(perm)
    a = AA.permuted_ages(ds, age, rng) / 10.0
    x = a - a.median()
    obs = C.on_host(d, k, h, x)
    sim = T.simulated_null(d, k, h, x, n_sim, rng, obs["x0"])[:, :len(STATS)]
    rows = []
    for j, s in enumerate(STATS):
        p, q = perm[s].dropna().to_numpy(), sim[:, j][np.isfinite(sim[:, j])]
        lo, hi = np.quantile(q, [0.025, 0.975])
        iqr = lambda v: np.subtract(*np.quantile(v, [0.75, 0.25]))  # noqa: E731
        rows.append(dict(design=name, compartment=k, stat=s, n_perm=len(p), n_sim=len(q), sd_perm=p.std(), sd_sim=q.std(),
                         sd_ratio=p.std() / q.std(), iqr_ratio=iqr(p) / iqr(q), med_perm=np.median(p), med_sim=np.median(q),
                         out95=float(np.mean((p < lo) | (p > hi)))))
    perm.insert(0, "compartment", k)
    perm.insert(0, "design", name)
    return pd.DataFrame(rows), perm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-perm", type=int, default=100)
    ap.add_argument("--n-sim", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--only", default=None, help="one compartment, e.g. lung:AT1 (quick check)")
    args = ap.parse_args()
    limit_threads(1)
    warnings.simplefilter("ignore")
    out = HERE / "output"
    out.mkdir(exist_ok=True)
    comp, fits = [], []
    for name in AA.DESIGNS:
        pb, donors, age, anc = AA.design(name, None, args.seed)
        ds = cs.Dataset(pb, anchors=anc, donors=donors, name=name)
        ks, c = AA.picks(ds, age, AA.CA.MIN_AGE_DONORS)
        for k in ks:
            if args.only and k != args.only:
                continue
            r, f = compartment(name, k, ds, age, args.n_perm, args.n_sim, args.seed)
            print(r.round(3).to_string(index=False), flush=True)
            comp.append(r)
            fits.append(f)
    pd.concat(comp).to_csv(out / "null_compare.csv", index=False)
    pd.concat(fits).to_csv(out / "perm_fits.csv", index=False)


if __name__ == "__main__":
    main()
