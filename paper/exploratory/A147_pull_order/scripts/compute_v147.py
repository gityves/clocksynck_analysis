"""A147: validation of the pull-versus-order readout on the A133 generative model (no GTEx donor-level data).

Simulator: A133 `simulate` (fast-death offsets phi, ventilated amplitudes and Poisson variances, 12 tissues, c = 0.45,
common displacement direction), 200 replicates per cell, n in {79, 298}. Reference group per replicate: independent
draw with f = 0, offsets intact, same n. Readouts: `pull_order` (draft clocksync function) on the raw cross-covariances,
all 66 pairs (all have >= 20 shared donors), 1,000 label permutations.
  pull arms: displacement fraction f in {0, 0.2, 0.44, 0.8}, offsets intact;
  loss arms: a share q in {0.25, 0.5, 1} of tissues (round(12 q), random per replicate) get offsets redrawn U[0, 2 pi),
             f in {0, 0.44}.
Verdicts V147.1 (AUC and f = 0.8 vs loss q = 0.5 f = 0 median) and V147.2 (ratio at f = 0, Spearman with f).

    python compute_v147.py [--subsample 1] [--seed 147] [--jobs 4]
Outputs: output/d147.csv (cells), output/reps147.csv (replicates), output/v147.csv (verdicts), output/failures.csv.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "A133_axial_decomposition" / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "A109_coupling_nulls" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr  # noqa: E402

import _helpers as H  # noqa: E402
import compute_v133_sims as S  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402
from pull_order import pull_order  # noqa: E402

FOLDER = "A147_pull_order"
SEED = 147
NS = (79, 298)
PULL_F = (0.0, 0.2, 0.44, 0.8)
LOSS_Q, LOSS_F = (0.25, 0.5, 1.0), (0.0, 0.44)


def cells():
    out = [dict(arm="pull", q=0.0, f=f, n=n) for n in NS for f in PULL_F]
    out += [dict(arm="loss", q=q, f=f, n=n) for n in NS for q in LOSS_Q for f in LOSS_F]
    return out


def cov(Z, M):
    return H.cross_moments(Z, M)[0]


def run_cell(k, cell, n_rep, n_perm, seed):
    a, V, phi = S.design()
    T = len(a)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 147, k]))
    rows, fails = [], []
    for rep in range(n_rep):
        try:
            Zr, Mr = S.simulate(cell["n"], 0.0, 0.0, "common", a, V, phi, rng)
            ph = phi.copy()
            m = int(round(cell["q"] * T))
            if m:
                idx = rng.choice(T, m, replace=False)
                ph[idx] = rng.uniform(0, 2 * np.pi, m)
            Zg, Mg = S.simulate(cell["n"], cell["f"], 0.0, "common", a, V, ph, rng)
            r = pull_order(cov(Zr, Mr), cov(Zg, Mg), n_perm=n_perm, rng=rng)
            if not np.isfinite(r["ratio"]):
                fails.append(dict(**cell, rep=rep, reason="ratio undefined (no same-pole pair with |ref dev| >= 1 h)"))
            rows.append(dict(**cell, rep=rep, **r))
        except Exception as e:  # noqa: BLE001
            fails.append(dict(**cell, rep=rep, reason=repr(e)))
    return rows, fails


def auc(x, y):
    """P(x > y) + 0.5 P(x = y)."""
    return float(mannwhitneyu(x, y, alternative="two-sided").statistic / (len(x) * len(y)))


def verdicts(reps):
    def g(arm, f, n, q=0.0):
        return reps[(reps.arm == arm) & (reps.f == f) & (reps.n == n) & (reps.q == q)]

    rows = []
    for n in NS:
        o_pull, o_loss = g("pull", 0.44, n).order, g("loss", 0.44 * 0 + 0.0, n, 0.5).order
        A = auc(o_pull.to_numpy(), o_loss.to_numpy())
        rows.append(dict(id="V147.1a", n=n, value=A, threshold=">= 0.9", passed=A >= 0.9, registered_n=(n == 79),
                         detail="AUC of order, pull f=0.44 vs loss q=0.5 f=0"))
        for fl in (0.0, 0.44):
            med = g("loss", fl, n, 0.5).order.median()
            share = float((g("pull", 0.8, n).order > med).mean())
            rows.append(dict(id="V147.1b" + ("" if fl == 0 else "_companion"), n=n, value=share, threshold=">= 0.80",
                             passed=share >= 0.8, registered_n=(n == 79 and fl == 0),
                             detail=f"share of pull f=0.8 replicates with order > median of loss q=0.5 f={fl} ({med:.3f})"))
        med = g("pull", 0.0, n).ratio.median()
        rows.append(dict(id="V147.2a", n=n, value=med, threshold="[0.85, 1.15]", passed=bool(0.85 <= med <= 1.15),
                         registered_n=True, detail="median shrink ratio at f=0"))
        cm = [g("pull", f, n).ratio.median() for f in PULL_F]
        rho = float(spearmanr(PULL_F, cm)[0])
        rows.append(dict(id="V147.2b", n=n, value=rho, threshold="<= -0.9", passed=rho <= -0.9, registered_n=True,
                         detail="Spearman of cell median ratio with f over pull arms; medians " + ", ".join(f"{c:.3f}" for c in cm)))
    v = pd.DataFrame(rows)
    reg = v[v.registered_n.astype(bool) | v.id.isin(["V147.2a", "V147.2b"])]
    v1 = bool(v[(v.id.isin(["V147.1a", "V147.1b"])) & (v.n == 79)].passed.all())
    v2 = bool(v[v.id.isin(["V147.2a", "V147.2b"])].passed.all())
    v = pd.concat([v, pd.DataFrame([
        dict(id="V147.1", n=79, value=np.nan, threshold="1a and 1b at n 79", passed=v1, registered_n=True, detail=""),
        dict(id="V147.2", n="79 and 298", value=np.nan, threshold="2a and 2b at both n", passed=v2, registered_n=True, detail=""),
        dict(id="PASS_BOTH", n="", value=np.nan, threshold="V147.1 and V147.2", passed=v1 and v2, registered_n=True, detail="")])],
        ignore_index=True)
    return v


def describe(reps):
    keys = ["arm", "q", "f", "n"]
    rows = []
    for key, d in reps.groupby(keys):
        rows.append(dict(zip(keys, key), reps=len(d), order_mean=d.order.mean(), order_sd=d.order.std(),
                         order_q05=d.order.quantile(.05), order_median=d.order.median(), order_q95=d.order.quantile(.95),
                         ratio_median=d.ratio.median(), ratio_q05=d.ratio.quantile(.05), ratio_q95=d.ratio.quantile(.95),
                         same_pole=d.same_pole.mean(), n_ratio_median=d.n_ratio.median(),
                         order_p_lt_005=float((d.p_order < 0.05).mean()), order_p_median=d.p_order.median()))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.set_defaults(seed=SEED)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = exploratory_dir(FOLDER, run["quick"])
    n_rep, n_perm = (10, 100) if run["quick"] else (200, 1000)
    t0 = time.time()
    cl = cells()
    res = Parallel(n_jobs=min(args.jobs, 4), verbose=5)(delayed(run_cell)(k, c, n_rep, n_perm, args.seed) for k, c in enumerate(cl))
    reps = pd.DataFrame([r for rows, _ in res for r in rows])
    fails = pd.DataFrame([f for _, fl in res for f in fl], columns=["arm", "q", "f", "n", "rep", "reason"])
    reps.to_csv(out / "reps147.csv", index=False)
    fails.to_csv(out / "failures.csv", index=False)
    d = describe(reps)
    d.to_csv(out / "d147.csv", index=False)
    v = verdicts(reps)
    v.to_csv(out / "v147.csv", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 90)
    print(d.round(3).to_string(index=False))
    print(v.round(4).to_string(index=False))
    print(f"failures: {len(fails)}; runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
