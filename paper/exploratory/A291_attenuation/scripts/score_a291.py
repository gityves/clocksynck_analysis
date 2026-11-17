"""A291 (registered in `paper/docs/registered_experiments.md`): is "stronger clocks are more coupled" (A104) what
measurement alone produces?

Rows and fit as A104 (`results/figure_2/scripts/backbone/compute_oscillator_gtex.py`: GTEx fast deaths, tissues with
>= 40 donors, periphery and brain blocks, every tissue an anchor, leave-one-out host). Per tissue,
`clocksync.simulate.like` draws data at its real design (host means and posterior covariances, measurement
covariances, strata, intercepts): 100 replicates at true coupling 1 (perfect locking) and 100 at the estimate.
c_ceil = median fitted coupling at truth 1. H291.1: Spearman rho(A, c / c_ceil) > 0 over the periphery tissues,
one-sided permutation p < 0.05. Aggregates only (GTEx controlled access).

    python score_a291.py [--n-rep 100] [--jobs 4] [--check]
Outputs (`../output/`): tissues.csv (one row per tissue), tests.csv.
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
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_PERM = 291, 10000
BLOCKS = ("periphery", "brain")


def datasets():
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    out = {}
    for block in BLOCKS:
        comps = sorted(c for c in pb.compartment.unique() if c.startswith("Brain") == (block == "brain"))
        out[block] = cs.Dataset(pb[pb.compartment.isin(comps)], anchors=comps, name=f"gtex_{block}")
    return out


def simulate_tissue(ds, block, k, idx, n_rep):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    like = cs.simulate.like(ds, k)
    row = dict(block=block, tissue=k)
    for truth, tag, j in ((1.0, "ceil", 1), (float(like.estimate["coupling"]), "at_est", 2)):
        r = like.replicates(truth, n=n_rep, seed=int(np.random.SeedSequence([SEED, BLOCKS.index(block), idx, j]).generate_state(1)[0]))
        row[f"c_{tag}"] = float(r.coupling.median())
        row[f"c_{tag}_q05"], row[f"c_{tag}_q95"] = (float(x) for x in r.coupling.quantile([0.05, 0.95]))
        row[f"boundary_{tag}"] = float(r.boundary.mean())
    row["truth_at_est"] = float(like.estimate["coupling"])
    return row


def perm_rho(x, y, seed):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    r0 = spearmanr(x, y).statistic
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(rng.permutation(x), y).statistic for _ in range(N_PERM)])
    return float(r0), float((1 + np.sum(null >= r0)) / (1 + N_PERM)), int(ok.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-rep", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    dss = datasets()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    obs = {b: A104.per_tissue(ds) for b, ds in dss.items()}
    rep = {b: A104.rho(t) for b, t in obs.items()}
    print(f"reproduction: periphery rho {rep['periphery']:.3f} (A104 0.734), brain {rep['brain']:.3f} (-0.115) "
          f"{time.time() - t0:.0f}s", flush=True)
    if round(rep["periphery"], 3) != 0.734 or round(rep["brain"], 3) != -0.115:
        raise SystemExit("A104 not reproduced: stop (registered check)")
    tasks = [(b, k, i) for b in BLOCKS for i, k in enumerate(obs[b].index)]
    if a.check:
        tasks, a.n_rep = tasks[:1], 5
    rows = Parallel(n_jobs=a.jobs)(delayed(simulate_tissue)(dss[b], b, k, i, a.n_rep) for b, k, i in tasks)
    sim = pd.DataFrame(rows)
    if a.check:
        print(sim.round(3).to_string(), flush=True)
        return
    t = pd.concat([o.assign(block=b) for b, o in obs.items()]).rename_axis("tissue").reset_index()
    t = t.merge(sim, on=["block", "tissue"])
    t["c_over_ceil"] = t.coupling / t.c_ceil
    t["design_bias"] = t.c_at_est - t.truth_at_est
    t.to_csv(OUT / "tissues.csv", index=False)
    tests = []
    for i, b in enumerate(BLOCKS):
        g = t[t.block == b]
        A = g.own_amplitude.to_numpy(float)
        r_corr, p_corr, n = perm_rho(A, g.c_over_ceil.to_numpy(float), SEED + i)
        r_ceil, p_ceil, _ = perm_rho(A, g.c_ceil.to_numpy(float), SEED + 10 + i)
        tests.append(dict(block=b, analysis="H291.1" if b == "periphery" else "descriptive", tissues=n,
                          rho_obs=rep[b], rho_corrected=r_corr, p_one_sided=p_corr, rho_ceiling=r_ceil,
                          p_ceiling=p_ceil, median_c_ceil=float(g.c_ceil.median()), min_c_ceil=float(g.c_ceil.min()),
                          median_design_bias=float(g.design_bias.median())))
    tests = pd.DataFrame(tests)
    tests["passes"] = (tests.analysis == "H291.1") & (tests.rho_corrected > 0) & (tests.p_one_sided < 0.05)
    tests.to_csv(OUT / "tests.csv", index=False)
    print(tests.round(3).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
