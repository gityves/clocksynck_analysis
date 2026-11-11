"""Adult power calibration of A242 (replication of the fetal adrenal clock on open data), before registration.

Read from the fetal data: the sample sheet (cohort, block, age, depth) and the cohort medians of the clock-gene counts
(`../output/clock_detection.csv`, written by `prep_adrenal.py`). No fetal clock statistic is computed.

Every fetus of a test (`score_a242.tests`) is stood for by one distinct GTEx adrenal gland (death at 08-18 h, any
Hardy class), which receives its cohort, age and depth. The adult clock-gene counts are thinned binomially gene by
gene to what the cohort captures: target count = cohort median count of the gene x depth of the fetus / median depth
of the cohort; where the adult libraries hold fewer reads than that, the targets of the gene are lowered by one
factor over the cohort (`cap`; it cancels with the cohort indicator). The adjustment and the test are the registered
ones (`common.development`: cohort indicator and age trend; signed coherence C; per-gene permutations within cohort).
For the Schroder design the surrogate components are those of the adult libraries, thinned to the fetal depths.

Settings: `a1`, `a0.75`, `a0.5`: the adult clock at relative amplitude a (`clocksync.simulate.scale_clock` over the
pool); `shuffled`: every clock gene permuted across the drawn donors on its own (no clock; the false-positive rate).
`p_global` (settings a1 and shuffled): the same test with permutations across all fetuses, as in A237-A239.
Controlled access: GTEx time of death is used in memory only; outputs are aggregates, without donor identifiers.

    python feasibility_a242.py [--n-sub 200] [--jobs 8] [--check]

Outputs (`../output/`): feasibility_design.csv, feasibility_subsamples.csv, feasibility_passrates.csv (`_check`).
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import score_a242 as a242  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_PERM, WF = 242, 500, a242.WF
GENES = list(WF.genes)
SCALES = (1.0, 0.75, 0.5)
POWERED = ["late", "late:deep", "early", "all", "kameneva", "wang", "schroder", "schroder:pc0", "schroder:pc2"]
S = {}  # state inherited by the workers


def setup():
    """Fetal sheet and tests, cohort targets, adult pool (clock genes at every scale; genome for the surrogates)."""
    obs = a242.sheet()
    obs["pc1"], obs["pc2"] = 0.0, 0.0
    full = pd.read_csv(a242.SRC / "samples.csv")
    pb, d = gtex.load(["Adrenal Gland"], hardy=(0, 1, 2, 3, 4), genes="all", depth="library_size")
    pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))].reset_index(drop=True)
    depth, rng = pb.depth.to_numpy(float), np.random.default_rng([SEED, 0])
    clock = pb[["donor", "compartment", "depth", *GENES]]
    genome = pb.drop(columns=["donor", "compartment", "depth", *GENES])
    genome = genome.loc[:, genome.mean(axis=0) >= 50]
    S.update(obs=obs, tests={k: v for k, v in a242.tests(obs)[0].items() if k in POWERED}, depth=depth,
             det=pd.read_csv(OUT / "clock_detection.csv", index_col=0)[GENES],
             med_depth=full.groupby("cohort").depth.median(),
             x={a: scale_clock(clock, WF, a, depth=clock.depth, rng=rng)[GENES].to_numpy(np.int64) for a in SCALES},
             genome=genome.to_numpy(np.int64), genome_genes=list(genome.columns),
             fmed=np.median(clock[GENES].to_numpy(float) / depth[:, None], axis=0), dref=depth.min())
    return pd.DataFrame([dict(test=k, n_fetuses=int(r.sum()), n_blocks=obs[r].block.nunique(),
                              n_covariates=c(obs[r]).shape[1], pool_donors=len(pb),
                              median_cap=np.median(cap(obs[r])), share_capped=np.mean(cap(obs[r]) < 1))
                         for k, (r, c, _) in S["tests"].items()])


def target(f):
    """Clock-gene counts the cohort would capture in each fetus of `f` (fetuses x genes)."""
    return S["det"].loc[f.cohort].to_numpy(float) * (f.depth / f.cohort.map(S["med_depth"])).to_numpy(float)[:, None]


def cap(f):
    """Factor (<= 1) on the targets of each gene, one per cohort: the adult libraries cannot be thinned upwards."""
    top = pd.DataFrame(target(f)).groupby(f.cohort.to_numpy()).transform("max").to_numpy()
    return np.minimum(1.0, S["fmed"] * S["dref"] / top)


def stand_in(f, x, depth, rng):
    """Table of the adult counts `x` (library sizes `depth`) thinned to the design `f`."""
    p = cap(f) * target(f) / (S["fmed"] * depth[:, None])
    t = f.reset_index(drop=True).copy()
    t[GENES] = rng.binomial(x, np.minimum(p, 1.0)).astype(float)
    return t


def one(i):
    rng, out = np.random.default_rng([SEED, 1, i]), []
    for j, (name, (rows, cov, strata)) in enumerate(S["tests"].items()):
        f = S["obs"][rows]
        pick = rng.choice(len(S["depth"]), len(f), replace=False)
        depth = S["depth"][pick]
        if name.startswith("schroder"):  # surrogate components of the adult libraries at the fetal depths
            g = rng.binomial(S["genome"][pick], np.minimum(1.0, f.depth.to_numpy(float) / depth)[:, None])
            f = f.copy()
            f[["pc1", "pc2"]] = dev.surrogates(pd.DataFrame(g), f.depth, a242.blocks(f), 2)
        null = np.column_stack([np.round(rng.permutation(c / depth) * depth) for c in S["x"][1.0][pick].T])
        for setting, x in [*((f"a{a:g}", S["x"][a][pick]) for a in SCALES), ("shuffled", null.astype(np.int64))]:
            tab = dev.adjust(stand_in(f, x, depth, rng), GENES, cov)[0]
            c, p, _ = dev.coherence(tab, WF, N_PERM, np.random.default_rng([SEED, 2, i, j]), strata)
            pg = (dev.coherence(tab, WF, N_PERM, np.random.default_rng([SEED, 3, i, j]))[1]
                  if setting in ("a1", "shuffled") else np.nan)
            out.append(dict(test=name, setting=setting, sub=i, C=c, p=p, p_global=pg))
    return out


def ci(k, n):
    r = binomtest(int(k), int(n)).proportion_ci(0.95, method="exact")
    return r.low, r.high


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    pd.set_option("display.width", 250)
    design = setup()
    design.to_csv(OUT / f"feasibility_design{sfx}.csv", index=False)
    print(design.round(3).to_string(index=False), flush=True)
    n, res = 4 if a.check else a.n_sub, []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, range(n))):
            res.extend(r)
            if (k + 1) % 10 == 0 or k + 1 == n:
                print(f"{k + 1}/{n} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res).sort_values(["test", "setting", "sub"])
    res.to_csv(OUT / f"feasibility_subsamples{sfx}.csv", index=False)
    rows = []
    for (test, setting), r in res.groupby(["test", "setting"], sort=False):
        k, m = int((r.p < .05).sum()), len(r)
        lo, hi = ci(k, m)
        rows.append(dict(test=test, setting=setting, n_subsamples=m, pass_rate=k / m, ci_lo=lo, ci_hi=hi,
                         pass_rate_global=(r.p_global < .05).mean() if r.p_global.notna().any() else np.nan,
                         median_C=r.C.median(), C_q05=r.C.quantile(.05)))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"feasibility_passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
