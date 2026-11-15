"""Adult power of the A243 readouts (fetal heart, kidney, lung, liver pooled over open cohorts), before any fetal clock
statistic. As `A242/scripts/feasibility_a242.py`: every fetus of a readout is stood for by one GTEx donor of the adult organ
(heart = Heart - Left Ventricle, kidney = Kidney - Cortex, lung = Lung, liver = Liver; deaths at 08-18 h, any Hardy class;
distinct donors within a cohort), which receives its cohort, age and depth; the adult clock-gene counts are thinned gene by
gene to what the cohort captures (`../output/clock_detection.csv` median count x depth of the fetus / median depth of the
cohort; capped where the adult library holds fewer reads, one factor per gene and cohort). Settings a1, a0.75, a0.5
(`simulate.scale_clock`) and shuffled (clock genes permuted across the drawn donors: false-positive rate). 200 subsamples,
500 permutations within cohort. Readouts: `<organ>:all`, `:late`, `:early`, `:all:new` and cohorts with >= 10 fetuses alone.
Controlled access: GTEx time of death is used in memory only; outputs are aggregates.

    python feasibility_a243.py [--n-sub 200] [--jobs 8] [--check]
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
sys.path.insert(0, str(HERE.parents[2] / "A242_fetal_adrenal_replication" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a243 as d  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from feasibility_a242 import ci  # noqa: E402

OUT, WF, SEED, N_PERM = d.OUT, d.WF, d.SEED, 500
GENES = list(WF.genes)
SCALES = (1.0, 0.75, 0.5)
TISSUE = {"heart": "Heart - Left Ventricle", "kidney": "Kidney - Cortex", "lung": "Lung", "liver": "Liver"}
S = {}


def setup():
    obs = d.sheet()
    pb, don = gtex.load(list(TISSUE.values()), hardy=(0, 1, 2, 3, 4), genes="clock", depth="library_size")
    pb = pb[pb.donor.isin(set(don.index[in_window(don.tod.to_numpy())]))]
    rng = np.random.default_rng([SEED, 0])
    pools = {}
    for o, t in TISSUE.items():
        x = pb[pb.compartment == t].reset_index(drop=True)
        depth = x.depth.to_numpy(float)
        pools[o] = dict(depth=depth, fmed=np.median(x[GENES].to_numpy(float) / depth[:, None], axis=0), dref=depth.min(),
                        x={a: scale_clock(x[["donor", "compartment", "depth", *GENES]], WF, a, depth=x.depth, rng=rng)[GENES]
                           .to_numpy(np.int64) for a in SCALES})
    det = pd.read_csv(OUT / "clock_detection.csv")
    det["key"] = det.organ + "|" + det.cohort
    obs["key"] = obs.organ + "|" + obs.cohort
    tests = {k: v for k, v in d.all_tests(obs).items() if ":-" not in k}
    S.update(obs=obs, tests=tests, pools=pools, det=det.set_index("key")[GENES], medd=obs.groupby("key").depth.median())
    rows = []
    for k, (r, c, _) in tests.items():
        f = obs[r]
        o = k.split(":")[0]
        rows.append(dict(test=k, n_fetuses=len(f), cohorts=f.cohort.nunique(), n_covariates=c(f).shape[1],
                         pool_donors=len(pools[o]["depth"]), max_fetuses_per_cohort=int(f.cohort.value_counts().max()),
                         share_capped=float(np.mean(cap(f, pools[o]) < 1))))
    return pd.DataFrame(rows)


def target(f):
    return S["det"].loc[f.key].to_numpy(float) * (f.depth / f.key.map(S["medd"])).to_numpy(float)[:, None]


def cap(f, pool):
    top = pd.DataFrame(target(f)).groupby(f.key.to_numpy()).transform("max").to_numpy()
    return np.minimum(1.0, pool["fmed"] * pool["dref"] / top)


def one(i):
    rng, out = np.random.default_rng([SEED, 1, i]), []
    for j, (name, (rows, cov, strata)) in enumerate(S["tests"].items()):
        f = S["obs"][rows]
        P = S["pools"][name.split(":")[0]]
        pick = np.empty(len(f), int)
        for c in f.cohort.unique():  # distinct donors within a cohort
            m = (f.cohort == c).to_numpy()
            pick[m] = rng.choice(len(P["depth"]), m.sum(), replace=False)
        depth = P["depth"][pick]
        null = np.column_stack([np.round(rng.permutation(c / depth) * depth) for c in P["x"][1.0][pick].T])
        for setting, x in [*((f"a{a:g}", P["x"][a][pick]) for a in SCALES), ("shuffled", null.astype(np.int64))]:
            p_ = cap(f, P) * target(f) / (P["fmed"] * depth[:, None])
            t = f.reset_index(drop=True).copy()
            t[GENES] = rng.binomial(x, np.minimum(p_, 1.0)).astype(float)
            tab = dev.adjust(t, GENES, cov)[0]
            c, p, _ = dev.coherence(tab, WF, N_PERM, np.random.default_rng([SEED, 2, i, j]), strata)
            out.append(dict(test=name, setting=setting, sub=i, C=c, p=p))
    return out


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
    print(design.to_string(index=False), flush=True)
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
                         median_stat=r.C.median(), stat_q05=r.C.quantile(.05)))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"feasibility_passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
