"""Adult power of the A269 readouts (fetal lung and liver: core loop, GC adjustment), before any fetal statistic.

As `A243/scripts/feasibility_a243.py`: every fetus of a readout is stood for by one GTEx donor of the adult organ (Lung,
Liver; deaths at 08-18 h, any Hardy class; distinct donors within a cohort), its clock-gene counts thinned gene by gene to
what the fetus's cohort captures; adult clock at 1, 0.75, 0.5 (`simulate.scale_clock`) and shuffled (false-positive
rate); 200 subsamples, 500 permutations within cohort. New here (as A264): the stand-in also receives the adult donor's
own GC score (`design_a269.GC`, `activity` over the adult pool of the organ, unthinned), so that the GC adjustment removes
from the adult clock whatever the adult score shares with it; each readout is read on its own waveform (CORE or WF).

Also written (aggregates of the adult pool, no donor row): the share of the variance of the adult GC score explained by
the hour of death (cosinor R^2, 08-18 h) and the mean pairwise correlation of the GC genes, per organ.
Controlled access: GTEx time of death is used in memory only; outputs are aggregates.

    python feasibility_a269.py [--n-sub 200] [--jobs 8] [--check]
Outputs (`../output/`): feasibility_design.csv, feasibility_subsamples.csv, feasibility_passrates.csv, adult_scores.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import design_a269 as d  # noqa: E402

sys.path.insert(0, str(HERE.parents[2] / "A242_fetal_adrenal_replication" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import feasibility_a243 as F  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from design_a264 import activity  # noqa: E402
from feasibility_a242 import ci  # noqa: E402

OUT, SEED, N_PERM = d.OUT, d.SEED, F.N_PERM
GENES, SCALES = F.GENES, F.SCALES
TISSUE = {o: F.TISSUE[o] for o in d.ORGANS}
S = F.S


def setup():
    counts = pd.read_parquet(d.d243.SRC / "pseudobulk.parquet", columns=list(d.GC)).loc[d.d243.sheet().index]
    obs = d.sheet(counts)
    pb, don = gtex.load(list(TISSUE.values()), hardy=(0, 1, 2, 3, 4), genes=[*GENES, *d.GC], depth="library_size")
    keep = in_window(don.tod.to_numpy())
    pb = pb[pb.donor.isin(set(don.index[keep]))]
    rng = np.random.default_rng([SEED, 0])
    pools, scores = {}, []
    for o, t in TISSUE.items():
        x = pb[pb.compartment == t].reset_index(drop=True)
        depth = x.depth.to_numpy(float)
        gc = activity(x, depth, d.GC)
        h = don.tod.reindex(x.donor).to_numpy(float) * 2 * np.pi / 24
        X = np.column_stack([np.ones(len(h)), np.cos(h), np.sin(h)])
        res = gc - X @ np.linalg.lstsq(X, gc, rcond=None)[0]
        r = np.corrcoef(np.log((x[list(d.GC)].to_numpy(float) + 0.5) / depth[:, None]).T)
        scores.append(dict(organ=o, donors=len(x), gc_r2_hour=float(1 - res.var() / gc.var()),
                           gc_mean_pairwise_r=float(r[np.triu_indices(len(d.GC), 1)].mean())))
        pools[o] = dict(depth=depth, gc=gc, fmed=np.median(x[GENES].to_numpy(float) / depth[:, None], axis=0),
                        dref=depth.min(),
                        x={a: scale_clock(x[["donor", "compartment", "depth", *GENES]], d.WF, a, depth=x.depth, rng=rng)[GENES]
                           .to_numpy(np.int64) for a in SCALES})
    pd.DataFrame(scores).to_csv(OUT / "adult_scores.csv", index=False)
    print(pd.DataFrame(scores).round(4).to_string(index=False), flush=True)
    det = pd.read_csv(d.d243.OUT / "clock_detection.csv")
    det["key"] = det.organ + "|" + det.cohort
    obs["key"] = obs.organ + "|" + obs.cohort
    tests = d.tests(obs)
    S.update(obs=obs, tests=tests, pools=pools, det=det.set_index("key")[GENES], medd=obs.groupby("key").depth.median())
    rows = []
    for k, (r, c, _, wf) in tests.items():
        f = obs[r]
        o = k.split(":")[0]
        rows.append(dict(test=k, n_fetuses=len(f), cohorts=f.cohort.nunique(), n_covariates=c(f).shape[1],
                         n_genes=len(wf.genes), pool_donors=len(pools[o]["depth"]),
                         share_capped=float(np.mean(F.cap(f, pools[o]) < 1))))
    return pd.DataFrame(rows)


def one(i):
    rng, out = np.random.default_rng([SEED, 1, i]), []
    for j, (name, (rows, cov, strata, wf)) in enumerate(S["tests"].items()):
        f = S["obs"][rows]
        P = S["pools"][name.split(":")[0]]
        pick = np.empty(len(f), int)
        for c in f.cohort.unique():  # distinct donors within a cohort
            m = (f.cohort == c).to_numpy()
            pick[m] = rng.choice(len(P["depth"]), m.sum(), replace=False)
        depth = P["depth"][pick]
        null = np.column_stack([np.round(rng.permutation(c / depth) * depth) for c in P["x"][1.0][pick].T])
        genes = list(wf.genes)
        for setting, x in [*((f"a{a:g}", P["x"][a][pick]) for a in SCALES), ("shuffled", null.astype(np.int64))]:
            p_ = F.cap(f, P) * F.target(f) / (P["fmed"] * depth[:, None])
            t = f.reset_index(drop=True).copy()
            t[GENES] = rng.binomial(x, np.minimum(p_, 1.0)).astype(float)
            t["gc"] = P["gc"][pick]
            tab = dev.adjust(t, genes, cov)[0]
            c, p, _ = dev.coherence(tab, wf, N_PERM, np.random.default_rng([SEED, 2, i, j]), strata)
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
    OUT.mkdir(exist_ok=True)
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
