"""Adult power calibration of A244 (cortex and medulla of the fetal adrenal), before any fetal clock statistic.

Read from the fetal data: the class sample sheet (cohort, age, depth, cells) and the cohort medians of the clock-gene
counts of each class (`../output/clock_detection_classes.csv`). No fetal clock statistic is computed.

L1 (`cortex`, `medulla`; `stroma` and the windows `<class>:late` / `:early` descriptive): as A242's feasibility. Every
fetus is stood for by one distinct GTEx adrenal gland (death at 08-18 h, any Hardy class) that receives its cohort,
age and depth; the adult clock-gene counts are thinned binomially, gene by gene, to what the class of the cohort
captures (median count of the gene in the class pseudobulks of the cohort x depth of the fetus / median depth of those
pseudobulks; where the adult library holds fewer reads the targets of the gene are lowered over the class and cohort,
`cap`). Statistic: `development.coherence` with cohort indicator and age trend, per-gene permutations within cohort.
Settings: a1, a0.75, a0.5 (adult clock scaled, `simulate.scale_clock`), shuffled (every clock gene permuted across the
drawn donors: the false-positive rate).

L2 (`two_classes`, H244.2): each fetus with both pseudobulks is stood for by one GTEx donor with an adrenal gland
(the cortex) and `SECOND` (the medulla), both thinned to the counts of their class and cohort; `development.shared`
on the cohort- and age-adjusted counts with the cohorts as strata (pairing permuted within cohort); pass = p_pairing
and p_acrophase < 0.05. Settings a1, a0.75, a0.5 and `mispaired` (adult clock at 1, the second tissue taken from
another drawn donor: the false-positive rate, whose 95 % CI must contain 0.05).
Controlled access: GTEx time of death is used in memory only; outputs are aggregates, without donor identifiers.

    python feasibility_a244.py [--n-sub 200] [--jobs 8] [--check]

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

import design_a244 as d  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from feasibility_a242 import ci  # noqa: E402

OUT, WF, SEED, N_PERM = d.OUT, d.WF, d.SEED, 500
GENES = list(WF.genes)
SCALES = (1.0, 0.75, 0.5)
SECOND = "Muscle - Skeletal"  # second adult tissue of H244.2 (the reason is in the design addendum)
S = {}  # state inherited by the workers


def pool_of(table):
    """Adult pool of one tissue: depths, clock-gene counts at every scale, median fractions, smallest depth."""
    depth, rng = table.depth.to_numpy(float), np.random.default_rng([SEED, 0])
    clock = table[["donor", "compartment", "depth", *GENES]]
    return dict(depth=depth, fmed=np.median(clock[GENES].to_numpy(float) / depth[:, None], axis=0), dref=depth.min(),
                x={a: scale_clock(clock, WF, a, depth=clock.depth, rng=rng)[GENES].to_numpy(np.int64) for a in SCALES})


def setup():
    obs = d.sheet()
    pb, don = gtex.load(["Adrenal Gland", SECOND], hardy=(0, 1, 2, 3, 4), genes="clock", depth="library_size")
    pb = pb[pb.donor.isin(set(don.index[in_window(don.tod.to_numpy())]))]
    both = set(pb[pb.compartment == "Adrenal Gland"].donor) & set(pb[pb.compartment == SECOND].donor)
    pools = {}
    for name, tissue, donors in (("adrenal", "Adrenal Gland", None), ("adrenal2", "Adrenal Gland", both),
                                 ("second", SECOND, both)):
        t = pb[pb.compartment == tissue]
        t = (t[t.donor.isin(donors)].sort_values("donor") if donors is not None else t).reset_index(drop=True)
        pools[name] = pool_of(t)
    det = pd.read_csv(OUT / "clock_detection_classes.csv")
    det["key"] = det.cohort + "|" + det.compartment
    medd = obs.groupby("key").depth.median()
    l1, l2 = d.tests(obs)
    l1 = {k: v for k, v in l1.items() if v[0].sum() >= 3}
    S.update(obs=obs, l1=l1, l2=l2, pools=pools, det=det.set_index("key")[GENES], medd=medd,
             n_both=pools["second"]["depth"].size)
    rows = []
    for k, (r, c, _) in l1.items():
        f = obs[r]
        rows.append(dict(test=k, level="L1", n_fetuses=f.donor.nunique(), n_samples=len(f), cohorts=f.cohort.nunique(),
                         n_covariates=c(f).shape[1], pool_donors=len(pools["adrenal"]["depth"]),
                         share_capped=np.mean(cap(f, pools["adrenal"]) < 1)))
    f = obs[l2["two_classes"][0]]
    rows.append(dict(test="two_classes", level="L2", n_fetuses=f.donor.nunique(), n_samples=len(f),
                     cohorts=f.cohort.nunique(), n_covariates=l2["two_classes"][1](f[f.compartment == "cortex"]).shape[1],
                     pool_donors=S["n_both"],
                     share_capped=np.mean(np.r_[cap(f[f.compartment == "cortex"], pools["adrenal2"]),
                                                cap(f[f.compartment == "medulla"], pools["second"])] < 1)))
    return pd.DataFrame(rows)


def target(f):
    """Clock-gene counts the class of the cohort would capture in each fetus of `f` (fetuses x genes)."""
    return S["det"].loc[f.key].to_numpy(float) * (f.depth / f.key.map(S["medd"])).to_numpy(float)[:, None]


def cap(f, pool):
    """Factor (<= 1) on the targets of each gene, one per class and cohort: adult libraries cannot be thinned upwards."""
    top = pd.DataFrame(target(f)).groupby(f.key.to_numpy()).transform("max").to_numpy()
    return np.minimum(1.0, pool["fmed"] * pool["dref"] / top)


def stand_in(f, x, depth, pool, rng):
    """Table of the adult counts `x` (library sizes `depth`) thinned to the design `f`."""
    p = cap(f, pool) * target(f) / (pool["fmed"] * depth[:, None])
    t = f.reset_index(drop=True).copy()
    t[GENES] = rng.binomial(x, np.minimum(p, 1.0)).astype(float)
    return t


def one(i):
    rng, out = np.random.default_rng([SEED, 1, i]), []
    P = S["pools"]
    for j, (name, (rows, cov, strata)) in enumerate(S["l1"].items()):
        f = S["obs"][rows]
        pick = rng.choice(len(P["adrenal"]["depth"]), len(f), replace=False)
        depth = P["adrenal"]["depth"][pick]
        null = np.column_stack([np.round(rng.permutation(c / depth) * depth) for c in P["adrenal"]["x"][1.0][pick].T])
        for setting, x in [*((f"a{a:g}", P["adrenal"]["x"][a][pick]) for a in SCALES), ("shuffled", null.astype(np.int64))]:
            tab = dev.adjust(stand_in(f, x, depth, P["adrenal"], rng), GENES, cov)[0]
            c, p, _ = dev.coherence(tab, WF, N_PERM, np.random.default_rng([SEED, 2, i, j]), strata)
            out.append(dict(test=name, setting=setting, sub=i, C=c, p=p, p_pairing=np.nan, p_acrophase=np.nan))
    rows, cov = S["l2"]["two_classes"]
    f = S["obs"][rows].sort_values(["donor", "compartment"])
    fe = sorted(f.donor.unique())
    pick = rng.choice(S["n_both"], len(fe), replace=False)
    other = np.roll(rng.permutation(pick), 1 + rng.integers(len(fe) - 1))  # no donor keeps its own second tissue
    for setting in [*(f"a{a:g}" for a in SCALES), "mispaired"]:
        a = 1.0 if setting == "mispaired" else float(setting[1:])
        parts = []
        for cls, pool, idx in (("cortex", "adrenal2", pick), ("medulla", "second", other if setting == "mispaired" else pick)):
            g = f[f.compartment == cls].set_index("donor").loc[fe].reset_index()
            t = stand_in(g, P[pool]["x"][a][idx], P[pool]["depth"][idx], P[pool], rng)
            t["donor"] = [f"s{k}" for k in range(len(fe))]
            parts.append(t)
        tab = pd.concat(parts, ignore_index=True)
        tab = dev.adjust(tab, GENES, cov)[0]
        g = dev.shared(tab, WF, N_PERM, SEED + i)
        out.append(dict(test="two_classes", setting=setting, sub=i, C=g.pairing_z if g else np.nan,
                        p=max(g.p_pairing, g.p_acrophase) if g else np.nan, p_pairing=g.p_pairing if g else np.nan,
                        p_acrophase=g.p_acrophase if g else np.nan))
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
        rows.append(dict(test=test, setting=setting, n_subsamples=m, n_estimable=int(r.p.notna().sum()), pass_rate=k / m,
                         ci_lo=lo, ci_hi=hi, median_stat=r.C.median(), stat_q05=r.C.quantile(.05)))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"feasibility_passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
