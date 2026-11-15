"""Adult power of the A264 readouts (the fetal adrenal cortex under three further adjustments), before any fetal statistic.

As `A244/scripts/feasibility_a244.py` (L1 part): every fetus of the cortex readout is stood for by one distinct GTEx adrenal
gland (death at 08-18 h, any Hardy class) that receives its cohort and age; the adult clock-gene counts are thinned gene by
gene to what the cortex of the cohort captures (`feasibility_a244.stand_in`). New here: the stand-in also receives the
adult donor's own activity scores (`design_a264.activity` over the adult pool, unthinned), so that an adjustment for
steroidogenic activity removes from the adult clock whatever the adult score shares with it. Settings a1, a0.75, a0.5
(adult clock scaled) and shuffled (clock genes permuted across the drawn donors: the false-positive rate).

Also written (aggregates of the adult pool, no donor row): the mean pairwise correlation of the STEROID genes, and the
share of the variance of each score explained by the hour of death (cosinor, R^2) in the 08-18 h window.
Controlled access: GTEx time of death is used in memory only.

    python feasibility_a264.py [--n-sub 200] [--jobs 8] [--check]

Outputs (`../output/`): feasibility_design.csv, feasibility_subsamples.csv, feasibility_passrates.csv,
adult_scores.csv (`_check` with --check).
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2] / "A244_fetal_adrenal_classes" / "scripts"))
sys.path.insert(0, str(HERE.parents[2] / "A242_fetal_adrenal_replication" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import design_a264 as d  # noqa: E402
import feasibility_a244 as F  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from feasibility_a242 import ci  # noqa: E402

OUT, WF, SEED, N_PERM = d.OUT, d.WF, d.SEED, 500
GENES = list(WF.genes)
SCALES = (1.0, 0.75, 0.5)
W = 2 * np.pi / 24


def setup():
    import score_a264 as S264
    counts, obs = S264.load([*GENES, *d.STEROID, *d.ACUTE])  # the sheet only: no fetal statistic
    pb, don = gtex.load(["Adrenal Gland"], hardy=(0, 1, 2, 3, 4), genes=[*GENES, *d.STEROID, *d.ACUTE],
                        depth="library_size")
    pb = pb[pb.donor.isin(set(don.index[in_window(don.tod.to_numpy())]))].sort_values("donor").reset_index(drop=True)
    depth, rng = pb.depth.to_numpy(float), np.random.default_rng([SEED, 0])
    clock = pb[["donor", "compartment", "depth", *GENES]]
    pool = dict(depth=depth, fmed=np.median(clock[GENES].to_numpy(float) / depth[:, None], axis=0), dref=depth.min(),
                x={a: scale_clock(clock, WF, a, depth=clock.depth, rng=rng)[GENES].to_numpy(np.int64) for a in SCALES},
                steroid=d.activity(pb, depth, d.STEROID), acute=d.activity(pb, depth, d.ACUTE))
    # aggregates of the adult pool
    y = np.log((pb[list(d.STEROID)].to_numpy(float) + 0.5) / depth[:, None])
    r = np.corrcoef(y.T)
    tod = don.tod.reindex(pb.donor).to_numpy(float)
    X = np.column_stack([np.ones(len(tod)), np.cos(W * tod), np.sin(W * tod)])
    agg = [dict(quantity="steroid_mean_pairwise_r", value=float(r[np.triu_indices_from(r, 1)].mean()), n=len(pb))]
    for k in ("steroid", "acute"):
        s = pool[k]
        res = s - X @ np.linalg.lstsq(X, s, rcond=None)[0]
        agg.append(dict(quantity=f"{k}_R2_hour_of_death", value=float(1 - res.var() / s.var()), n=len(pb)))
    det = pd.read_csv(HERE.parents[2] / "A244_fetal_adrenal_classes" / "output" / "clock_detection_classes.csv")
    det["key"] = det.cohort + "|" + det.compartment
    tests = {k: v for k, v in d.tests(obs).items()}
    F.S.update(det=det.set_index("key")[GENES], medd=obs.groupby("key").depth.median())
    F.S.update(obs=obs, tests=tests, pool=pool)
    rows = []
    for k, (r_, c, _, wf) in tests.items():
        f = obs[r_]
        rows.append(dict(test=k, n_fetuses=f.donor.nunique(), cohorts=f.cohort.nunique(), n_covariates=c(f).shape[1],
                         n_genes=len(wf.genes), pool_donors=len(depth), share_capped=np.mean(F.cap(f, pool) < 1)))
    return pd.DataFrame(rows), pd.DataFrame(agg)


def one(i):
    rng, out, P = np.random.default_rng([SEED, 1, i]), [], F.S["pool"]
    for j, (name, (rows, cov, strata, wf)) in enumerate(F.S["tests"].items()):
        f = F.S["obs"][rows]
        pick = rng.choice(len(P["depth"]), len(f), replace=False)
        depth = P["depth"][pick]
        null = np.column_stack([np.round(rng.permutation(c / depth) * depth) for c in P["x"][1.0][pick].T])
        for setting, x in [*((f"a{a:g}", P["x"][a][pick]) for a in SCALES), ("shuffled", null.astype(np.int64))]:
            t = F.stand_in(f, x, depth, P, rng)
            t["steroid"], t["acute"] = P["steroid"][pick], P["acute"][pick]
            tab = dev.adjust(t, GENES, cov)[0]
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
    pd.set_option("display.width", 250)
    design, agg = setup()
    design.to_csv(OUT / f"feasibility_design{sfx}.csv", index=False)
    agg.to_csv(OUT / f"adult_scores{sfx}.csv", index=False)
    print(design.round(3).to_string(index=False), agg.round(3).to_string(index=False), sep="\n", flush=True)
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
