"""A264 scoring (registered in `paper/docs/registered_experiments.md`; NOT BLIND): does the clock-shaped covariation of
the fetal adrenal cortex (A244) survive (i) further age and cohort adjustment, (ii) adjustment for steroidogenic activity,
(iii) removal of PER1 and NR1D1?

Readouts and covariates: `design_a264.tests`. Statistic as A244: L1 signed coherence C (`common.development.l1_test`),
10,000 per-gene permutations within cohort, 200 expression-matched and 200 responsive companion sets adjusted in the same
way (and on the same waveform). A readout is scored only if its adult pass rate at the full clock is >= 0.8
(`feasibility_a264.py`); otherwise it is written as "not estimable" and not computed.

    python score_a264.py [--jobs 6] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design.csv, l1.csv, verdicts.csv (`_check` with --check: 200 permutations, 10 sets).
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
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import design_a264 as d  # noqa: E402
import prep_a244 as PR  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import development as dev  # noqa: E402

OUT = d.OUT
S = {}


def load(genes=None):
    """Class pseudobulks (all genes, or `genes`) and the A244 sheet with the activity scores."""
    pb = pd.read_parquet(PR.PB, columns=None if genes is None else list(dict.fromkeys(genes)))
    obs = d.sheet(pb.loc[d.d244.sheet().index])
    return pb.loc[obs.index], obs


def task(name):
    rows, cov, strata, wf = S["tests"][name]
    j = list(S["tests"]).index(name)
    return name, dev.l1_test(S["counts"], S["obs"][rows], cov, wf, S["n_perm"], S["n_pseudo"], [d.SEED, 1, j], strata)


def verdict(r, adult):
    if adult < 0.8:
        return "not estimable"
    if r.p < 0.05:
        return "PASS" if r.clock_specific else "not clock-specific"
    return "fail"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--design", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    t0, sfx = time.time(), "_check" if a.check else ""
    counts, obs = load()
    tests = d.tests(obs)
    design = pd.DataFrame([dict(test=k, n=int(r.sum()), cohorts=obs[r].cohort.nunique(), n_covariates=c(obs[r]).shape[1],
                                n_genes=len(wf.genes)) for k, (r, c, _, wf) in tests.items()])
    design.to_csv(OUT / "design.csv", index=False)
    print(design.to_string(index=False), flush=True)
    if a.design:
        return
    rates = pd.read_csv(OUT / "feasibility_passrates.csv")
    adult = rates[rates.setting == "a1"].set_index("test").pass_rate
    run = [k for k in tests if adult.get(k, 0.0) >= 0.8]
    S.update(counts=counts, obs=obs, tests=tests, n_perm=200 if a.check else a.n_perm, n_pseudo=10 if a.check else a.n_pseudo)
    res = {}
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for name, r in pool.imap_unordered(task, run):
            res[name] = r
            print(f"{name} done {time.time() - t0:.0f}s", flush=True)
    l1 = pd.DataFrame([dict(test=k, adult_pass_a1=adult.get(k, np.nan), **res.get(k, {})) for k in tests])
    l1.to_csv(OUT / f"l1{sfx}.csv", index=False)
    v = pd.DataFrame([dict(test=k, adult_pass_a1=adult.get(k, np.nan), n=r.get("n"), C=r.get("C"), p=r.get("p"),
                           matched_pct=r.get("matched_pct"), responsive_pct=r.get("responsive_pct"),
                           registered=k in d.REGISTERED,
                           verdict=verdict(pd.Series(r), adult.get(k, 0.0)) if r else "not estimable")
                      for k in tests for r in [res.get(k, {})]])
    v.to_csv(OUT / f"verdicts{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(l1.round(4).to_string(index=False), v.round(4).to_string(index=False), f"runtime {time.time() - t0:.0f}s",
          sep="\n", flush=True)


if __name__ == "__main__":
    main()
