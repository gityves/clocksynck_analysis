"""A243 scoring (registered in `paper/docs/registered_experiments.md`, with its design addendum): is a clock of the size read
in the fetal adrenal present in the heart, kidney, lung or liver before birth when open cohorts are pooled?

One whole-organ pseudobulk per fetus, organ and cohort (`prep_a243.py`; eligible rows of `samples.csv`). L1 signed coherence
C across fetuses (`common.development`), pooled over cohorts with >= 2 fetuses in the readout: log fractions residualised
gene by gene on a cohort indicator and a common age polynomial; per-gene permutations within cohort; matched and responsive
companions among the genes present in every cohort of the readout. As A242.

Readouts per organ: `<organ>:all` (quadratic age), `<organ>:late` (>= 12 weeks, linear), `<organ>:early` (< 12 weeks,
quadratic); descriptive `<organ>:all:new` (without cao2020 and cardoso_moreira), `<organ>:all:-<cohort>`, `<organ>:<cohort>`
(cohorts with >= 10 fetuses, quadratic age). Only readouts with an adult pass rate >= 0.8 at the full clock are scored
(`feasibility_a243.py`); the others are written as not estimable.

    python score_a243.py [--jobs 8] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]
Outputs (`../output/`): design.csv, l1.csv, verdicts.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SRC = DATA / "natexp" / "fetal_organs_pooled"
SEED = 243
WF = cs.gtex()
ORGANS = ("heart", "kidney", "lung", "liver")
OLD = ("cao2020", "cardoso_moreira")  # scored alone before (A189 / A241, A237)
LATE, MIN_FETUSES, MIN_ALONE = 12.0, 2, 10


def sheet() -> pd.DataFrame:
    s = pd.read_csv(SRC / "samples.csv")
    s = s[s.eligible].set_index("id")
    return pd.DataFrame(dict(donor=s.fetus.astype(str) + "@" + s.cohort, compartment=s.organ, organ=s.organ,
                             cohort=s.cohort, block=s.cohort, age=s.age_weeks.astype(float), depth=s.depth.astype(float),
                             stratum="all"))


def load():
    obs = sheet()
    counts = pd.read_parquet(SRC / "pseudobulk.parquet").loc[obs.index]
    return counts, obs


def cov(df: int):
    return lambda g: dev.design(g, "age" if df else None, df, ("block",))


def pooled(obs, rows, without=()):
    rows = rows & ~obs.cohort.isin(list(without))
    n = obs[rows].groupby("cohort").donor.transform("size").reindex(obs.index, fill_value=0)
    return rows & (n >= MIN_FETUSES)


def all_tests(obs) -> dict:
    l1 = {}
    for o in ORGANS:
        org = obs.organ == o
        l1[f"{o}:all"] = (pooled(obs, org), cov(2), "block")
        l1[f"{o}:late"] = (pooled(obs, org & (obs.age >= LATE)), cov(1), "block")
        l1[f"{o}:early"] = (pooled(obs, org & (obs.age < LATE)), cov(2), "block")
        l1[f"{o}:all:new"] = (pooled(obs, org, OLD), cov(2), "block")
        cohorts = sorted(obs.cohort[l1[f"{o}:all"][0]].unique())
        for c in cohorts:
            l1[f"{o}:all:-{c}"] = (pooled(obs, org, [c]), cov(2), "block")
            alone = org & (obs.cohort == c)
            if alone.sum() >= MIN_ALONE:
                l1[f"{o}:{c}"] = (alone, cov(2), "block")
    return {k: v for k, v in l1.items() if v[0].sum() >= 3}


def families(scored) -> dict:
    return {"H243.1": [f"{o}:all" for o in ORGANS if f"{o}:all" in scored],
            "H243.2": [f"{o}:{w}" for o in ORGANS for w in ("late", "early") if f"{o}:{w}" in scored]}


def tests(obs):
    rates = OUT / "feasibility_passrates.csv"
    l1 = all_tests(obs)
    if rates.exists():
        r = pd.read_csv(rates)
        a1 = r[r.setting == "a1"].set_index("test").pass_rate
        power_of = lambda k: a1.get(k, a1.get(k.split(":-")[0], 0.0))  # leave-one-cohort-out: power of `<organ>:all`
        l1 = {k: v for k, v in l1.items() if power_of(k) >= 0.8}
    return l1, {}


if __name__ == "__main__":
    counts, obs = load()
    scored = list(tests(obs)[0])
    dev.main(__doc__, load, tests, WF, SEED, OUT, families(scored))
