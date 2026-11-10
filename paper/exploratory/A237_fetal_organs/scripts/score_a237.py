"""A237 scoring (registered in `paper/docs/registered_experiments.md`): a clock in human prenatal organs, 4-20 weeks
post conception (Cardoso-Moreira 2019, E-MTAB-6814, bulk RNA-seq).

L1 per organ and L2 across the organs of one individual, on log fractions residualised gene by gene within each organ
on age (polynomial) and sex, with expression-matched and responsive pseudo-clock companions (`common.development`).
Counts are recovered from the released CPM (CPM / smallest non-zero CPM of the library); one row per gene symbol (the
Ensembl gene with the largest mean count). Individuals are those of `prep_fetal_bulk/scripts/link_cardoso_individuals.py`
(`organs`: linked libraries, without the one whose sex differs from that of its individual; `organs_all`: every
library, the unlinked ones at their guess).

    python score_a237.py [--jobs 12] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design.csv, l1.csv, l2.csv, l2_pseudo.csv, verdicts.csv (`_check` with --check: 200
permutations, 10 sets).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SRC = DATA / "natexp" / "cardoso_moreira2019"
SEED = 237
WF = cs.gtex()
ORGANS = ["forebrain", "hindbrain", "heart", "kidney", "liver", "ovary", "testis"]
RANGES = {"4-20": (4, 20), "9-20": (9, 20), "13-20": (13, 20), "4-8": (4, 8)}  # weeks post conception
POSTNATAL = ["newborn", "infant", "toddler", "school", "youngTeenager", "teenager", "oldTeenager", "youngAdult",
             "youngMidAge", "olderMidAge", "senior"]
CONTROLS = ["forebrain", "hindbrain", "heart", "kidney", "liver"]  # postnatal L1 with an adult pass rate >= 0.8
FAMILIES = {"H237.1": ["organs:9-20wpc"],  # Holm within a hypothesis
            "H237.2": [f"{o}:4-20wpc" for o in ("heart", "liver", "forebrain", "hindbrain", "kidney", "ovary")],
            "H237.3": [f"{o}:4-8wpc" for o in ("heart", "liver", "hindbrain", "kidney")]
                      + [f"{o}:13-20wpc" for o in ("heart", "liver", "kidney")],
            **{f"control {o}": [f"{o}:postnatal"] for o in CONTROLS}}


def load():
    """Counts (libraries x gene symbols) and the library sheet (donor, compartment, depth, age, sex, stratum)."""
    s = pd.read_csv(SRC / "individuals.csv").set_index("sample")
    cpm = pd.read_parquet(SRC / "expr_cpm.parquet")[s.index]
    x = cpm.to_numpy(float)
    x = x / np.where(x > 0, x, np.inf).min(0)  # counts: CPM / smallest non-zero CPM of the library
    sym = pd.read_csv(SRC / "genes.csv").set_index("gene_id").symbol.reindex(cpm.index).dropna().map(canonical)
    top = pd.Series(x.mean(1), index=cpm.index)[sym.index].groupby(sym).idxmax()
    counts = pd.DataFrame(x[cpm.index.get_indexer(top.to_numpy())].T, index=s.index, columns=top.index)
    linked = s.prenatal & (s.link != "open")
    major = s[linked].groupby("individual").sex.agg(lambda x: x.mode().iloc[0])
    linked &= s.sex == s.individual.map(major)  # an individual has one sex: a library of the other sex is not linked
    obs = pd.DataFrame(dict(donor=s.individual.map(lambda v: f"i{v:.0f}"), compartment=s.organ,
                            depth=counts.sum(axis=1), prenatal=s.prenatal, linked=linked, sex=s.sex, rqn=s.rqn,
                            age=np.where(s.prenatal, s.age_wpc, s.stage.str.lower().map(
                                {k.lower(): i for i, k in enumerate(POSTNATAL)})),
                            stratum=s.age_bin.astype(str)))
    obs.loc[~obs.prenatal, "donor"] = obs.index[~obs.prenatal]  # postnatal individuals are not used
    return counts, obs


def cov(df, rqn=False):
    """Design of a set of libraries: polynomial of age (df), sex, and the RNA quality number when asked."""
    def f(g):
        Z = dev.design(g, "age", df, ("sex",))
        return np.column_stack([Z, g.rqn - g.rqn.mean()]) if rqn else Z
    return f


def tests(obs):
    """L1 and L2 tests {name: (rows, covariates)}; names are `<readout>:<range>` as in the registration."""
    pre, l1, l2 = obs.prenatal, {}, {}
    for o in ORGANS:
        for tag, (lo, hi) in RANGES.items():
            rows = pre & (obs.compartment == o) & obs.age.between(lo, hi)
            if rows.sum() >= 8:
                l1[f"{o}:{tag}wpc"] = (rows, cov(3 if hi - lo > 8 else 2))
        rows = pre & (obs.compartment == o) & obs.rqn.notna()
        l1[f"{o}:4-20wpc:rqn"] = (rows, cov(3, rqn=True))
        rows = ~pre & (obs.compartment == o)
        if rows.sum() >= 8:
            l1[f"{o}:postnatal"] = (rows, cov(2))
    for name, base in (("organs", pre & obs.linked), ("organs_all", pre)):
        for tag, (lo, hi) in RANGES.items():
            rows = base & obs.age.between(lo, hi)
            rows &= obs[rows].groupby("donor").compartment.transform("nunique").reindex(obs.index, fill_value=0) >= 2
            l2[f"{name}:{tag}wpc"] = (rows, cov(3 if hi - lo > 8 else 2))
    rows = l2["organs:9-20wpc"][0] & obs.rqn.notna()
    l2["organs:9-20wpc:rqn"] = (rows, cov(3, rqn=True))
    return l1, l2


if __name__ == "__main__":
    dev.main(__doc__, load, tests, WF, SEED, OUT, FAMILIES)
