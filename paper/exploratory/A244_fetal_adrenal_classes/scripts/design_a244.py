"""Design of the A244 readouts (shared by the feasibility and the scoring): the class pseudobulk sheet and the tests.

`sheet()`: one row per fetus and class (>= 200 cells; `prep_a244.py`), with the columns `common.development` expects
(donor = fetus, compartment = class, block = cohort, stratum = cohort, age in weeks as stated by the source, depth).
`tests(obs)`: the L1 readouts {name: (rows, covariates, strata)} of A242's `all`, per class; the windows
`<class>:late` (from 12 weeks) and `<class>:early` as in A242 (descriptive); and the rows of the two-class gate
(`two_classes`: fetuses with both a cortex and a medulla pseudobulk). Cohorts enter with >= 2 fetuses.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A242_fetal_adrenal_replication" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a242 as a242  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "natexp" / "fetal_adrenal_open"
OUT = HERE.parents[1] / "output"
SEED = 244
WF = a242.WF
CLASSES = ["cortex", "medulla"]  # scored; stroma is descriptive
LATE = a242.LATE


def sheet() -> pd.DataFrame:
    s = pd.read_csv(SRC / "a244_class_samples.csv").set_index("id")
    return pd.DataFrame(dict(donor=s.donor, compartment=s.compartment, cohort=s.cohort, block=s.cohort, stratum=s.cohort,
                             age=s.age, depth=s.depth, cells=s.cells, key=s.cohort + "|" + s.compartment))


def pooled(obs, rows):
    """`rows` restricted to the cohorts with >= 2 of them (a single fetus has no residual)."""
    n = obs[rows].groupby("cohort").donor.transform("size").reindex(obs.index, fill_value=0)
    return rows & (n >= a242.MIN_FETUSES) & obs.cohort.isin(a242.POOLED)


def two_classes(obs):
    """Rows of the fetuses that have both a cortex and a medulla pseudobulk, in cohorts with >= 2 such fetuses."""
    both = obs.groupby("donor").compartment.transform(lambda c: set(CLASSES) <= set(c))
    rows = obs.compartment.isin(CLASSES) & both
    n = obs[rows].groupby("cohort").donor.nunique()
    return rows & obs.cohort.map(n).fillna(0).ge(a242.MIN_FETUSES)


def tests(obs):
    """L1 readouts and the L2 readout (rows, covariates, strata) of A244."""
    l1 = {}
    for k in [*CLASSES, "stroma"]:
        c = obs.compartment == k
        l1[k] = (pooled(obs, c), a242.cov(2), "block")
        l1[f"{k}:late"] = (pooled(obs, c & (obs.age >= LATE)), a242.cov(1), "block")
        l1[f"{k}:early"] = (pooled(obs, c & (obs.age < LATE)), a242.cov(2), "block")
    return l1, {"two_classes": (two_classes(obs), a242.cov(2))}
