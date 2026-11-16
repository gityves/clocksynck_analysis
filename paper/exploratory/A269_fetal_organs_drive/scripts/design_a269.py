"""Design of A269 (fetal lung and liver of A243: entrained loop or driven response?), shared by the feasibility and the
scoring. The A266 and A264 checks of the fetal adrenal cortex applied to the two organs that passed H243.1.

Rows and base covariates: `<organ>:all` of A243 (`score_a243`: cohorts with >= 2 fetuses, cohort indicator, common
quadratic age). Waveforms:
  WF     the 14-gene GTEx waveform (A243);
  CORE   WF without PER1, PER2 and NR1D1 (the clock genes with a documented direct acute response to glucocorticoid or
         cAMP signalling; as A266), 11 genes.
GC score: per fetus, the mean over GC of the gene's log fraction standardised over the rows of the organ (equal weights;
`design_a264.activity`). GC: four canonical direct glucocorticoid-receptor targets that are not clock genes (FKBP5,
TSC22D3 / GILZ, ZBTB16, KLF9), fixed before any statistic of this amendment.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A243_fetal_organs_pooled" / "scripts"))
sys.path.insert(0, str(HERE.parents[2] / "A264_fetal_adrenal_covariates" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a243 as d243  # noqa: E402
from common import development as dev  # noqa: E402
from design_a264 import activity  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 269
WF = d243.WF
SIGNAL_INDUCED = ("PER1", "PER2", "NR1D1")
CORE = WF.drop(list(SIGNAL_INDUCED))
GC = ("FKBP5", "TSC22D3", "ZBTB16", "KLF9")
ORGANS = ("lung", "liver")
FAMILIES = {"H269.1": ["lung:core", "liver:core"], "H269.2": ["lung:gc", "liver:gc"]}


def sheet(counts: pd.DataFrame) -> pd.DataFrame:
    """The A243 sheet with the GC score (standardised within each organ)."""
    obs = d243.sheet()
    obs["gc"] = np.nan
    for _, g in obs.groupby("organ"):
        obs.loc[g.index, "gc"] = activity(counts.loc[g.index], g.depth, GC)
    return obs


def cov(score: str | None = None):
    """A243 covariates (cohort indicator, common quadratic age), optionally with a score column."""
    def f(g):
        Z = [dev.design(g, "age", 2, ("block",))]
        if score:
            Z.append(g[[score]].to_numpy(float))
        return np.column_stack(Z)
    return f


def tests(obs) -> dict:
    """{name: (rows, covariates, strata, waveform)}: H269.1, H269.2 and the descriptive readouts."""
    t = {}
    for o in ORGANS:
        rows = d243.pooled(obs, obs.organ == o)
        t[f"{o}:core"] = (rows, cov(), "block", CORE)
        t[f"{o}:gc"] = (rows, cov("gc"), "block", WF)
        t[f"{o}:core:gc"] = (rows, cov("gc"), "block", CORE)
    t["lung:early:core"] = (d243.pooled(obs, (obs.organ == "lung") & (obs.age < d243.LATE)), cov(), "block", CORE)
    return t
