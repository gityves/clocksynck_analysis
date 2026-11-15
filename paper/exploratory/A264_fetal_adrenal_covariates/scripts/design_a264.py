"""Design of the A264 readouts (shared by the feasibility and the scoring): the fetal adrenal cortex of A244 under three
further adjustments.

Rows: the `cortex` readout of A244 (class pseudobulks with >= 200 cells, cohorts with >= 2 fetuses; `design_a244`).
Adjustments, each on top of the A244 design (cohort indicator, common quadratic age):
  (i)   `cortex:age`      common cubic age and a linear age slope of its own for every cohort with >= SLOPE_MIN fetuses
                          in the readout (cohort x age);
  (ii)  `cortex:steroid`  a steroidogenic activity score: per fetus, the mean over STEROID of the gene's log fraction
                          standardised over the rows of the compartment (equal weights);
  (iii) `cortex:-PER1-NR1D1`  the waveform without PER1 and NR1D1 (12 clock genes).
Descriptive: `cortex` (the A244 design, reproduced), `cortex:all3` ((i) to (iii) together), `cortex:acute` (the A244
design plus the score of ACUTE in place of STEROID).

STEROID: the 20 genes whose expression rose >= 4-fold after 48 h of ACTH in both adult and fetal primary human
adrenocortical cells (Xing, Parker, Edwards and Rainey 2010, J Mol Endocrinol 45:59, Table 3; PMC3725639), fixed before
any statistic of this amendment; none is a clock gene. ACUTE: the NR4A family (immediate-early cAMP targets; the
candidate NR4A1 of the prompt, with its two paralogues), descriptive only.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A244_fetal_adrenal_classes" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import design_a244 as d244  # noqa: E402
from common import development as dev  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 264
WF = d244.WF
DROP = ("PER1", "NR1D1")
WF12 = WF.drop(list(DROP))
STEROID = ("MRAP", "CYP11B1", "HSD3B2", "CYP21A2", "MVD", "SCARB1", "MT3", "IDI1", "MC2R", "AMDHD1", "ACSS2", "TMEM97",
           "GSTA4", "STAR", "CYP11A1", "EBP", "FDX1", "CYP17A1", "EPB41L1", "INHA")
ACUTE = ("NR4A1", "NR4A2", "NR4A3")
SLOPE_MIN = 5  # fetuses of a cohort in the readout for a slope of its own


def activity(counts: pd.DataFrame, depth, genes) -> np.ndarray:
    """Equal-weight score: mean over `genes` of the log fraction standardised over the rows given."""
    y = np.log((counts[list(genes)].to_numpy(float) + 0.5) / np.asarray(depth, float)[:, None])
    return np.nanmean((y - y.mean(0)) / y.std(0), axis=1)


def sheet(counts: pd.DataFrame) -> pd.DataFrame:
    """The A244 sheet with the two scores (standardised within each class)."""
    obs = d244.sheet()
    obs["steroid"], obs["acute"] = np.nan, np.nan
    for k, g in obs.groupby("compartment"):
        obs.loc[g.index, "steroid"] = activity(counts.loc[g.index], g.depth, STEROID)
        obs.loc[g.index, "acute"] = activity(counts.loc[g.index], g.depth, ACUTE)
    return obs


def cov(age_df: int = 2, slopes: bool = False, score: str | None = None):
    """Covariates of a set of rows: cohort indicator and common age polynomial (A244), optional cohort-specific linear
    age slopes (cohorts with >= SLOPE_MIN rows) and an activity score column."""
    def f(g):
        Z = [dev.design(g, "age", age_df, ("block",))]
        if slopes:
            a = ((g.age - g.age.mean()) / g.age.std()).to_numpy(float)
            Z += [(a * (g.cohort == c).to_numpy())[:, None] for c, n in g.cohort.value_counts().items() if n >= SLOPE_MIN]
        if score:
            Z.append(g[[score]].to_numpy(float))
        return np.column_stack(Z)
    return f


def tests(obs) -> dict:
    """{name: (rows, covariates, strata, waveform)}."""
    rows = d244.pooled(obs, obs.compartment == "cortex")
    return {"cortex": (rows, cov(), "block", WF),
            "cortex:age": (rows, cov(3, slopes=True), "block", WF),
            "cortex:steroid": (rows, cov(score="steroid"), "block", WF),
            "cortex:-PER1-NR1D1": (rows, cov(), "block", WF12),
            "cortex:all3": (rows, cov(3, slopes=True, score="steroid"), "block", WF12),
            "cortex:acute": (rows, cov(score="acute"), "block", WF)}


REGISTERED = ("cortex:age", "cortex:steroid", "cortex:-PER1-NR1D1")
