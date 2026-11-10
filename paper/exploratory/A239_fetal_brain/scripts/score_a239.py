"""A239 scoring (registered in `paper/docs/registered_experiments.md`): a clock shared by the regions of the human
brain before and after birth (Kang 2011, GSE25219, Affymetrix exon arrays, 57 brains from 5.7 weeks post conception
to 82 years).

L2 across the regions of one brain and L1 on the mean over the regions of a brain, on log intensities residualised
gene by gene within each region on age (polynomial of log age) and sex, with expression-matched and responsive
pseudo-clock companions (`common.development`). Arrays: the samples of a brain x region are averaged in log2, the
intensity 2^x stands for a count at a constant depth, and the gene-level table lacks PER3 and NR1D1 (12 waveform
genes). L2 strata: stages 1-3 | 4-6 | 7 | 8-9 | 10-11 | 12-15. The postnatal brains are the positive control.

    python score_a239.py [--jobs 12] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design.csv, l1.csv, l2.csv, l2_pseudo.csv, verdicts.csv (`_check` with --check).
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
SRC = DATA / "natexp" / "kang2011"
SEED = 239
WF = cs.gtex().drop(["PER3", "NR1D1"])
REGIONS = ["A1C", "AMY", "CBC", "DFC", "HIP", "IPC", "ITC", "M1C", "MD", "MFC", "OFC", "S1C", "STC", "STR", "V1C",
           "VFC"]  # the regions sampled from stage 3 onwards
GROUPS = {"prenatal": (1, 7), "pcw13-22": (4, 6), "infant+child": (8, 11), "adult": (12, 15), "postnatal": (8, 15)}
BINS = [3.5, 6.5, 7.5, 9.5, 11.5]  # strata: stages 1-3 | 4-6 | 7 | 8-9 | 10-11 | 12-15
FAMILIES = {"control L2": ["regions:postnatal"], "control L1": ["brain_mean:postnatal"],
            "H239.1": ["regions:prenatal"], "H239.2": ["brain_mean:prenatal"]}
MIN_REGIONS, MIN_N, DEPTH = 6, 8, 1e6  # regions per brain; brains per test; nominal depth of an array


def load():
    """Intensities 2^x (brain x region, and one `brain_mean` row per brain with >= MIN_REGIONS regions) and the sheet
    (donor, compartment, depth, age = log age, sex, stage, rin, stratum)."""
    s = pd.read_csv(SRC / "samples.csv")
    s = s[s.region.isin(REGIONS)]
    x = pd.read_parquet(SRC / "expr_gene.parquet")[s.gsm].T
    key = (s.brain + "|" + s.region).to_numpy()
    x.columns = x.columns.map(canonical)
    m = x.groupby(key).mean()
    o = s.assign(key=key).groupby("key").agg(donor=("brain", "first"), compartment=("region", "first"),
                                              age=("age_log", "first"), sex=("sex", "first"),
                                              stage=("stage", "first"), rin=("rin", "mean")).loc[m.index]
    full = o.groupby("donor").compartment.transform("nunique") >= MIN_REGIONS
    m, o = m[full.to_numpy()], o[full]
    mean = m.groupby(o.donor.to_numpy()).mean()
    om = o.groupby("donor").agg(age=("age", "first"), sex=("sex", "first"), stage=("stage", "first"),
                                rin=("rin", "mean")).loc[mean.index]
    om = om.assign(donor=om.index, compartment="brain_mean")
    mean.index = om.index = om.index + "|mean"
    obs = pd.concat([o, om[o.columns]])
    obs["depth"] = DEPTH
    obs["stratum"] = "s" + pd.Series(np.digitize(obs.stage, BINS), index=obs.index).astype(str)
    return 2.0 ** pd.concat([m, mean]), obs


def cov(df, rin=False):
    """Design of a set of rows: polynomial of log age (df), sex, and the RNA integrity number when asked."""
    def f(g):
        Z = dev.design(g, "age", df, ("sex",))
        return np.column_stack([Z, g.rin - g.rin.mean()]) if rin else Z
    return f


def tests(obs):
    """L1 (`brain_mean:<group>`) and L2 (`regions:<group>`) tests {name: (rows, covariates)}; >= MIN_N brains."""
    l1, l2 = {}, {}
    mean = obs.compartment == "brain_mean"
    for tag, (lo, hi) in GROUPS.items():
        s = obs.stage.between(lo, hi)
        if obs[s].donor.nunique() < MIN_N:
            continue
        df = 2 if tag in ("prenatal", "postnatal") else 1
        l1[f"brain_mean:{tag}"], l2[f"regions:{tag}"] = (s & mean, cov(df)), (s & ~mean, cov(df))
        if tag in ("prenatal", "postnatal"):
            l1[f"brain_mean:{tag}:rin"], l2[f"regions:{tag}:rin"] = (s & mean, cov(df, True)), (s & ~mean, cov(df, True))
    return l1, l2


if __name__ == "__main__":
    dev.main(__doc__, load, tests, WF, SEED, OUT, FAMILIES)
