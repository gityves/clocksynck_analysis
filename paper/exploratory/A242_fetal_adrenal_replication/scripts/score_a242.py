"""A242 scoring (registered in `paper/docs/registered_experiments.md`): does the clock of the human fetal adrenal
(A189, Cao 2020, 7 fetuses) replicate in open data from other laboratories?

One whole-adrenal pseudobulk per fetus and cohort (`prep_adrenal.py`: all cells of the fetus summed; bulk libraries as
released). L1 signed coherence C across the fetuses (`common.development`), pooled over cohorts: log fractions are
residualised gene by gene on a cohort indicator and a common polynomial of the age stated by the source, and the
per-gene permutations of the null are done within cohort (cohorts differ in platform, depth and noise, so their
fetuses are not exchangeable). Companions: expression-matched and responsive pseudo-clock sets among the genes present
in every cohort of the test, adjusted and scored in the same way.

Tests
- `late`:        fetuses from 12 weeks (the age range of A189), cohorts with >= 2 such fetuses; linear age trend.
- `late:deep`:   the same without the cohorts whose median clock-gene count is < 100 (han, roost).
- `late:-<c>`:   `late` without cohort c (leave one cohort out).
- `early`:       fetuses before 12 weeks, cohorts with >= 2; quadratic age trend.
- `all`:         every fetus of the cohorts with >= 2; quadratic age trend.
- `kameneva`, `wang`: one cohort, every age; quadratic age trend.
- `schroder`:    13 bulk fetal adrenals without recorded age (Schroder 2022, targeted transcriptome): sequencing batch
                 and the first principal component of the non-clock genes (`development.surrogates`), permutations
                 within batch; `schroder:pc0` and `schroder:pc2` with no and two components.

    python score_a242.py [--jobs 12] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design.csv, l1.csv, verdicts.csv (`_check` with --check: 200 permutations, 10 sets).
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
SRC = DATA / "natexp" / "fetal_adrenal_open"
SEED = 242
WF = cs.gtex()
POOLED = ["encode_pa", "kameneva", "dong", "wang", "han", "sasaki", "roost", "kildisiute"]
SHALLOW = ["han", "roost"]  # median clock-gene count < 100 (`clock_detection.csv`)
LATE = 12.0  # weeks, as stated by each source
MIN_FETUSES = 2  # per cohort and test: a single fetus has no residual
FAMILIES = {"H242.1": ["late", "schroder"], "H242.2": ["early"]}  # Holm within a hypothesis


def sheet() -> pd.DataFrame:
    """One row per fetus and cohort: donor, compartment, cohort, block (cohort; sequencing batch for Schroder), age
    (weeks, as stated by the source; none for Schroder) and depth."""
    s = pd.read_csv(SRC / "samples.csv").set_index("id")
    s = s[s.cohort.isin([*POOLED, "schroder"])]
    block = s.cohort.where(s.cohort != "schroder", "schroder:" + s.convention.str.replace("batch ", ""))
    return pd.DataFrame(dict(donor=s.index, compartment="adrenal", cohort=s.cohort, block=block, age=s.age,
                             depth=s.depth, stratum="all"))


def blocks(obs) -> np.ndarray:
    """Indicators of the blocks of a set of rows."""
    return dev.design(obs, None, 0, ("block",))


def load():
    """Counts (fetuses x gene symbols; NaN = symbol absent from the source) and the sheet, with the surrogate
    components of the Schroder libraries (`pc1`, `pc2`: non-clock genes, batch regressed out)."""
    obs = sheet()
    counts = pd.read_parquet(SRC / "pseudobulk.parquet").loc[obs.index]
    s = obs.cohort == "schroder"
    pcs = dev.surrogates(counts[s].dropna(axis=1), obs.depth[s], blocks(obs[s]), 2, WF.genes)
    obs["pc1"], obs["pc2"] = np.nan, np.nan
    obs.loc[s, ["pc1", "pc2"]] = pcs
    return counts, obs


def cov(df: int, pcs: int = 0):
    """Design of a set of rows: block indicators, a common polynomial of age (df) and `pcs` surrogate components."""
    def f(g):
        Z = dev.design(g, "age" if df else None, df, ("block",))
        return np.column_stack([Z, g[[f"pc{i + 1}" for i in range(pcs)]].to_numpy(float)]) if pcs else Z
    return f


def pooled(obs, rows, without=()):
    """`rows` restricted to the pooled cohorts (minus `without`) that have >= MIN_FETUSES fetuses among them."""
    rows = rows & obs.cohort.isin([c for c in POOLED if c not in without])
    n = obs[rows].groupby("cohort").donor.transform("size").reindex(obs.index, fill_value=0)
    return rows & (n >= MIN_FETUSES)


def tests(obs):
    """L1 tests {name: (rows, covariates, permutation strata)}; no L2 (one compartment)."""
    late, every = obs.age >= LATE, obs.age.notna()
    l1 = {"late": (pooled(obs, late), cov(1), "block"),
          "late:deep": (pooled(obs, late, SHALLOW), cov(1), "block"),
          "early": (pooled(obs, obs.age < LATE), cov(2), "block"),
          "all": (pooled(obs, every), cov(2), "block")}
    for c in obs.cohort[l1["late"][0]].unique():
        l1[f"late:-{c}"] = (pooled(obs, late, [c]), cov(1), "block")
    for c in ("kameneva", "wang"):
        l1[c] = (every & (obs.cohort == c), cov(2), "block")
    s = obs.cohort == "schroder"
    l1.update({"schroder": (s, cov(0, 1), "block"), "schroder:pc0": (s, cov(0), "block"),
               "schroder:pc2": (s, cov(0, 2), "block")})
    return l1, {}


if __name__ == "__main__":
    dev.main(__doc__, load, tests, WF, SEED, OUT, FAMILIES)
