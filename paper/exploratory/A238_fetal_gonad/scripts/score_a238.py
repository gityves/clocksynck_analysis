"""A238 scoring (registered in `paper/docs/registered_experiments.md`): a clock in the cell classes of the human fetal
gonad, 6-21 weeks post fertilisation (Garcia-Alonso 2022, CELLxGENE collection 661a402a, 10x 5').

L1 per cell class across donors and L2 across the cell classes of one gonad, on log fractions of the donor x class
pseudobulks (>= 200 cells; `prep_fetal_gonad/scripts/prep_gonad.py`) residualised gene by gene within each class on
age (quadratic, one trend per sex), sex and assay, with expression-matched and responsive pseudo-clock companions
(`common.development`). L2 strata: sex x (before | from 13 weeks).

    python score_a238.py [--jobs 12] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design.csv, l1.csv, l2.csv, l2_pseudo.csv, verdicts.csv (`_check` with --check).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SRC = DATA / "natexp" / "garcia_alonso2022"
SEED = 238
WF = cs.gtex()
MIN_CELLS, MIN_N, SPLIT = 200, 10, 13  # cells per pseudobulk; donors per test; weeks post fertilisation
CLASSES = ["mesenchymal_stromal", "supporting", "endothelial", "epithelial"]
SETS = {"classes4": CLASSES, "classes3": CLASSES[:3], "pair": CLASSES[:2]}
FAMILIES = {"H238.1": [f"{c}:all" for c in CLASSES], "H238.2": ["classes4:all"]}  # Holm within a hypothesis
SUBSETS = {"all": lambda o: o.age > 0, "female": lambda o: o.sex == "female", "male": lambda o: o.sex == "male",
           "13+wpf": lambda o: o.age >= SPLIT}


def load():
    """Counts (pseudobulks x gene symbols) and the sheet (donor, compartment, depth, age, sex, assay, stratum)."""
    a = ad.concat([ad.read_h5ad(SRC / f"garcia_{s}_pseudobulk_class.h5ad") for s in ("female", "male")], merge="same")
    a = a[(a.obs.gonad & (a.obs.n_cells >= MIN_CELLS) & a.obs.compartment.isin(CLASSES)).to_numpy()]
    x = np.asarray(a.X.toarray() if hasattr(a.X, "toarray") else a.X, float)
    top = pd.Series(x.mean(0), index=np.arange(a.n_vars)).groupby(a.var.symbol.astype(str).map(canonical).to_numpy()).idxmax()
    counts = pd.DataFrame(x[:, top.to_numpy()], index=a.obs_names, columns=top.index)
    o = a.obs
    obs = pd.DataFrame(dict(donor=o.donor.astype(str), compartment=o.compartment.astype(str), depth=o.depth.astype(float),
                            age=o.wpf.astype(float), sex=o.sex.astype(str), assay=o.assay.astype(str)))
    obs["stratum"] = obs.sex + np.where(obs.age >= SPLIT, "|late", "|early")
    return counts, obs


def cov(g):
    """Design of a set of pseudobulks: quadratic age trend of each sex, sex and assay (constant columns dropped)."""
    age = [dev.basis(g.age, 2) * (g.sex == s).to_numpy()[:, None] for s in ("female", "male")]
    Z = np.column_stack([*age, dev.design(g, None, 0, ("sex", "assay"))])
    return Z[:, np.ptp(Z, 0) > 0]


def tests(obs):
    """L1 and L2 tests {name: (rows, covariates)}, named `<class or class set>:<subset>`; >= MIN_N donors."""
    l1, l2 = {}, {}
    for sub, rule in SUBSETS.items():
        s = rule(obs)
        for c in CLASSES:
            rows = s & (obs.compartment == c)
            if rows.sum() >= MIN_N:
                l1[f"{c}:{sub}"] = (rows, cov)
        for tag, use in SETS.items():
            rows = s & obs.compartment.isin(use)
            rows &= obs[rows].groupby("donor").compartment.transform("nunique").reindex(obs.index, fill_value=0) >= 2
            if obs[rows].donor.nunique() >= MIN_N:
                l2[f"{tag}:{sub}"] = (rows, cov)
    return l1, l2


if __name__ == "__main__":
    dev.main(__doc__, load, tests, WF, SEED, OUT, FAMILIES)
