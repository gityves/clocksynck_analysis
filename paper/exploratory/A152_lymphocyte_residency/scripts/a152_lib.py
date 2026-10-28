"""Shared helpers of A152 (lymphocyte residency): A6 per-cell residency rule, lymphocyte cell types of the two Census
atlases, pseudobulk assembly on the A17 grouping, the A17 fit and the within-donor label permutation.

Registered: paper/docs/registered_experiments.md (A152). Legacy per-cell rule: legacy/desync/common/cell_states.py (A6).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "atlas_fits"))
sys.path.insert(0, str(PAPER / "common" / "datasets" / "_atlas"))
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = PAPER / "exploratory" / "A152_lymphocyte_residency" / "output"
SEED = 152
RESIDENT, CIRCULATING = ["CD69", "ITGAE", "CXCR6"], ["SELL", "CCR7", "KLF2", "S1PR1"]
MARKERS = RESIDENT + CIRCULATING
GENES = list(cs.gtex().genes)
STATES = np.array(["circulating_like", "intermediate", "resident_like"])
# tissue -> A17 target, dataset folder, Census dataset id, A17 classes replaced by the two states, class of the CD8 split
CFG = {
    "heart": dict(target="heart", atlas="heart_global", dataset="364bd0c7-f7fd-48ed-99c1-ae26872b1042",
                  cells="a152_heart_global", replaced=["heart:lymphocyte"], cd8=None, prefix="heart"),
    "breast": dict(target="organ_breast", atlas="breast_kumar", dataset="842c6f5d-4a94-4eef-8510-8c792d1124bc",
                   cells="a152_breast_kumar", replaced=[f"breast:{c}" for c in ("CD4_T", "CD8_T", "NK", "T_other", "Treg")],
                   cd8="breast:CD8_T", prefix="breast"),
}
NOT_LYMPH = r"mast cell"  # matches the A17 regexes ("t cell" inside "mast cell"); not a lymphocyte


def classes_of(tissue: str, cell_type: pd.Series) -> pd.Series:
    """A17 class of a Census cell type: `<organ>:<class>` by the builder's own ordered regexes."""
    import census as C
    import heart as H
    import organ as O
    cl = H.CLASSES if tissue == "heart" else O.COARSE
    return CFG[tissue]["prefix"] + ":" + C.classify(cell_type.astype(str), cl)


def lymph_cell_types(tissue: str) -> list:
    """Census cell types of the replaced classes (incl. mast cells that the A17 regex puts there), from the atlas obs."""
    import census as C
    o, _ = C.read(C.ATLAS / CFG[tissue]["atlas"], "clock")
    cls = classes_of(tissue, o.cell_type)
    return sorted(set(o.cell_type.astype(str)[cls.isin(CFG[tissue]["replaced"]).to_numpy()]))


# ---- A6 residency state -------------------------------------------------------------------------------------------
def tertile(score: np.ndarray, block: np.ndarray, rng) -> np.ndarray:
    """Tertile state within each block; ties broken at random (A6 `tertile`, unchanged)."""
    out = np.empty(len(score), dtype=object)
    for b in np.unique(block):
        i = np.flatnonzero(block == b)
        order = np.lexsort((rng.random(len(i)), score[i]))
        rank = np.empty(len(i))
        rank[order] = np.arange(len(i))
        out[i] = STATES[np.minimum((3 * rank / len(i)).astype(int), 2)]
    return out


def score_cells(cells: pd.DataFrame) -> pd.Series:
    """mean x of RESIDENT minus mean x of CIRCULATING, x = log1p(1e4 count / raw_sum)."""
    x = np.log1p(1e4 * cells[MARKERS].to_numpy(float) / cells.raw_sum.to_numpy(float)[:, None])
    x = pd.DataFrame(x, columns=MARKERS, index=cells.index)
    return x[RESIDENT].mean(1) - x[CIRCULATING].mean(1)


# ---- pseudobulk assembly ------------------------------------------------------------------------------------------
def H_frame(tissue: str) -> pd.DataFrame:
    import heart as H
    import organs as M
    import organ as O
    return O.groups(H.spec_global() if tissue == "heart" else M.spec("breast"))


def pseudobulk(cells: pd.DataFrame, label: pd.Series, tissue: str, min_cells: int = 20) -> pd.DataFrame:
    """Donor x label pseudobulks (standard-table columns) of the cells with a label, in the donor's A17 stratum."""
    pre = CFG[tissue]["atlas"]
    key = ["donor", "compartment"]
    c = cells.assign(donor=pre + ":" + cells.donor_id.astype(str), compartment=label.to_numpy())
    c = c[c.compartment.notna()]
    num = GENES + ["raw_sum"]
    g = c.groupby(key)[num].sum()
    g["n_cells"] = c.groupby(key).size()
    g = g.rename(columns={"raw_sum": "depth"}).reset_index()
    g["stratum"] = g.donor.map(donor_batch_cached(tissue))
    g = g[g.n_cells >= min_cells]
    return g[["donor", "compartment", "stratum", "n_cells", "depth"] + GENES].reset_index(drop=True)


_CACHE = {}


def donor_batch_cached(tissue: str) -> pd.Series:
    if tissue not in _CACHE:
        _CACHE[tissue] = H_frame(tissue).drop_duplicates("resection").set_index("resection").batch
    return _CACHE[tissue]


# ---- A17 fit ------------------------------------------------------------------------------------------------------
def anchors(target: str) -> list:
    from compute_atlas_fits import anchor_sets
    return anchor_sets(target)["learned"]


def fit_table(pb: pd.DataFrame, donors: pd.DataFrame, anc: list, name: str, targets=None, n_boot: int = 0) -> pd.DataFrame:
    ds = cs.Dataset(pb, anchors=anc, donors=donors, name=name)
    return ds.fit(n_boot=n_boot, gate=False, targets=targets, seed=SEED).table


def permute_pairs(pb: pd.DataFrame, a: str, b: str, rng) -> pd.DataFrame:
    """Swap the labels `a` / `b` within each donor that has both, each donor with probability 1/2."""
    both = pb[pb.compartment.isin([a, b])].groupby("donor").compartment.nunique()
    sw = both.index[(both == 2).to_numpy()]
    flip = set(sw[rng.random(len(sw)) < 0.5])
    out = pb.copy()
    m = out.donor.isin(flip) & out.compartment.isin([a, b])
    out.loc[m, "compartment"] = out.loc[m, "compartment"].map({a: b, b: a})
    return out


def compartments_of(target: str) -> list:
    from _shared import compartments
    return compartments(target)
