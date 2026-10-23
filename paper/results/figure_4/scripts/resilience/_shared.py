"""Loaders and small helpers shared by the scripts of resilience.

Method calls are clocksync's. Everything marked "local: not in clocksync" is a short helper for something the package
has no public call for (listed in the folder README).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import anchors as AN  # noqa: E402
from common import data  # noqa: E402


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
TARGETS = DATASETS / "_atlas" / "a17_targets.csv"
GROUPS = ("control", "case")
RULES = ("learned", "fixed")
STROMAL_VASCULAR = ("fibroblast_stromal", "endothelial", "mural")
IMMUNE = ("lymphocyte", "myeloid")


def targets(kind: str | None = None) -> pd.DataFrame:
    """The A17 registry (`data/_atlas/a17_targets.csv`): name, dataset, kind, within, primary, ..."""
    t = pd.read_csv(TARGETS)
    return t if kind is None else t[t.kind == kind].reset_index(drop=True)


def family_map() -> dict:
    """compartment -> family of the A17 anchor table (one family per compartment name)."""
    t = pd.read_csv(AN.TABLE)
    return t.drop_duplicates("compartment").set_index("compartment").family.to_dict()


def anchor_set(name: str, rule: str):
    """Anchor compartments of an A17 target, or None when the rule gives fewer than 3 ("no anchor")."""
    t = pd.read_csv(AN.TABLE)
    t = t[t.target == name]
    if t.empty:
        return None
    anc = AN.anchors(name, rule)
    if len(anc) < 3 or (rule == "learned" and str(t.status.iloc[0]).startswith("no anchor")):
        return None
    return anc


def subsample_donors(lab: pd.Series, n: int | None, seed: int) -> pd.Series:
    """`n` random labelled donors, both groups kept (at least n // 2 of each when available)."""
    if n is None or n >= len(lab):
        return lab
    rng = np.random.default_rng(seed)
    keep = []
    for g in lab.unique():
        idx = lab.index[lab == g].to_numpy()
        keep += list(rng.permutation(idx)[: max(n // lab.nunique(), 5)])
    return lab.loc[keep]


def load_contrast(row, subsample: int | None = None, seed: int = 0):
    """(pseudobulk of the contrast's donors, donor table, labels, blocks argument for `Dataset.contrast`).

    Labels: the donor-table column named like the target ("control" / "case"). Blocks follow the legacy `within`
    setting: "stratum" -> None (permute within strata), empty -> "none" (global permutation), "block:<col>" -> that
    donor-table column."""
    pb, donors = data.load(row.dataset)
    lab = donors[row.name].dropna()
    lab = lab[lab.isin(GROUPS)]
    lab = subsample_donors(lab, subsample, seed)
    pb = pb[pb.donor.isin(lab.index)].reset_index(drop=True)
    w = row.within
    blocks = None if w == "stratum" else "none" if not isinstance(w, str) else w.split(":", 1)[1]
    return pb, donors.loc[donors.index.isin(lab.index)], lab, blocks


def holm_within(df: pd.DataFrame, p: str, by: str = "target") -> pd.Series:
    """Holm adjustment of column `p` within each value of `by` (NaN where p is NaN)."""
    out = pd.Series(np.nan, index=df.index)
    for _, g in df.groupby(by):
        ok = g[p].notna()
        if ok.any():
            out.loc[g.index[ok]] = cs.stats.holm(g.loc[ok, p].to_numpy())
    return out
