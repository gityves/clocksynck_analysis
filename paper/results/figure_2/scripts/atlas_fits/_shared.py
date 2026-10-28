"""Helpers shared by atlas_fits, host_identity and backbone: the A17 target registry,
the table of every healthy fit and contrast of the re-analysis (donors, compartments, labels, permutation blocks) and the
host-identity diagnostics that clocksync has no public call for ("local: not in clocksync").

A17 target = (dataset, kind, reference, label column, blocks). The standard table of a dataset is the union of the rows
of every target that reads it, so a target is selected here: healthy fits keep the donors of the healthy fit
(`healthy` column) and, for every target, the compartments listed in the frozen anchor table (`common/anchors_a17.csv`,
which holds the compartments of the legacy fit); contrasts keep the donors with a case / control label.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from common import data  # noqa: E402
from common.anchors import TABLE as ANCHOR_TABLE  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
REGISTRY = DATASETS / "_atlas" / "a17_targets.csv"
FOLDER_5 = "atlas_fits"
FOLDER_4 = "backbone"
FOLDER_S04 = "host_identity"


def registry() -> pd.DataFrame:
    return pd.read_csv(REGISTRY).set_index("name", drop=False)


def frozen() -> pd.DataFrame:
    return pd.read_csv(ANCHOR_TABLE)


def compartments(name: str) -> list:
    a = frozen()
    return a.loc[a.target == name, "compartment"].tolist()


class Target:
    """One A17 target: standard table restricted to the fit's donors and compartments."""

    def __init__(self, name: str, subsample: int | None = None, seed: int = 0, new_classes: bool = False):
        r = registry().loc[name]
        self.name, self.kind, self.row = name, r.kind, r
        self.study, self.tissue, self.ref = r.study, r.tissue, r.ref
        pb, donors = data.load(r.dataset)
        # new_classes: the anchor learning (compute_coherence) also reads the compartments that the class-rule fix added
        # to the table (common/class_fix_new_compartments.csv); elsewhere the anchor table defines the compartments
        comps = compartments(name)
        if new_classes:
            add = pd.read_csv(ANCHOR_TABLE.with_name("class_fix_new_compartments.csv"))
            comps = comps + [c for c in add.loc[add.dataset == r.dataset, "compartment"] if c not in comps]
        if r.kind == "healthy":
            keep = donors.index[donors["healthy"].astype(bool)] if "healthy" in donors else donors.index
            lab = None
        else:
            lab = donors[name].dropna().astype(str)
            keep = lab.index
        pb = pb[pb.donor.isin(keep) & pb.compartment.isin(comps)].reset_index(drop=True)
        if subsample is not None:
            sel = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
            pb = pb[pb.donor.isin(sel)].reset_index(drop=True)
        self.pb = pb
        self.donors = donors.loc[donors.index.isin(pb.donor.unique())]
        self.labels = None if lab is None else lab.reindex(self.donors.index).dropna()
        self.blocks = None
        if r.kind == "contrast":
            w = r.within
            if isinstance(w, str) and w.startswith("block:"):
                self.blocks = self.donors[w.split(":", 1)[1]].astype(str)
            elif not isinstance(w, str):
                self.blocks = "none"  # one global permutation (labels cross studies)

    def dataset(self, anchors=None, **kw) -> cs.Dataset:
        return cs.Dataset(self.pb, anchors=anchors, donors=self.donors, name=self.name, **kw)


# ---- identity diagnostics (local: not in clocksync) ---------------------------------------------------------------

def similarity(fa: pd.DataFrame, fb: pd.DataFrame, strata: pd.Series | None = None) -> float:
    """|complex correlation| of two host states (x, y columns) over their common donors, centred within strata
    (the legacy `anchored.similarity`; `clocksync.stats.ccorr` takes the pairs of states, this aligns donors)."""
    common = fa.index.intersection(fb.index)
    a, b = np.array(fa.loc[common, ["x", "y"]], float), np.array(fb.loc[common, ["x", "y"]], float)
    if strata is not None:
        s = strata.reindex(common).astype(str).to_numpy()
        for g in np.unique(s):
            m = s == g
            a[m], b[m] = a[m] - a[m].mean(0), b[m] - b[m].mean(0)
    else:
        a, b = a - a.mean(0), b - b.mean(0)
    za, zb = a[:, 0] + 1j * a[:, 1], b[:, 0] + 1j * b[:, 1]
    den = np.sqrt(np.sum(np.abs(za) ** 2) * np.sum(np.abs(zb) ** 2))
    return float(np.abs(np.sum(za * np.conj(zb))) / den) if den > 0 else np.nan


def host_of(ds: cs.Dataset, comp: str | None = None) -> pd.DataFrame:
    """Full host state (x, y per donor) of the dataset's anchors (one cheap fit: a single target compartment)."""
    comp = comp or ds.compartments[0]
    return ds.fit(n_boot=0, gate=False, targets=[comp]).host


def donor_strata(ds: cs.Dataset) -> pd.Series:
    t = ds.table
    return t.groupby("resection").stratum.first().astype(str) if "stratum" in t else None


def agreement_free(ds: cs.Dataset, host: pd.DataFrame | None = None) -> float:
    """Similarity of the free host (every compartment an anchor) and the anchored host over their common donors."""
    free = cs.Dataset(ds.table.rename(columns={"resection": "donor"}), anchors=None, donors=ds.donors, name=ds.name)
    return similarity(host_of(free), host if host is not None else host_of(ds), donor_strata(ds))


def stability(ds: cs.Dataset, n_boot: int = 200, seed: int = 0, threshold: float = 0.7) -> dict:
    """Share of donor-bootstrap refits of the anchored host whose state correlates > threshold with the full-data host
    (first copy of each distinct donor drawn); legacy `anchored.stability` (without its warm start)."""
    base = host_of(ds)
    comp = ds.compartments[0]
    strata = donor_strata(ds)
    rng = np.random.default_rng(seed)
    sims = []
    for _ in range(n_boot):
        try:
            hb = ds.resample(rng).fit(n_boot=0, gate=False, targets=[comp]).host
        except (ValueError, np.linalg.LinAlgError):
            continue
        orig = pd.Series(hb.index.astype(str).str.split("#").str[0], index=hb.index)
        first = hb[~orig.duplicated().to_numpy()]
        first.index = orig[first.index]
        sims.append(similarity(first, base, strata))
    sims = np.array(sims)
    return dict(stability=float(np.mean(sims > threshold)), median_similarity=float(np.median(sims)), n=len(sims))
