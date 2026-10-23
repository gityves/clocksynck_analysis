"""Helpers shared by the scripts of gtex_cues (GTEx: covariates, cue sets, internal vs external
time). Clock quantities come from `clocksync`; everything here is data access or a small local statistic
("local: not in clocksync", listed in the README).

GTEx time of death and protected phenotypes are controlled access: callers hold them in memory and write aggregates.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import gtex  # noqa: E402

W = 2 * np.pi / 24

# v2_gtex.TISSUES (legacy `scripts/v2_gtex.py`): the 15 tissues whose clock genes form the host of a target tissue
TISSUES = (
    "Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
    "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
    "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum", "Brain - Cortex", "Whole Blood",
)
REF = "Muscle - Skeletal"
BODY = tuple(t for t in TISSUES if t != "Pituitary" and not t.startswith("Brain"))  # A81 / A92 body tissues
FAST, VENT = (1, 2), (0,)
GC = ("FKBP5", "TSC22D3", "ZBTB16", "KLF9", "KLF15", "DDIT4")  # A27 glucocorticoid score genes
CAUSE_FILE_NOTE = "cause-of-death groups are built from the protected DTHCOD field (aggregates only)"


# ---- genome-wide expression of one tissue (local: not in clocksync; gtex.load has no RIN / centre / all-gene CPM) ------

def expression(tissue: str, hardy=(1, 2), tod: bool = True, median_cpm: float | None = 5.0):
    """(Y, genes, obs) of one GTEx tissue: log(CPM + 0.5) of every gene (library-size CPM), the genes with median CPM
    >= `median_cpm` (None: all) and the sample annotations (subject, hardy, tod mod 24, sex, age bracket and midpoint,
    ischaemic time in minutes, rin, centre). One sample per donor (first). Only the rows of this tissue are read."""
    import anndata as ad

    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        keep = (o.DTHHRDY.isin(list(hardy)) & (o.tissue.astype(str) == tissue)).to_numpy()
        if tod:
            keep = keep & o.tod.notna().to_numpy()
        idx = np.flatnonzero(keep)
        X = a.X[idx]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        genes = np.asarray(a.var_names)
        ob = o.iloc[idx].copy()
    finally:
        a.file.close()
    cpm = X / ob.library_size.to_numpy(float)[:, None] * 1e6
    first = ~ob.subject.astype(str).duplicated().to_numpy()
    cpm, ob = cpm[first], ob[first]
    ok = np.ones(cpm.shape[1], bool) if median_cpm is None else np.median(cpm, 0) >= median_cpm
    obs = pd.DataFrame({"subject": ob.subject.astype(str).to_numpy(), "hardy": ob.DTHHRDY.astype(int).to_numpy(),
                        "tod": ob.tod.to_numpy(float) % 24, "sex": ob.sex.astype(str).to_numpy(),
                        "age_bracket": ob.age_bracket.astype(str).to_numpy(),
                        "age": ob.age_bracket.astype(str).map(gtex.AGE).astype(float).to_numpy(),
                        "ischemic_min": ob.ischemic_time.astype(float).to_numpy(),
                        "rin": ob.rin.astype(float).to_numpy(), "center": ob.center.astype(str).to_numpy()})
    return np.log(cpm[:, ok] + 0.5), genes[ok], obs


def d2_covariates(ob: pd.DataFrame) -> np.ndarray:
    """The D2 covariates of the gene-level GTEx tests (legacy `diag2_degrade.covariates`): sex, age bracket and Hardy
    class as dummies (first level dropped), ischaemic time (z-scored, missing -> median) and its square."""
    Z = pd.get_dummies(pd.DataFrame(dict(sex=ob.sex.astype(str), age=ob.age_bracket.astype(str),
                                         hardy=ob.hardy.astype(str))), drop_first=True).to_numpy(float)
    isch = ob.ischemic_min.astype(float).fillna(ob.ischemic_min.astype(float).median()).to_numpy()
    isch = (isch - isch.mean()) / (isch.std() + 1e-9)
    return np.c_[Z, isch, isch ** 2]


def lstsq_fit(Y: np.ndarray, cols: np.ndarray, Z: np.ndarray):
    """Least squares of every gene (columns of Y) on [1, Z, cols]; returns (B, design)."""
    D = np.c_[np.ones(len(Y)), Z, cols]
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    return B, D


def partial_f(Y: np.ndarray, D: np.ndarray, idx) -> tuple[np.ndarray, np.ndarray]:
    """Coefficients and partial-F p of the columns `idx` of D for every gene (columns of Y)."""
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    rss1 = ((Y - D @ B) ** 2).sum(0)
    D0 = np.delete(D, idx, axis=1)
    B0, *_ = np.linalg.lstsq(D0, Y, rcond=None)
    rss0 = ((Y - D0 @ B0) ** 2).sum(0)
    dfr = len(Y) - D.shape[1]
    return B, stats.f.sf(((rss0 - rss1) / len(idx)) / (rss1 / dfr), len(idx), dfr)


def timed_cosinor(Y: np.ndarray, t_h: np.ndarray, Z: np.ndarray):
    """Per-gene cosinor on the time of death with covariates: complex loading c (cos + i sin), partial-F p."""
    D = np.c_[np.ones(len(Y)), Z, np.cos(W * t_h), np.sin(W * t_h)]
    B, p = partial_f(Y, D, [D.shape[1] - 2, D.shape[1] - 1])
    return B[-2] + 1j * B[-1], p


# ---- the host (internal time) of a donor from the clock genes of other tissues --------------------------------------

def clock_table(hardy=(1, 2), tissues=TISSUES, require_tod: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """GTEx clock-gene pseudobulk (donor x tissue, depth = library size) of the given Hardy classes and tissues
    (through `common.gtex.load`) and the donor table (hardy, tod, sex, age, ischaemic_h; controlled access)."""
    return gtex.load(tissues=list(tissues), hardy=hardy, genes="clock", depth="library_size", require_tod=require_tod)


def host_states(pb: pd.DataFrame, tissues, donors=None, name: str = "") -> pd.DataFrame:
    """Host clock state (x, y; phase_h; radius) of every donor from the clock genes of `tissues` (all anchors, as the
    legacy `host_f`): `Dataset.fit(...).host`."""
    tissues = list(tissues)
    sub = pb[pb.compartment.isin(tissues)]
    if donors is not None:
        sub = sub[sub.donor.isin(set(donors))]
    ds = cs.Dataset(sub.reset_index(drop=True), anchors=tissues, name=name)
    return ds.fit(n_boot=0, gate=False, targets=[tissues[0]]).host
