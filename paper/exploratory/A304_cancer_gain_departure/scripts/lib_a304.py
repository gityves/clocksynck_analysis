"""A304: tumour gain and departure on the patient's own normal sample (TCGA paired, 10 cancers).

Registered in paper/docs/registered_experiments.md (A304). Shared by `power_a304.py` and `score_a304.py`.

Model: normal N = h + e_N, tumour T = g h + e_T. Per arm, on centred complex clock-plane states (one projector per cancer,
as A138), Poisson-corrected:
    V_N, V_T   clock-plane variances;   C = sum z_T conj(z_N) / (n - 1)
    G = |C| / V_N                  gain of the tumour on its normal (= |g| rel)
    R = (V_T - |C|^2 / V_N) / V_N^all   departure (unitless; V_N^all = the cancer's normal variance over all pairs)
Pooled over cancers (weights n_c): Q_G = log(sum n G_case / sum n G_ref), Q_R likewise, Q_N for V_N (host check).
Every function takes label matrices (rows = label draws) so observed values, permutation nulls and simulations share
one vectorised path.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "paper" / "exploratory" / "A138_cancer_estimability" / "scripts"))
import q3_lib as q  # noqa: E402  (tcga_states, sh = figure_5 cancer _shared)

sh = q.sh
OUT = HERE.parent / "output"
EXP = REPO / "paper" / "exploratory"
CANCERS = sh.CANCERS
MIN_ARM = 8
CONTRASTS = {"H304.1": "proliferation high vs low", "H304.2": "ER negative vs positive (BRCA)",
             "H304.3": "TP53 altered vs intact", "H304.4": "luminal A vs other PAM50 (BRCA)"}
# GTEx tissue of each cancer's normal (A104 table); BLCA and UCEC: median of the mapped tissues
TISSUE = {"BRCA": "Breast - Mammary Tissue", "COAD": "Colon - Transverse", "KIRC": "Kidney - Cortex",
          "KIRP": "Kidney - Cortex", "LIHC": "Liver", "LUAD": "Lung", "PRAD": "Prostate", "THCA": "Thyroid"}
A104 = REPO / "paper" / "results" / "figure_2" / "output" / "backbone" / "a104_gtex_tissues.csv"


def reliability() -> dict:
    """Cancer -> GTEx coupling of the matching tissue (S / (S + 2 psi), A104); BLCA, UCEC: median of the mapped ones."""
    a = pd.read_csv(A104).set_index("tissue").coupling
    rel = {c: float(a[t]) for c, t in TISSUE.items()}
    med = float(np.median([rel[c] for c in rel]))
    return {c: rel.get(c, med) for c in CANCERS}


# ---- design --------------------------------------------------------------------------------------------------------

def annotations(samples: pd.DataFrame) -> dict:
    """Cancer -> DataFrame (index = patient, pair order) with prolif, purity, er, tp53, subtype."""
    st = pd.read_csv(EXP / "A215_checkpoint_coherence" / "output" / "status_tcga.csv", index_col=0).tp53
    sub = pd.read_csv(EXP / "A217_subtype_matched" / "output" / "subtypes_tcga.csv", keep_default_na=False)
    out = {}
    for c in CANCERS:
        case, _ = sh.pair_frames(samples, c)
        t = pd.read_parquet(sh.group_out(sh.FOLDER) / "cache" / f"tcga_tumours_{c}.parquet").set_index("unit")
        t = t.reindex(case.index)
        assert (t.barcode.to_numpy() == case.barcode.to_numpy()).all(), c
        s = sub[sub.project == c].set_index("patient").subtype.replace("", np.nan)
        out[c] = pd.DataFrame(dict(prolif=t.prolif, purity=t.purity, er=t.er,
                                   tp53=case.barcode.str[:15].map(st).to_numpy(),
                                   subtype=case.index.map(s)), index=case.index)
    return out


def labels(ann: pd.DataFrame, contrast: str, cancer: str) -> np.ndarray:
    """1 = case, 0 = reference, -1 = not in the contrast (pair order)."""
    lab = np.full(len(ann), -1)
    if contrast == "H304.1":
        x = ann.prolif.to_numpy(float)
        lab = (x > np.median(x)).astype(int)
    elif contrast == "H304.2" and cancer == "BRCA":
        lab = np.where(ann.er == "Negative", 1, np.where(ann.er == "Positive", 0, -1))
    elif contrast == "H304.3":
        lab = np.where(ann.tp53 == 1, 1, np.where(ann.tp53 == 0, 0, -1))
    elif contrast == "H304.4" and cancer == "BRCA":
        lab = np.where(ann.subtype == "BRCA_LumA", 1, np.where(ann.subtype.notna(), 0, -1))
    if min((lab == 1).sum(), (lab == 0).sum()) < MIN_ARM:
        return np.full(len(ann), -1)
    return lab


def purity_half(ann: pd.DataFrame) -> np.ndarray:
    x = ann.purity.to_numpy(float)
    return (x > np.median(x)).astype(int)


def cell_from_raw(cancer: str, raw, ann_c: pd.DataFrame, contrast: str) -> dict | None:
    """One cell from (bt, vt, bn, vn) of all pairs of a cancer (pair order): states, variances and labels of the pairs
    in the contrast; V_N^all over all pairs. None if the cancer is not in the contrast."""
    lab = labels(ann_c, contrast, cancer)
    keep = lab >= 0
    if not keep.any():
        return None
    bt, vt, bn, vn = raw
    zt, zn = q.S.complex_state(bt), q.S.complex_state(bn)
    zn_all = zn - zn.mean()
    vn_all = float(np.sum(np.abs(zn_all) ** 2) / (len(zn) - 1) - vn.mean())
    return dict(cancer=cancer, zt=zt[keep], vt=vt[keep], zn=zn[keep], vn=vn[keep], lab=lab[keep],
                n=int(keep.sum()), vn_all=vn_all, ann=ann_c[keep])


def design(states: dict, ann: dict, contrast: str) -> list[dict]:
    """Cells of one contrast: per cancer the complex states, Poisson variances and labels of the pairs in the contrast."""
    cells = [cell_from_raw(c, states[c]["raw"], ann[c], contrast) for c in CANCERS if c in states]
    return [c for c in cells if c is not None]


# ---- statistics (label matrices: P x n, values 0 / 1) ----------------------------------------------------------------

def arm_moments(zt, vt, zn, vn, L):
    """V_N, V_T, C of the arm L (P x n masks), centred within arm, Poisson-corrected. Returns three (P,) arrays."""
    L = L.astype(float)
    m = L.sum(1)
    st, sn = L @ zt, L @ zn
    vN = (L @ np.abs(zn) ** 2 - np.abs(sn) ** 2 / m) / (m - 1) - (L @ vn) / m
    vT = (L @ np.abs(zt) ** 2 - np.abs(st) ** 2 / m) / (m - 1) - (L @ vt) / m
    C = (L @ (zt * np.conj(zn)) - st * np.conj(sn) / m) / (m - 1)
    return vN, vT, C


def cell_values(cell: dict, labs: np.ndarray, zt=None, zn=None) -> dict:
    """Per arm (case 1, ref 0) G, R, V_N, C for each label row. zt / zn override the cell's states (simulation)."""
    zt = cell["zt"] if zt is None else zt
    zn = cell["zn"] if zn is None else zn
    out = {}
    for a in (1, 0):
        vN, vT, C = arm_moments(zt, cell["vt"], zn, cell["vn"], labs == a)
        with np.errstate(divide="ignore", invalid="ignore"):
            G = np.where(vN > 0, np.abs(C) / vN, np.nan)
            R = np.where(vN > 0, (vT - np.abs(C) ** 2 / vN) / cell["vn_all"], np.nan)
        out[a] = dict(G=G, R=R, VN=vN / cell["vn_all"], C=C)
    return out


def pooled(cells: list, vals: list) -> dict:
    """Q_G, Q_R, Q_N (P,) pooled over cells with weights n; a cell with a non-finite value leaves that row."""
    res = {}
    for key, name in (("G", "Q_G"), ("R", "Q_R"), ("VN", "Q_N")):
        num = np.zeros_like(vals[0][1][key]); den = np.zeros_like(num)
        for c, v in zip(cells, vals):
            a, b = v[1][key], v[0][key]
            ok = np.isfinite(a) & np.isfinite(b)
            num += np.where(ok, c["n"] * a, 0.0); den += np.where(ok, c["n"] * b, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            res[name] = np.where((num > 0) & (den > 0), np.log(num / den), np.nan)
    return res


def perm_labels(lab: np.ndarray, n_perm: int, rng, strata=None) -> np.ndarray:
    """n_perm permutations of lab (within strata when given)."""
    P = np.empty((n_perm, len(lab)), int)
    if strata is None:
        strata = np.zeros(len(lab), int)
    for s in np.unique(strata):
        i = np.flatnonzero(strata == s)
        P[:, i] = lab[i][np.argsort(rng.random((n_perm, len(i))), axis=1)]
    return P


def p_less(null, obs):
    x = np.asarray(null)[np.isfinite(null)]
    return float((1 + np.sum(x <= obs)) / (1 + len(x)))


def p_greater(null, obs):
    x = np.asarray(null)[np.isfinite(null)]
    return float((1 + np.sum(x >= obs)) / (1 + len(x)))
