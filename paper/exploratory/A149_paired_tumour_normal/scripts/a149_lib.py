"""A149 library: cohort tables, paired states and permutation nulls of H149.1-H149.3.

Registered in paper/docs/registered_experiments.md (A149, notes 1-2). Statistics as A138 Q3 (`q3_lib`): states of two
compartments on one projector, centred within cohort x chemistry, Poisson-corrected |complex correlation|; nulls permute
the rows of the second compartment across patients within cohort.

S1 (H149.1)  malignant (or tumour epithelium) vs adjacent-normal epithelium, unreferenced.
S2 (H149.2)  the same after regressing both on the patient's tumour stromal host (H138.5 statistic).
S3 (H149.3)  malignant vs the tumour stromal host state; companion S3n: adjacent-normal epithelium vs the same host.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER / "exploratory" / "A138_cancer_estimability" / "scripts"))
sys.path.insert(0, str(PAPER / "common" / "datasets"))
import q3_lib as q  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402

from cancer_common import studies as ST  # noqa: E402
from cancer_sc import build as B  # noqa: E402

AMENDMENT = "A149_paired_tumour_normal"
OUT = PAPER / "exploratory" / AMENDMENT / "output"
DATA = ST.DATA
CNV_A149 = DATA / "cnv_a149"
SEED = 149
MIN_PAIRED, MIN_HOST = 4, 5
MAL, NORM, EPI = q.MAL, q.NORM, q.EPI
STROMA = q.STROMA
NEW = {"GSE164789": "lung", "GSE229253": "lung", "GSE153935": "lung", "GSE310916": "lung",
       "GSE200972": "lung", "GSE117570": "lung", "GSE150290": "gastric", "GSE167297": "gastric",
       "GSE199654": "esophageal", "GSE160269": "esophageal", "GSE178481": "kidney", "GSE242299": "kidney",
       "GSE282701": "liver", "GSE242889": "liver", "GSE293189": "bladder", "GSE144240": "skin",
       "GSE193337": "prostate", "GSE176031": "prostate"}
CANCER_ALIAS = {"lung cancer": "lung", "colorectal cancer": "colorectal", "breast cancer": "breast",
                "esophageal cancer": "esophageal", "bladder cancer": "bladder", "liver/biliary cancer": "liver",
                "head and neck cancer": "head and neck", "ovarian cancer": "ovarian", "cervical cancer": "cervical"}
KEEPING = ("lung", "liver", "kidney", "prostate", "thyroid", "breast")  # note 1: TCGA clock-keeping
AT_NULL = ("colorectal", "bladder", "endometrial", "ovarian")
ST.GEO_CANCER.update(NEW)


def cancer_of(study: str) -> str:
    c = NEW.get(study) or (ST.GEO_CANCER.get(study, "") if study.startswith("GSE") else B._inventory_cancer(study))
    return CANCER_ALIAS.get(c, c)


def available_new() -> list:
    return [s for s in NEW if (ST.GEO / s / f"{s}.h5ad").exists()]


def cohorts() -> list:
    return list(B.STUDIES) + available_new()


# ---- cohort tables ---------------------------------------------------------------------------------------------

def new_calls(study: str) -> pd.Series | None:
    """Copy-number calls of a new cohort (`cancer_common.studies.calls`), cached in `$SCCC_DATA/cnv_a149/` (note 1)."""
    f = CNV_A149 / f"{study}.csv.gz"
    if not f.exists():
        c = ST.calls(study)
        if c is None:
            return None
        CNV_A149.mkdir(exist_ok=True)
        c.to_csv(f, index=False)
    return pd.read_csv(f).set_index(["sample", "barcode"]).cnv_call


def table(study: str) -> tuple[pd.DataFrame, str]:
    """Cohort table and malignant compartment name: stored tables for existing cohorts (`q3_lib.study_table`), built from
    counts for new ones with the same rule (called malignant cells where >= 5 patients are called, else tumour
    epithelium). Cached in output/pb/."""
    if study not in NEW:
        return q.study_table(study)
    f = OUT / "pb" / f"{study}.parquet"
    if f.exists():
        pb = pd.read_parquet(f)
        return pb.drop(columns="_mal"), str(pb["_mal"].iloc[0])
    pb, _ = B.build(study, cnv=False)
    mal = EPI
    calls = new_calls(study)
    if calls is not None:
        pc, _ = B.build(study, cnv=True, calls=calls)
        m = pc[pc.compartment == MAL]
        if m.donor.nunique() >= 5:
            pb, mal = pd.concat([pb[pb.compartment != EPI], m], ignore_index=True), MAL
    f.parent.mkdir(parents=True, exist_ok=True)
    pb.assign(_mal=mal).to_parquet(f)
    return pb, mal


def genome_epithelium(study: str, mal: str, pts: list) -> pd.DataFrame:
    """Note 2 (iv): genome-wide counts of the malignant (or tumour epithelium) and adjacent-normal epithelium
    pseudobulks of `pts`, recomputed from the counts with the `cancer_sc` labels (rows: donor x compartment; columns:
    donor, compartment, stratum, depth, genes)."""
    us, _ = ST.units(study)
    us = [u for u in ST.majority_chemistry(us) if u[1] in set(pts)]
    calls = None
    if mal == MAL:
        calls = new_calls(study) if study in NEW else ST.read_calls(study)
    sums, depth, strat, genes = {}, {}, {}, None
    for sample, patient, kind, load, chem, _ in us:
        X, _g, umi, comp, barcodes, sym = load()
        lab = B._labels(comp, kind, sample, barcodes, calls)
        if genes is None:
            genes = pd.Index(sym)
        for c in (mal, NORM):
            i = np.flatnonzero(lab == c)
            if not len(i):
                continue
            k = (patient, c)
            v = np.asarray(X[i].sum(0)).ravel()
            sums[k] = sums.get(k, 0) + v
            depth[k] = depth.get(k, 0) + float(umi[i].sum())
            strat[k] = chem
    keys = sorted(sums)
    G = pd.DataFrame(np.vstack([sums[k] for k in keys]), columns=genes)
    G.columns = [canonical(str(g)) for g in G.columns]
    G = G.T.groupby(level=0).sum().T  # aliases and duplicated symbols summed
    meta = pd.DataFrame({"donor": [k[0] for k in keys], "compartment": [k[1] for k in keys],
                         "stratum": [strat[k] for k in keys], "depth": [depth[k] for k in keys]})
    return pd.concat([meta, G.reset_index(drop=True)], axis=1)


def host_comps(study: str, pb: pd.DataFrame) -> list:
    """Note 2 (i): tumour-sample stromal compartments; learned A17 anchors among them where >= 2, else all of them."""
    pre = [c for c in pb.compartment.unique() if c.startswith("tumour:") and c.split(":")[1] in STROMA]
    try:
        a = [c for c in q.anchors(f"cancer_{study}", "learned") if c in pre]
    except KeyError:
        a = []
    return sorted(a) if len(a) >= 2 else sorted(pre)


# ---- states ----------------------------------------------------------------------------------------------------

def _rows(pb, comp):
    return pb[pb.compartment == comp].drop_duplicates("donor").set_index("donor")


def pair_states(pb: pd.DataFrame, a: str, b: str, wf=q.WF, pts=None) -> dict | None:
    """Centred complex states (and Poisson traces) of compartments a and b on one projector, for the patients having
    both (or `pts`), centred within stratum. None if fewer than MIN_PAIRED patients."""
    A_, B_ = _rows(pb, a), _rows(pb, b)
    common = sorted(set(A_.index) & set(B_.index) if pts is None else set(pts) & set(A_.index) & set(B_.index))
    if len(common) < MIN_PAIRED:
        return None
    A_, B_ = A_.loc[common], B_.loc[common]
    L = projector(pd.concat([A_, B_]), wf, counts=True)
    stratum = A_.stratum.astype(str).to_numpy() if "stratum" in A_ else np.zeros(len(common), int)

    def z(x):
        beta, v = states(x, L, wf, depth=x.depth.to_numpy(float))
        zz = S.complex_state(beta)
        return zz - pd.Series(zz).groupby(stratum).transform("mean").to_numpy(), v

    (za, va), (zb, vb) = z(A_), z(B_)
    return dict(pts=common, stratum=stratum, za=za, vA=va, zb=zb, vB=vb, n=len(common))


def host_states(study: str, pb: pd.DataFrame, pts) -> tuple[np.ndarray, np.ndarray, list] | None:
    """Tumour stromal host per patient (centred within stratum later), its posterior trace, and the patients."""
    hc = host_comps(study, pb)
    if not hc:
        return None
    h = q.host_of(pb, hc)
    if h is None:
        return None
    keep = [p for p in pts if p in h.index]
    if len(keep) < MIN_HOST:
        return None
    H = h.loc[keep]
    return (H.x + 1j * H.y).to_numpy(), (H.xx + H.yy).to_numpy(), keep


def _centre(z, stratum):
    return z - pd.Series(z).groupby(stratum).transform("mean").to_numpy()


def units(study: str) -> dict:
    """S1, S2, S3 and S3n units of one cohort (entries None when not estimable at the design)."""
    pb, mal = table(study)
    out = dict(study=study, cancer=cancer_of(study), malignant=mal, s1=None, s2=None, s3=None, s3n=None)
    if mal not in set(pb.compartment) or NORM not in set(pb.compartment):
        return out
    u = pair_states(pb, mal, NORM)
    if u is None:
        return out
    out["s1"] = dict(za=u["za"], vA=u["vA"], zb=u["zb"], vB=u["vB"], n=u["n"])
    hs = host_states(study, pb, u["pts"])
    if hs is None:
        return out
    h, vh, keep = hs
    idx = [u["pts"].index(p) for p in keep]
    st = u["stratum"][idx]
    za, va, zb, vb = (_centre(u["za"][idx], st), u["vA"][idx], _centre(u["zb"][idx], st), u["vB"][idx])
    hc = _centre(h, st)
    ra, vra = q.residual(za, va, hc, vh)
    rb, vrb = q.residual(zb, vb, hc, vh)
    out["s2"] = dict(za=ra, vA=vra, zb=rb, vB=vrb, n=len(keep))
    out["s3"] = dict(za=za, vA=va, zb=hc, vB=vh, n=len(keep))
    out["s3n"] = dict(za=zb, vA=vb, zb=hc, vB=vh, n=len(keep))
    # S3 for all malignant patients with a host (not only the paired ones) is not used: one patient set per cohort
    return out


# ---- statistics ------------------------------------------------------------------------------------------------

def null_of(u: dict, n_perm: int, rng) -> np.ndarray:
    return q.perm_null(u["za"], u["vA"], u["zb"], u["vB"], n_perm, rng)


def obs_of(u: dict) -> float:
    return q.corr_stat(u["za"], u["vA"], u["zb"], u["vB"])


def estimable(u: dict, null: np.ndarray) -> bool:
    return q.estimable(u["za"], u["vA"], u["zb"], u["vB"], null)


def boot_of(u: dict, n_boot: int, rng) -> np.ndarray:
    """Statistic on patients resampled with replacement (note 2 v)."""
    n, out = u["n"], np.empty(n_boot)
    for b in range(n_boot):
        i = rng.integers(0, n, n)
        za, zb = u["za"][i] - u["za"][i].mean(), u["zb"][i] - u["zb"][i].mean()
        out[b] = q.corr_stat(za, u["vA"][i], zb, u["vB"][i])
    return out


pool = q.pool
summarise = q.summarise
