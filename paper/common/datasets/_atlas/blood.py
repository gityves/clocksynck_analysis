"""Blood PBMC atlases: OneK1K, AIDA and the SLE cohort of Perez et al. 2022: shared builder of `blood_onek1k`,
`blood_aida`, `blood_lupus`. Amendments A6 (d), A9.

Legacy: desync/09_blood_onek1k/scripts/blood.py (OneK1K: exact cell-type map, stratum = processing pool) and
blood_a9.py (AIDA, lupus: regex classes, `BLOCKS`). A17 targets `blood_onek1k`, `blood_aida` (healthy fits) and `lupus`
(SLE vs control, donor column `blood_lupus`).
Raw inputs $SCCC_DATA/atlas/{onek1k,aida_a7,lupus_a7} (Census group pseudobulks).

Choices. OneK1K: cell types mapped by `ONEK1K` (CD4 T, CD8 T, NK, B, CD14 and CD16 monocytes; others dropped), no
assay filter, stratum = `pool_number`; pseudobulks >= 50 cells; classes with >= 20 donors. AIDA / lupus: classes
`CLASSES` (case-insensitive, first match), compartment `blood:<class>`, majority assay (stratum = assay | site columns of
`BLOCKS` found in donors.csv: none in these files); pseudobulks >= 50 cells; AIDA classes with >= 20 donors; lupus keeps all
classes for donors with Census disease "normal" (control) or "systemic lupus erythematosus" (case).
"""

from __future__ import annotations

import pandas as pd

import census as C
import legacy as L
import organ as O

ONEK1K = {
    "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4_T",
    "central memory CD4-positive, alpha-beta T cell": "CD4_T",
    "effector memory CD4-positive, alpha-beta T cell": "CD4_T",
    "CD4-positive, alpha-beta cytotoxic T cell": "CD4_T", "CD4-positive, alpha-beta T cell": "CD4_T",
    "regulatory T cell": "CD4_T",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "CD8_T",
    "central memory CD8-positive, alpha-beta T cell": "CD8_T",
    "effector memory CD8-positive, alpha-beta T cell": "CD8_T", "CD8-positive, alpha-beta T cell": "CD8_T",
    "natural killer cell": "NK", "CD16-negative, CD56-bright natural killer cell, human": "NK",
    "naive B cell": "B", "memory B cell": "B", "transitional stage B cell": "B",
    "CD14-positive monocyte": "CD14_monocyte", "CD14-low, CD16-positive monocyte": "CD16_monocyte"}
CLASSES = (
    ("CD4_T", r"CD4-positive|regulatory T"), ("CD8_T", r"CD8-positive"),
    ("ncMono", r"CD16-positive monocyte|non-classical monocyte"), ("cMono", r"CD14-positive monocyte|classical monocyte"),
    ("NK", r"natural killer"), ("B", r"^(naive |memory |mature )?B cell$"),
    ("cDC", r"conventional dendritic|CD1c-positive myeloid dendritic|CD141-positive myeloid dendritic"),
)
BLOCKS = {"aida_a7": ("Country", "institute", "country", "site"), "lupus_a7": ("Processing_Cohort",)}
DS = {"blood_onek1k": "onek1k", "blood_aida": "aida_a7", "blood_lupus": "lupus_a7"}
LEGACY = {"blood_onek1k": ["blood_onek1k"], "blood_aida": ["blood_aida"], "blood_lupus": ["lupus"]}


def _block(d: pd.DataFrame, ds: str):
    cols = [c for c in BLOCKS[ds] if c in d]
    return d[cols].astype(str).agg("|".join, axis=1) if cols else None


def _lupus_labels(d):
    return O.labels_from(d, "disease", ["normal"], ["systemic lupus erythematosus"])


def build(name: str, donors=None, genes="clock"):
    ds = DS[name]
    path = C.ATLAS / ds
    d = C.census_donors(path, ds)
    if name == "blood_onek1k":
        o, X = C.read(path, genes, donors)
        cls = o.cell_type.astype(str).map(ONEK1K)
        F = C.frame(o, X, ds, "blood:" + cls.astype(str), pd.Series("", index=o.index), keep=cls.notna().to_numpy())
        pb = C.donor_table(F, None, 50, stratum=None)
        pb["stratum"] = pb.resection.str.removeprefix(f"{ds}:").map(pd.read_csv(path / "donors.csv", dtype={"donor_id": str})
                                                                    .set_index("donor_id").pool_number).astype(str).to_numpy()
        if donors is None:
            n = pb.groupby("compartment").resection.nunique()
            pb = pb[pb.compartment.isin(n[n >= 20].index)]
        extra = {"healthy": O.flag(pb.resection.unique())}
        return O.finish(pb, d, extra)
    F = O.groups(O.Spec(path, ds, "blood", CLASSES), donors, genes)
    block = _block(d, ds)
    if name == "blood_aida":
        pb = O.table(F, F.resection.unique(), 50, block, min_donors=None if donors is not None else 20)
        return O.finish(pb, d, {"healthy": O.flag(pb.resection.unique())})
    lab = _lupus_labels(d)
    lab = lab[lab.index.isin(F.resection.unique())]
    pb = O.table(F, lab.index, 50, block)
    return O.finish(pb, d, {"lupus": lab.rename("lupus")})


def from_legacy(name: str):
    ds = DS[name]
    d = C.census_donors(C.ATLAS / ds, ds)
    return L.from_pickles(LEGACY[name], covariates=C.donors_frame(d, d.index))
