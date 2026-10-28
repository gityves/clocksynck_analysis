"""Kidney single-cell atlas (KPMP, scRNA-seq; healthy reference and chronic kidney disease): donor x class pseudobulks.

Source: CELLxGENE Census KPMP single-cell dataset, raw input $SCCC_DATA/atlas/kpmp_sc (group pseudobulk). Amendment A8.
Legacy builder replaced: desync/11_kidney/scripts/kidney.py (`kgroups`, `labels`, `block`, `ckd`, `healthy`).
A17 targets: `kidney_sc` (healthy fit: control donors, classes with >= 20 of them) and `ckd_sc` (donor column `ckd_sc`).

Choices: classes `CLASSES` (case-insensitive except the T pattern, first match), compartment `kidney:<class>`; majority assay; stratum =
assay | `dataSource` (KPMP / UMICH processing site, donor-constant); pseudobulks >= 20 cells; donors = Census disease
"normal" (control) and "chronic kidney disease" (case); one table for both (the healthy table is its control part).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import census as C  # noqa: E402
import legacy as L  # noqa: E402
import organ as O  # noqa: E402

DATASET, DS = "kidney_sc", "kpmp_sc"
CLASSES = (
    ("proximal_tubule", r"proximal tubule"), ("loop_of_Henle", r"loop of Henle"),
    ("distal", r"distal convoluted|connecting tubule"), ("collecting_duct", r"collecting duct"),
    ("podocyte", r"podocyte|parietal epithelial"), ("endothelial", r"endothelial"),
    ("myeloid", r"macrophage|monocyte|dendritic|mononuclear phagocyte"),
    ("interstitial", r"interstitial fibroblast|^kidney interstitial cell$"), ("pericyte", r"pericyte"),
    ("T", r"(?-i:\bT cell)"), ("B", r"^B cell$"),  # case-sensitive: "mast cell" is not a T cell
)
SPEC = O.Spec(C.ATLAS / DS, DS, "kidney", CLASSES)
REF = "kidney:endothelial"


def build(donors=None, genes="clock"):
    d = C.census_donors(C.ATLAS / DS, DS)
    lab = O.labels_from(d, "disease", ["normal"], ["chronic kidney disease"])
    if donors is not None:
        lab = lab[lab.index.isin([f"{DS}:{x}" for x in donors])]
    pb = O.table(O.groups(SPEC, donors, genes), lab.index, 20, block=d["dataSource"].astype(str))
    extra = {"ckd_sc": lab.rename("ckd_sc")}
    if donors is None:
        extra["healthy"] = O.flag(O.healthy_members(pb, lab.index[lab == "control"], 20))
    return O.finish(pb, d, extra, cols=("age", "sex", "disease", "dataSource"))


def from_legacy():
    d = C.census_donors(C.ATLAS / DS, DS)
    return L.from_pickles(["kidney_sc", "ckd_sc"], covariates=C.donors_frame(d, d.index, cols=("age", "sex", "disease", "dataSource")))


def check_legacy(ds, pb, donors):
    d = C.census_donors(C.ATLAS / DS, DS)
    lab = O.labels_from(d, "disease", ["normal"], ["chronic kidney disease"])
    h = O.healthy_members(pb, lab.index[lab == "control"], 20)
    return dict(healthy_rule_ok=set(h) == set(donors.index[donors.healthy]))
