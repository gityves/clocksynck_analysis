"""Skin atopic dermatitis (GSE328048, 10x 5') and systemic sclerosis (GSE195452, MARS-seq): shared builder of `skin_ad`,
`skin_ad_nonlesional` and `skin_ssc`. Amendment A11 parts 1-2.

Legacy: desync/14_skin/scripts/skin.py (`X`, `healthy`, `ssc`, `ad_vs_control`), grouping by
desync/14_skin/scripts/build_ssc.py and build_ad.py (GEO matrices; the group files are the raw input here:
$SCCC_DATA/geo/GSE328048/pseudobulk and $SCCC_DATA/geo/GSE195452/pseudobulk).
A17 targets: `skin_ad_control` (healthy), `skin_ad` (lesional vs control), `skin_ad_nonlesional`, `skin_ssc_control`
(healthy), `skin_ssc`. Donor columns: `healthy`, `skin_ad` / `skin_ad_nonlesional` / `skin_ssc` (case/control).

Choices: classes `CLASSES` (exact names), compartment `skin:<class>`; majority assay; pseudobulks >= 20 cells; AD:
control-tissue groups of normal donors plus lesional (or nonlesional) groups of AD donors, majority assay per sample
type; SSc: all donors, control = Census disease "normal", case = any other label; the healthy fits keep classes with
>= 20 normal donors. A donor's lesional and nonlesional samples are two different tables (the paired layout of A54 is
not a standard table).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import census as C
import legacy as L
import organ as O

CLASSES = tuple((c, f"^{c}$") for c in ("fibroblast", "pericyte", "endothelial", "T", "NK", "myeloid", "mast",
                                        "keratinocyte", "B_plasma", "melanocyte"))
AD = ("gse328048_ad", C.GEO / "GSE328048" / "pseudobulk")
SSC = ("gur2022_ssc", C.GEO / "GSE195452" / "pseudobulk")
LEGACY = {"skin_ad": ["skin_ad_control", "skin_ad"], "skin_ad_nonlesional": ["skin_ad_nonlesional"],
          "skin_ssc": ["skin_ssc_control", "skin_ssc"]}


def _spec(ds, path, cond=None):
    keep = None if cond is None else (lambda o: o.tissue.astype(str).to_numpy() == cond)
    return O.Spec(path, ds, "skin", CLASSES, keep=keep)


def _labels(d):
    return pd.Series(np.where(d.disease.astype(str) == "normal", "control", "case"), index=d.index)


def build(name: str, donors=None, genes="clock"):
    ds, path = SSC if name == "skin_ssc" else AD
    d = C.census_donors(path, ds)
    lab = _labels(d)
    if name == "skin_ssc":
        F = O.groups(_spec(ds, path), donors, genes)
        lab = lab[lab.index.isin(F.resection.unique())] if donors is not None else lab
        pb = O.table(F, lab.index, 20)
    else:
        case = "lesional" if name == "skin_ad" else "nonlesional"
        F = pd.concat([O.groups(_spec(ds, path, c), donors, genes) for c in ("control", case)], ignore_index=True)
        lab = lab[lab.index.isin(F.resection.unique())]
        pb = O.table(F, F.resection.unique(), 20)
    extra = {name: lab.rename(name)}
    if donors is None and name != "skin_ad_nonlesional":
        extra["healthy"] = O.flag(O.healthy_members(pb, d.index[d.disease.astype(str) == "normal"], 20))
    return O.finish(pb, d, extra, cols=("disease",))


def from_legacy(name: str):
    ds, path = SSC if name == "skin_ssc" else AD
    d = C.census_donors(path, ds)
    return L.from_pickles(LEGACY[name], covariates=C.donors_frame(d, d.index, cols=("disease",)))


def check_legacy(name, pb, donors):
    ds, path = SSC if name == "skin_ssc" else AD
    d = C.census_donors(path, ds)
    if "healthy" not in donors:
        return {}
    return dict(healthy_rule_ok=set(O.healthy_members(pb, d.index[d.disease.astype(str) == "normal"], 20)) ==
                set(donors.index[donors.healthy]))
