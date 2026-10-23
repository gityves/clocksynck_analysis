"""Pancreatic islets, type 2 diabetes vs control (Census islet atlas e51bae9a..., "islets_t2d"): donor x class pseudobulks.

Raw input $SCCC_DATA/atlas/islets_t2d (group pseudobulk). Amendment A13. A17 targets: `islets` (healthy fit: control
donors, classes with >= 15 of them; donor column `healthy`) and `islets_t2d` (donor column `islets_t2d`).
Legacy builder replaced: desync/15_islets/scripts/islets.py (`X`, `labels`, `healthy`, `t2d`).

Choices: classes `CLASSES` (case-insensitive, first match), compartment `islet:<class>`; majority assay; pseudobulks
>= 20 cells; donors = Census disease "normal" (control) or "type 2 diabetes mellitus" (case; prediabetic donors are
not used).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import census as C  # noqa: E402
import legacy as L  # noqa: E402
import organ as O  # noqa: E402

DATASET = DS = "islets_t2d"
CLASSES = (("beta", r"type B pancreatic"), ("alpha", r"pancreatic A cell"), ("delta", r"pancreatic D cell"),
           ("PP", r"pancreatic PP cell"), ("acinar", r"acinar"), ("ductal", r"ductal"), ("endothelial", r"endothelial"),
           ("stellate", r"stellate"), ("immune", r"leukocyte|mast cell"))
SPEC = O.Spec(C.ATLAS / DS, DS, "islet", CLASSES)
REF = "islet:alpha"


def _lab(d):
    return O.labels_from(d, "disease", ["normal"], ["type 2 diabetes mellitus"])


def build(donors=None, genes="clock"):
    d = C.census_donors(C.ATLAS / DS, DS)
    lab = _lab(d)
    if donors is not None:
        lab = lab[lab.index.isin([f"{DS}:{x}" for x in donors])]
    pb = O.table(O.groups(SPEC, donors, genes), lab.index, 20)
    extra = {"islets_t2d": lab.rename("islets_t2d")}
    if donors is None:
        extra["healthy"] = O.flag(O.healthy_members(pb, lab.index[lab == "control"], 15))
    return O.finish(pb, d, extra)


def from_legacy():
    d = C.census_donors(C.ATLAS / DS, DS)
    return L.from_pickles(["islets", "islets_t2d"], covariates=C.donors_frame(d, d.index))


def check_legacy(ds, pb, donors):
    d = C.census_donors(C.ATLAS / DS, DS)
    lab = _lab(d)
    return dict(healthy_rule_ok=set(O.healthy_members(pb, lab.index[lab == "control"], 15)) ==
                set(donors.index[donors.healthy]))
