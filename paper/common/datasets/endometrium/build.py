"""Human endometrium (GSE266265; endometriosis vs control, scRNA-seq of biopsies): donor x class pseudobulks.

Source: GEO GSE266265 (Seurat object parsed without R by desync/19_endometrium/scripts/stream_rds.py, then grouped by
desync/19_endometrium/scripts/build_endometrium.py into the group pseudobulk), raw input
$SCCC_DATA/geo/GSE266265/pseudobulk. Amendment A11 part 3. A17 targets: `endometrium` (healthy fit: control donors,
classes with >= 20 of them), `endometriosis` (case = Census-style disease "endometriosis"), `endometriosis_severe`
(case = grade "severe"); donor columns `healthy`, `endometriosis`, `endometriosis_severe`.
Legacy builder replaced: desync/19_endometrium/scripts/endometrium.py.

Choices: classes `CLASSES` (exact names), compartment `endometrium:<class>`; stratum = `10x|<menstrual phase>` (the
group file's assay column); pseudobulks >= 20 cells; donors with disease "normal" are controls.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import census as C  # noqa: E402
import legacy as L  # noqa: E402
import organ as O  # noqa: E402

DATASET = DS = "gse266265"
PATH = C.GEO / "GSE266265" / "pseudobulk"
CLASSES = tuple((c, f"^{c}$") for c in ("stromal", "epithelial", "endothelial", "myeloid", "lymphocyte"))
SPEC = O.Spec(PATH, DS, "endometrium", CLASSES)
REF = "endometrium:stromal"


def _donors():
    d = C.census_donors(PATH, DS)
    return d


def _labels(d):
    ctrl = d.index[d.disease == "normal"]
    out = {}
    for name, case in (("endometriosis", d.disease == "endometriosis"), ("endometriosis_severe", d.grade == "severe")):
        lab = pd.Series("control", index=ctrl).combine_first(pd.Series("case", index=d.index[case]))  # case wins
        out[name] = lab[(lab == "control") | case.reindex(lab.index).fillna(False).to_numpy()]
    return out


def build(donors=None, genes="clock"):
    d = _donors()
    labs = _labels(d)
    keep = labs["endometriosis"].index
    if donors is not None:
        keep = keep[keep.isin([f"{DS}:{x}" for x in donors])]
    pb = O.table(O.groups(SPEC, donors, genes), keep, 20)
    extra = {k: v.rename(k) for k, v in labs.items()}
    if donors is None:
        extra["healthy"] = O.flag(O.healthy_members(pb, d.index[d.disease == "normal"], 20))
    return O.finish(pb, d, extra, cols=("disease", "grade"))


def from_legacy():
    d = _donors()
    return L.from_pickles(["endometrium", "endometriosis", "endometriosis_severe"],
                          covariates=C.donors_frame(d, d.index, cols=("disease", "grade")))


def check_legacy(ds, pb, donors):
    d = _donors()
    return dict(healthy_rule_ok=set(O.healthy_members(pb, d.index[d.disease == "normal"], 20)) ==
                set(donors.index[donors.healthy]))
