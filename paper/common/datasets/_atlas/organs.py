"""Healthy organ atlases of A9 part 2 (breast, bone marrow, adipose, muscle, ileum, colon, optic nerve, anterior eye):
shared builder of the eight `organ_<name>/build.py`.

Legacy: desync/13_healthy_organs/scripts/organs.py (`ORGANS`, `run`) with `run_atlas.COARSE` classes, captured as A17
targets `organ_<name>` (healthy fits). All donors of the dataset are `healthy`; ileum_kong and colon_crc keep only
Census disease "normal" groups.

Choices: classes `COARSE` (case-insensitive, first match), compartment `<organ>:<class>` (with a region map
`<organ>:<region>_<class>` for adipose, optic nerve, anterior eye); majority assay; pseudobulks >= 20 cells; classes
with >= 20 donors.
"""

from __future__ import annotations

import census as C
import legacy as L
import organ as O


def _normal(o):
    return o.disease.astype(str).to_numpy() == "normal"


ORGANS = {  # name: (data dir, region map or None, group filter or None)
    "breast": ("breast_kumar", None, None),
    "bone_marrow": ("bone_marrow", None, None),
    "adipose": ("adipose_emont", {"subcutaneous adipose tissue": "subcutaneous", "omental fat pad": "omental"}, None),
    "muscle": ("muscle", None, None),
    "ileum": ("ileum_kong", None, _normal),
    "colon": ("colon_crc", None, _normal),
    "optic_nerve": ("optic_nerve", {"cranial nerve II": "nerve", "optic disc": "head"}, None),
    "eye_anterior": ("eye_anterior", {"ciliary body": "ciliary", "eye trabecular meshwork": "trabecular",
                                      "uvea": "uvea"}, None),
}


def spec(organ: str) -> O.Spec:
    ds, region, keep = ORGANS[organ]
    return O.Spec(C.ATLAS / ds, ds, organ, O.COARSE, keep=keep, region=region)


def build(name: str, donors=None, genes="clock"):
    organ = name.removeprefix("organ_")
    s = spec(organ)
    F = O.groups(s, donors, genes)
    pb = O.table(F, F.resection.unique(), 20, min_donors=None if donors is not None else 20)
    return O.finish(pb, C.census_donors(s.path, s.prefix), {"healthy": O.flag(pb.resection.unique())})


def from_legacy(name: str):
    s = spec(name.removeprefix("organ_"))
    d = C.census_donors(s.path, s.prefix)
    return L.from_pickles([name], covariates=C.donors_frame(d, d.index))
