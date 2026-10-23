"""Pancreatic islets, T2D vs control, with the A96 gene panel (core clock + cue + locked genes; panel of the h5ad): donor x class pseudobulks.

Raw input $SCCC_DATA/atlas/islets_t2d_a96 (the Census build of islets_t2d rerun with the A96 panel). Amendment A96 part 2.
Legacy builder replaced: desync/validation/masld/scripts/a96_counts.py (`islets`). Same rules as `islets` (classes, >= 20 cells,
Census disease normal / type 2 diabetes mellitus donors); donor column `arm` (control / case); gene columns are canonical symbols.
"""

import importlib.util
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "_atlas"))
import census as C  # noqa: E402
import organ as O  # noqa: E402

DATASET, DS = "islets_a96", "islets_t2d"
_spec = importlib.util.spec_from_file_location("build_islets", HERE.parent / "islets" / "build.py")
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)
PATH = C.ATLAS / "islets_t2d_a96"


sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/
from common.paths import ARCHIVE  # noqa: E402
def build(donors=None, genes="all"):
    from clocksync.genes import canonical

    d = C.census_donors(C.ATLAS / DS, DS)
    lab = _base._lab(d)
    if donors is not None:
        lab = lab[lab.index.isin([f"{DS}:{x}" for x in donors])]
    spec = O.Spec(PATH, DS, "islet", _base.CLASSES)
    pb = O.table(O.groups(spec, donors, genes), lab.index, 20)
    pb = pb.rename(columns=lambda c: canonical(c) if c not in ("resection", "compartment", "stratum", "n_cells", "depth") else c)
    return O.finish(pb, d, {"arm": lab.rename("arm")})


def from_legacy():
    f = ARCHIVE / "legacy" / "desync" / "validation" / "masld" / "output" / "a96_counts.csv"
    t = pd.read_csv(f)
    t = t[t.dataset == "islets"].dropna(axis=1, how="all").rename(columns={"resection": "donor"})
    arm = t.groupby("donor").arm.first()
    pb = C.standard(t.drop(columns=["dataset", "arm"]).rename(columns={"donor": "resection"}))
    d = C.census_donors(C.ATLAS / DS, DS)
    out = C.donors_frame(d, pb.donor.unique(), {"arm": arm})
    return pb, out
