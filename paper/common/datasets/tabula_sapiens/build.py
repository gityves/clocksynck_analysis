"""Tabula Sapiens (14 deceased donors, all organs; Census tabula_sapiens): donor x organ x cell class pseudobulks.

Raw input $SCCC_DATA/atlas/tabula_sapiens (Census group pseudobulk). A17 step 7 (confirmation only; excluded from the
anchor-learning set). A17 target `tabula_sapiens`.
Legacy builder replaced: desync/common/run_atlas.py::load(dataset, 50, "coarse") + a17_tables.py::tabula.

Choices: 10x groups; classes `COARSE` (case-sensitive, first match; unmatched dropped); compartment
`<tissue_general>:<class>`; majority assay per donor; stratum `tabula_sapiens|<assay>`; pseudobulks >= 50 cells; classes
with >= 5 donors.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import census as C  # noqa: E402
import legacy as L  # noqa: E402
import organ as O  # noqa: E402

DATASET = DS = "tabula_sapiens"
MIN_CELLS, MIN_DONORS = 50, 5


def build(donors=None, genes="clock"):
    o, X = C.read(C.ATLAS / DS, genes, donors)
    cell = C.classify(o.cell_type.astype(str), O.COARSE, ignore_case=False)
    umi = o.assay.astype(str).str.startswith("10x").to_numpy()
    F = C.keep_majority(C.frame(o, X, DS, o.tissue_general.astype(str) + ":" + cell.astype(str), o.assay,
                                keep=umi & cell.notna().to_numpy()))
    pb = C.donor_table(F, None, MIN_CELLS, stratum="batch", min_donors=None if donors is not None else MIN_DONORS)
    pb["stratum"] = DS + "|" + pb.stratum
    return O.finish(pb, C.census_donors(C.ATLAS / DS, DS), {"healthy": O.flag(pb.resection.unique())})


def from_legacy():
    d = C.census_donors(C.ATLAS / DS, DS)
    return L.from_pickles([DS], covariates=C.donors_frame(d, d.index))
