"""Tabula Sapiens stromal-vascular pseudobulks per donor x organ, all genes (A109b: replication of organ coupling).

Raw input: `$SCCC_DATA/atlas_wg/tabula_sapiens` if present (genome-wide group pseudobulk), else the all-gene group
pseudobulk `$SCCC_DATA/atlas/tabula_sapiens` (Census 2025-11-08, dataset 53d208b0..., 61497 genes; no download needed).
Choices: 10x groups only; cell classes `COARSE`; stromal-vascular = endothelial + fibroblast_stromal + mural_smooth_muscle
(the A17 anchor families), pooled into one compartment `<tissue_general>:stromal_vascular`; majority assay per donor;
pseudobulks >= `MIN_CELLS` = 50 cells (the Tabula Sapiens rule of `tabula_sapiens/build.py`); organs with >= 5 donors.
Genes: canonical symbols matched on Census `feature_name`; all genes with a count are kept (pseudo-clock null).
Output: `$SCCC_DATA/pseudobulk/tabula_sapiens_organs.parquet` (+ `_donors`), `stats.csv` next to the builder's log.

    python build.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import census as C  # noqa: E402
import organ as O  # noqa: E402
from common import data  # noqa: E402

DATASET = "tabula_sapiens_organs"
SV = {"endothelial", "fibroblast_stromal", "mural_smooth_muscle"}
MIN_CELLS, MIN_DONORS = 50, 5


def build(donors=None, min_cells: int = MIN_CELLS):
    path = C.ATLAS_WG / "tabula_sapiens"
    path = path if (path / "pseudobulk.h5ad").exists() else C.ATLAS / "tabula_sapiens"
    o, X, names = C.read_sparse(path, donors)
    cls = C.classify(o.cell_type.astype(str), O.COARSE, ignore_case=False)
    sv = cls.isin(SV).to_numpy() & o.assay.astype(str).str.startswith("10x").to_numpy()
    comp = o.tissue_general.astype(str) + ":stromal_vascular"
    F = C.keep_majority(C.frame(o, pd.DataFrame(index=range(len(o))), "tabula_sapiens", comp, o.assay, keep=sv))
    pb = C.sum_groups_sparse(F, X, names, None, min_cells, "batch", None if donors is not None else None)
    n = pb.groupby("compartment").resection.nunique()
    pb = pb[pb.compartment.isin(n[n >= MIN_DONORS].index)].reset_index(drop=True)
    pb["stratum"] = "tabula_sapiens|" + pb.stratum
    d = C.census_donors(path, "tabula_sapiens")
    pb, out = O.finish(pb, d, {"healthy": O.flag(pb.resection.unique())})
    missing = [g for g in C.GENES if g not in pb.columns]
    print(f"clock genes found: {len(C.GENES) - len(missing)}/{len(C.GENES)}; missing {missing}; genes {pb.shape[1] - 5}")
    return pb, out


if __name__ == "__main__":
    pb, d = build()
    data.save(DATASET, pb, d)
    print(pb.shape)
