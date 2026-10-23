"""Mouse colon, pre-publication chrono-atlas (8 animals x proximal / distal colon, ZT 1, 7, 13, 19), 14 clock genes. Input
`/mnt/max8tb/data/CoPhaser_data/ciradian_intestin/270525_colon_data.h5ad` (Seurat export; X = counts; override with $SCCC_COLON_PREPUB).
Used by: V13 (cell-type coupling against known time; `simulations/compute_celltype_truth.py`). The file holds the
regenerating colon (R-DC + R-PC), not healthy colon (A21 note); V13 is reported as such.

Table: donor = animal (first field of `orig.ident`), compartment = `<region>:<class>` (region proximal / distal from `cond` PC / DC;
author `celltype` mapped to the A17 class names, others dropped), `stratum` = region (one row per region, so a donor has two strata:
analyses use `Dataset(..., mixed_strata=True)`), n_cells, depth (nCount_RNA summed), genes. Donor table: `zt`.
All pseudobulks with a cell are kept (the V13 thresholds, >= 30 cells and a class in >= 75 % of the animals, are applied by the analysis).
Legacy builder replaced: legacy/scripts/v13_celltype_truth.py::load("colon").
"""

import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import clocksync as cs  # noqa: E402

SRC = Path(os.environ.get("SCCC_COLON_PREPUB", "/mnt/max8tb/data/CoPhaser_data/ciradian_intestin/270525_colon_data.h5ad"))
CLASSES = {  # author `celltype` -> A17 class name; unlisted types are dropped
    "Interstitial": "fibroblast", "Trophocytes": "fibroblast", "Telocytes": "fibroblast", "Reticular": "fibroblast",
    "MM": "smooth_muscle", "MP": "smooth_muscle", "LPM": "smooth_muscle", "Pericytes": "pericyte",
    "Endothelial": "endothelial", "Endothelial_proliferative": "endothelial", "Lymphatic": "lymphatic_EC",
    "T": "T", "T_proliferative": "T", "B": "B", "B_proliferative": "B", "Plasma": "plasma",
    "Macrophages": "macrophage", "Dendritic": "DC",
    "Colonocytes": "colonocyte", "Colonocytes_TA": "colonocyte", "Goblet": "secretory", "Goblet_TA": "secretory",
    "ISC": "stem", "ISC_proliferative": "stem",
}


def build(donors=None):
    o = ad.read_h5ad(SRC, backed="r").obs
    region = o["cond"].astype(str).map({"PC": "proximal", "DC": "distal"})
    cls = o["celltype"].astype(str).map(CLASSES)
    lab = pd.DataFrame({"donor": o["orig.ident"].astype(str).str.split("_").str[0], "compartment": region + ":" + cls,
                        "region": region, "nCount_RNA": o["nCount_RNA"].astype(float)}, index=o.index.astype(str))
    lab = lab[cls.notna().to_numpy()]
    pb = cs.pseudobulk(SRC, "donor", "compartment", labels=lab, depth="nCount_RNA", min_cells=1)
    pb["stratum"] = pb.compartment.str.split(":").str[0]
    zts = pd.Series(o["ZT"].astype(str).str[2:].astype(float).to_numpy(), index=o["orig.ident"].astype(str).str.split("_").str[0].to_numpy())
    dt = pd.DataFrame({"zt": zts.groupby(level=0).first()}).rename_axis("donor")
    if donors is not None:
        pb = pb[pb.donor.isin(donors)]
        dt = dt.loc[dt.index.isin(donors)]
    return pb.reset_index(drop=True), dt


if __name__ == "__main__":
    from common import data
    pb, dt = build()
    data.save("colon_prepub", pb, dt)
    print(pb.shape, dt.shape)
