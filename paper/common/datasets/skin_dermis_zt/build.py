"""Mouse dorsal skin single cells over the day (GSE223109, 18 animals, ZT 2..70 every 4 h), 14 clock genes. Input
`/mnt/max8tb/data/CoPhaser_data/training_set/circadian_cycle/skin.h5ad` (layer `counts`; override with $SCCC_SKIN_ZT).
Used by: V13 (cell-type coupling against known time; `simulations/compute_celltype_truth.py`).

Table: donor = `sample_name`, compartment = author `cell_type` (fibroblast and immune among others), stratum '', n_cells, depth (`n_counts`),
genes. Donor table: `zt` (`ZTmod` mod 24, h). All pseudobulks with a cell are kept. The split of the fibroblasts into exchangeable
anchor compartments (V13: random parts of each animal's fibroblasts) is made by the analysis, which needs cell-level labels:
`fibroblast_parts(k, seed)` returns the label table for `clocksync.pseudobulk`.
Legacy builder replaced: legacy/scripts/v13_celltype_truth.py::load("skin").
"""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import clocksync as cs  # noqa: E402

SRC = Path(os.environ.get("SCCC_SKIN_ZT", "/mnt/max8tb/data/CoPhaser_data/training_set/circadian_cycle/skin.h5ad"))


def labels(parts: int = 0, seed: int = 0) -> pd.DataFrame:
    """Cell labels (donor, compartment, n_counts, zt). `parts` >= 2: each fibroblast is assigned at random to one of
    `parts` compartments `fibroblast_a`, `fibroblast_b`, ... (exchangeable anchors)."""
    o = ad.read_h5ad(SRC, backed="r").obs
    cls = o["cell_type"].astype(str).to_numpy().copy()
    if parts >= 2:
        fib = cls == "fibroblast"
        rng = np.random.default_rng(seed)
        if parts == 2:  # as the legacy: one uniform draw per fibroblast, < 0.5 -> a
            cls[fib] = np.where(rng.random(fib.sum()) < 0.5, "fibroblast_a", "fibroblast_b")
        else:
            cls[fib] = np.array([f"fibroblast_{chr(97 + i)}" for i in range(parts)])[rng.integers(0, parts, fib.sum())]
    return pd.DataFrame({"donor": o["sample_name"].astype(str).to_numpy(), "compartment": cls,
                         "n_counts": o["n_counts"].astype(float).to_numpy(), "zt": o["ZTmod"].astype(float).to_numpy() % 24}, index=o.index.astype(str))


def build(donors=None, parts: int = 0, seed: int = 0):
    lab = labels(parts, seed)
    pb = cs.pseudobulk(SRC, "donor", "compartment", labels=lab, layer="counts", depth="n_counts", min_cells=1)
    pb["stratum"] = ""
    dt = pd.DataFrame({"zt": lab.groupby("donor").zt.first()})
    if donors is not None:
        pb, dt = pb[pb.donor.isin(donors)], dt.loc[dt.index.isin(donors)]
    return pb.reset_index(drop=True), dt


if __name__ == "__main__":
    pb, dt = build()
    from common import data
    data.save("skin_dermis_zt", pb, dt)
    print(pb.shape, dt.shape)
