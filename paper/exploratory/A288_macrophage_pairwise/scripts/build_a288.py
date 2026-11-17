"""A288 tables: IBDverse macrophage rows (registered author labels, healthy cells) and the HLCA macrophage split
(alveolar against monocyte-derived, Census `cell_type`). Registered in paper/docs/registered_experiments.md (A288).

    python -I build_a288.py ibd RAW_H5AD OUT_DIR     # raw IBDverse -> group pseudobulk per macrophage author label
    (the score script calls `ibd_macrophage(...)` and `hlca_split()`)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common/datasets/_atlas"))
sys.path.insert(0, str(PAPER / "exploratory/A284_residency/scripts"))

IBD_LABELS = {"mac_CD163hi": {"Macrophage CD163hi"}, "mac_CD163lo": {"Macrophage CD163lo"},
              "mac_CCL": {"Macrophage CCL3/4 CXCL9/10"}, "mac_MMP9": {"Macrophage MMP9"},
              "mac_prolif": {"Macrophage proliferating"}}
IBD_MAC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/atlas_wg/ibdverse_mac")
HLCA_SPLIT = {"alveolar macrophage": "macrophage_alveolar", "elicited macrophage": "macrophage_monocyte_derived"}
MIN_CELLS = 20


def ibd_macrophage(tissue, donors, genes):
    """Healthy donor rows of `<tissue>:macrophage` (all five labels summed), `<tissue>:mac_CD163hi` and
    `<tissue>:mac_CD163lo`, >= 20 cells each, in the genome table's format (stratum = assay, as A284)."""
    import anndata as ad
    a = ad.read_h5ad(IBD_MAC / "pseudobulk.h5ad")
    o = a.obs
    m = (o.tissue == f"{tissue}|healthy").to_numpy()
    X = pd.DataFrame(a.X[m].toarray(), columns=a.var.feature_name.astype(str))
    X = X.loc[:, ~X.columns.duplicated()]
    X = X.reindex(columns=[g for g in genes if g in X.columns]).fillna(0.0)
    meta = pd.DataFrame(dict(donor="ibdverse:" + o.donor_id[m].astype(str).to_numpy(), label=o.cell_type[m].astype(str).to_numpy(),
                             n_cells=o.n_cells[m].to_numpy(), depth=o.depth[m].to_numpy()))
    t = pd.concat([meta, X], axis=1)
    parts = []
    for name, sel in (("macrophage", list(IBD_LABELS)), ("mac_CD163hi", ["mac_CD163hi"]), ("mac_CD163lo", ["mac_CD163lo"])):
        g = t[t.label.isin(sel)].drop(columns="label").groupby("donor", as_index=False).sum()
        g.insert(1, "compartment", f"{tissue}:{name}")
        parts.append(g[g.n_cells >= MIN_CELLS])
    rows = pd.concat(parts, ignore_index=True)
    rows.insert(2, "stratum", "10x 3' v3")
    return rows[rows.donor.isin(donors)].reset_index(drop=True)


def _refine(o, cls):
    cls = cls.copy()
    for k, v in HLCA_SPLIT.items():
        cls[(o.cell_type.astype(str) == k).to_numpy() & (cls == "macrophage").to_numpy()] = v
    return cls


def hlca_split():
    """lung_hlca_genome with the macrophage class split by Census cell_type (other macrophages stay `lung:macrophage`)."""
    import genome as G
    return G.build("lung_hlca_genome", refine=_refine)


if __name__ == "__main__":
    if sys.argv[1] == "ibd":
        import build_a284 as B
        B.build_ibd(sys.argv[2], sys.argv[3], classes=IBD_LABELS,
                    title="IBDverse macrophage labels (A288; built with build_a284.build_ibd)")
