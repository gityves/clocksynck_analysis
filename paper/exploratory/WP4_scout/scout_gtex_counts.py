"""WP4 scout (A123, A131): aggregate counts of GTEx v10 fast-death donors in the 08-18 h window, per tissue, and SMCENTER codes.

Reads only `obs` of the labelled h5ad (no expression). Controlled access (time of death): aggregates only, no per-donor rows are
printed or written. Writes output/gtex_window_counts.csv (per tissue) and output/gtex_center_codes.csv (per code).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import anndata as ad
import numpy as np
import pandas as pd

from common.benchmark.data import in_window
from common.gtex import FAST, LABELLED

OUT = Path(__file__).parent / "output"
a = ad.read_h5ad(LABELLED, backed="r")
o = a.obs[["tissue", "subject", "DTHHRDY", "tod", "center"]].copy()
a.file.close()
o = o[~o.tissue.astype(str).str.startswith("Cells -")]
fast = o[o.DTHHRDY.isin(list(FAST)) & o.tod.notna()].drop_duplicates(["subject", "tissue"])
fast = fast.assign(win=in_window(fast.tod.astype(float) % 24))
don = fast.drop_duplicates("subject")
print("fast donors with tod:", don.subject.nunique(), "| in 08-18 h:", int(don.win.sum()))
t = fast.groupby("tissue").agg(n_fast=("subject", "nunique"), n_window=("win", "sum")).astype(int)
t["ge20"] = t.n_window >= 20
t.sort_values("n_window", ascending=False).to_csv(OUT / "gtex_window_counts.csv")
print("tissues:", len(t), "| with >=20 window donors:", int(t.ge20.sum()), "| >=20 fast:", int((t.n_fast >= 20).sum()))
print("benchmark 8 tissues:")
B = ("Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Nerve - Tibial", "Skin - Sun Exposed (Lower leg)",
     "Skin - Not Sun Exposed (Suprapubic)", "Adipose - Subcutaneous", "Lung")
print(t.reindex(B).to_string())
# centres (all labelled samples, not only fast)
c = o.drop_duplicates("subject").center.astype(str).value_counts()
c.rename_axis("center").rename("n_donors_all_hardy").to_frame().to_csv(OUT / "gtex_center_codes.csv")
print("center codes (donors, first sample):", c.to_dict())
f = fast.drop_duplicates("subject")
print("fast donors with tod per center:", f.center.astype(str).value_counts().to_dict(),
      "| in window:", f[f.win].center.astype(str).value_counts().to_dict())
