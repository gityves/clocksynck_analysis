"""C4 scout, Wu 2021 and Pal 2021 breast: malignant cells per patient, UMI depth, clock-gene detection (counts only).
Wu per-cell table prepared by a chunked read of the mtx (scratch pickle); Pal from A170 units (pseudobulk depth, no per-cell file)."""
import sys, pandas as pd, numpy as np
from pathlib import Path
OUT = Path(__file__).resolve().parent
SC = "/tmp/claude-1000/-home-maxine-Documents-paychere-scCircadianOverview/255baa01-84eb-4f87-a164-c3c898245f47/scratchpad/wu_cells.pkl"
G = "ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART".split()
md = pd.read_pickle(SC)
rec = pd.read_csv(sys.argv[1]).set_index("orig.ident") if len(sys.argv) > 1 else None
m = md[md.celltype_major == "Cancer Epithelial"]
r = {"cells": len(m), "patients": m.orig_ident.nunique() if "orig_ident" in m else m["orig.ident"].nunique(), "median_umi": m.umi.median(),
     "median_nCount_RNA_author": m.nCount_RNA.median()}
for g in G: r["det_" + g] = (m["c_" + g] >= 1).mean()
r["det_mean14"] = np.mean([r["det_" + g] for g in G]); r["frac_any_clock"] = (m[["c_" + g for g in G]].sum(1) >= 1).mean()
r["median_clock_umi"] = m[["c_" + g for g in G]].sum(1).median()
pd.DataFrame([r]).to_csv(OUT / "wu_cell_level_detection.csv", index=False)
p = m.groupby("orig.ident").agg(n_malignant=("umi", "size"), median_umi=("umi", "median"), subtype_author=("subtype", "first"))
if rec is not None: p = p.join(rec[[c for c in rec.columns if c in ("final_ER", "final_label", "treatment")]])
allp = md.groupby("orig.ident").size().rename("n_cells_all")
p = allp.to_frame().join(p).fillna({"n_malignant": 0}); p.to_csv(OUT / "wu_per_patient.csv")
print(pd.DataFrame([r]).round(3).T); print(p.to_string())
print("ge300", int((p.n_malignant >= 300).sum()), "ge500", int((p.n_malignant >= 500).sum()), "ge50", int((p.n_malignant >= 50).sum()))
