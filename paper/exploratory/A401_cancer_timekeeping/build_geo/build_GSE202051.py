"""GSE202051 (Hwang 2022 PDAC snRNA-seq, 10x): untreated tumours only (U1-U18); author labels (Level 1/2 Annotation)."""
import anndata as ad, h5py
from geobuild import *
G = "GSE202051"
f = OUT / G / "raw" / "GSE202051_totaldata-final-toshare.h5ad"
a = ad.read_h5ad(f, backed="r")
obs = a.obs[["sampleid", "pid", "treatment_status", "response", "Level 1 Annotation", "Level 2 Annotation"]].copy()
keep = (obs.treatment_status == "Untreated").to_numpy()
with h5py.File(f, "r") as h:
    X = ad.io.read_elem(h["layers/counts"]).tocsr()
sym = np.array(a.var_names)
X = X[np.flatnonzero(keep)]
obs = obs[keep]
l1, l2 = obs["Level 1 Annotation"].astype(str), obs["Level 2 Annotation"].astype(str)
comp = pd.Series("drop", index=obs.index)
comp[l1 == "Epithelial (malignant)"] = "tumour:malignant"
comp[l1 == "Cancer-associated fibroblast"] = "tumour:fibroblast"
comp[l1 == "Endothelial"] = "tumour:endothelial"
comp[l1.isin(["Pericyte", "Vascular smooth muscle"])] = "tumour:mural"
comp[l1 == "Myeloid"] = "tumour:myeloid"
comp[(l1 == "Lymphoid") & ~l2.isin(["B", "Plasma"])] = "tumour:T_NK"
comp[(l1 == "Lymphoid") & l2.isin(["B", "Plasma"])] = "tumour:B_plasma"
o = pd.DataFrame({"donor": G + "|" + obs.pid.astype(str).to_numpy(), "sample": obs.sampleid.astype(str).to_numpy(),
                  "compartment": comp.to_numpy(), "stratum": "10x_snRNA"})
pb = pseudobulk(X, sym, o, "counts")
donors = obs.groupby("pid", observed=True).agg(sample=("sampleid", lambda s: ";".join(sorted(set(map(str, s))))), cells=("pid", "size")).reset_index()
donors = donors.rename(columns={"pid": "patient"}).assign(subtype="PDAC", treatment="Untreated (author treatment_status)",
        response="", time_fields="none", source="author obs of GSE202051 h5ad")
cm = pd.crosstab(l1.to_numpy(), comp.to_numpy()).stack().reset_index()
cm.columns = ["author_level1_annotation", "compartment", "n_cells"]
cm = cm[cm.n_cells > 0]
cm["note"] = "Level 2 splits Lymphoid into B/Plasma (tumour:B_plasma) vs other (T_NK); other epithelial, endocrine, Schwann, adipocyte, neurons dropped"
finish(G, pb, donors, cm)
