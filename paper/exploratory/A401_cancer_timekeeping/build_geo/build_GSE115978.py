"""GSE115978 (Jerby-Arnon 2018 melanoma, Smart-seq2, 32 samples / 31 patients, 'counts.csv'): author cell types.
Treatment group (naive / post-treatment) from the GEO cell annotation file; response labels are NOT in GEO (paper Table S1A)."""
import re
from geobuild import *
G = "GSE115978"
raw = OUT / G / "raw"
a = pd.read_csv(raw / "GSE115978_cell.annotations.csv.gz").set_index("cells")
C = pd.read_csv(raw / "GSE115978_counts.csv.gz", index_col=0)
print(C.shape, flush=True)
a = a.loc[C.columns]
vals = C.to_numpy(dtype=np.float32)
print("non-integer fraction of nonzeros:", float((np.abs(vals[vals > 0] - np.round(vals[vals > 0])) > 1e-6).mean()))
X = sp.csr_matrix(vals.T); sym = np.array(C.index.astype(str)); del vals
M = {"Mal": "malignant", "CAF": "fibroblast", "Endo.": "endothelial", "Macrophage": "myeloid", "T.CD8": "T_NK", "T.CD4": "T_NK",
     "T.cell": "T_NK", "NK": "T_NK", "B.cell": "B_plasma", "?": "drop"}
pat = a.samples.str.replace(r"(pa|pb)$", "", regex=True).str.replace(r"\.\d+$", "", regex=True)
o = pd.DataFrame({"donor": G + "|" + pat.to_numpy(), "sample": a.samples.to_numpy(),
                  "compartment": ["tumour:" + M[t] if M[t] != "drop" else "drop" for t in a["cell.types"]], "stratum": "Smart-seq2"})
pb = pseudobulk(X, sym, o, "counts (author counts.csv; non-integer if RSEM-estimated)")
donors = (a.groupby("samples").agg(treatment_group=("treatment.group", "first"), cohort=("Cohort", "first"), cells=("samples", "size")).reset_index()
          .rename(columns={"samples": "sample"}))
donors.insert(0, "patient", donors["sample"].str.replace(r"(pa|pb)$", "", regex=True).str.replace(r"\.\d+$", "", regex=True))
donors["subtype"] = "melanoma"
donors["treatment"] = donors.treatment_group.map({"treatment.naive": "naive", "post.treatment": "post-treatment (type not in GEO; paper: ICI)"})
donors["treatment_source"] = "GSE115978_cell.annotations.csv.gz column treatment.group"
donors["response"] = ""; donors["response_source"] = "NOT AVAILABLE: paper Table S1A (Cell 2018), not retrievable from this machine"
donors["time_fields"] = "none"
cm = pd.crosstab(a["cell.types"].rename("a"), pd.Series([M[t] for t in a["cell.types"]], index=a.index, name="b")).stack().rename("n").reset_index()
cm = cm[cm.n > 0].rename(columns={"a": "author_cell_type", "b": "compartment"})
finish(G, pb, donors, cm)
