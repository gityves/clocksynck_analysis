"""GSE212527 (sarcoma UPS / leiomyosarcoma, 6 patients, 7 samples, 10x): author Cell.Annotations (annotated barcodes of the
raw matrices only)."""
from geobuild import *
G = "GSE212527"
raw = OUT / G / "raw"
a = pd.read_csv(raw / "GSE212527_Cell.Annotations.txt.gz", sep="\t")
M = {"Malignant fibroblasts": "malignant", "Malignant epithelial-like cells": "malignant", "Cancer-associated fibroblasts": "fibroblast",
     "Endothelial cells": "endothelial", "Monocytes and macrophages": "myeloid", "Dendritic cells": "myeloid", "CD8 T cells": "T_NK",
     "CD4 T cells": "T_NK", "NK cells": "T_NK", "B cells": "B_plasma", "Mast cells": "drop"}
sub = {"SRC164": "UPS", "SRC147": "leiomyosarcoma", "SRC141": "UPS", "SRC171": "UPS", "SRC241": "leiomyosarcoma",
       "SRC242_biceps": "leiomyosarcoma", "SRC242_pectoralis": "leiomyosarcoma"}
items = []
for s in sub:
    pre = next((raw / "x").glob(f"GSM*_{s}_barcodes.tsv.gz")).name.replace("barcodes.tsv.gz", "")
    X, sym, bc = read_10x(raw / "x" / pre)
    t = a[a.Sample == s].set_index("Cell_Barcode").Cell_Type
    assert t.index.isin(bc).all(), s
    pos = pd.Series(np.arange(len(bc)), index=bc)[t.index].to_numpy()
    items.append((s, s.split("_")[0], X[pos], sym, t.index.to_numpy(), np.array([M[x] for x in t], dtype=object)))
    print(s, X.shape, len(pos), flush=True)
X, sym, obs = by_sample(items, G, lambda n: "tumour", "10x_cells")
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame([dict(patient=k.split("_")[0], sample=k, subtype=v, treatment="not stated", response="", time_fields="none",
        note="biceps / pectoralis = two samples of SRC242 (inferred from names)" if "242" in k else "") for k, v in sub.items()])
cm = pd.crosstab(a.Cell_Type.rename("a"), a.Cell_Type.map(M).rename("b")).stack().rename("n").reset_index()
cm = cm[cm.n > 0].rename(columns={"a": "author_cell_type", "b": "compartment"})
finish(G, pb, donors, cm)
