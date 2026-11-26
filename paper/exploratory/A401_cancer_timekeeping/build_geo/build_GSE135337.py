"""GSE135337 (bladder cancer, 7 tumours + one paracancerous sample BCN, 10x gene x cell tables): marker-rule compartments."""
import re
from geobuild import *
G = "GSE135337"
raw = OUT / G / "raw" / "x"
age = {"BC1": 67, "BC2": 67, "BC3": 38, "BC4": 80, "BC5": 81, "BC6": 58, "BC7": 66, "BCN": 80}
stage = {"BC1": "Ta", "BC2": "T1", "BC3": "T1", "BC4": "T2", "BC5": "T3", "BC6": "pTa", "BC7": "pTa", "BCN": "pT2 (as listed; paracancerous)"}
items = []
for f in sorted(raw.glob("*_gene_cell_exprs_table.*gz")):
    s = re.match(r"GSM\d+_(BC\w+?)_gene", f.name).group(1)
    X, sym, bc = read_dense(f, sym_col="Symbol", skip_cols=2)
    items.append((s, s, X, sym, bc, None)); print(s, X.shape, flush=True)
X, sym, obs = by_sample(items, G, lambda n: "adjacent_normal" if n == "BCN" else "tumour", "10x_cells")
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame([dict(patient=k, sample=k, tissue="adjacent_normal (paracancerous)" if k == "BCN" else "tumour",
    subtype="bladder cancer (UC)", age=age[k], tumour_stage=stage[k],
    note=("BCN is not assigned to a tumour patient in GEO; age 80 and pT2 equal BC4, match unconfirmed" if k == "BCN" else ""),
    treatment="not stated", response="", time_fields="none") for k in age])
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": list(MARKER_MAP.values())})
finish(G, pb, donors, cm)
