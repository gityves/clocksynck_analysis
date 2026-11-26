"""GSE163203 (PTC +/- Hashimoto, Microwell-seq, 8 patients, 2 adjacent normal): marker-rule compartments."""
import re
from geobuild import *
G = "GSE163203"
raw = OUT / G / "raw" / "x"
items, meta = [], {}
for f in sorted(raw.glob("*_dge.txt.gz")):
    m = re.match(r"GSM\d+_(Adj_)?PTC(with|without)HT_(\d+)(?:_\d)?_dge", f.name)
    adj, ht, n = m.group(1), m.group(2), m.group(3)
    sample = ("Adj_" if adj else "") + f"PTC_{ht}HT_{n}"
    X, sym, bc = read_dense(f)
    meta[sample] = (f"P{n}", "adjacent_normal" if adj else "tumour", ht + " Hashimoto thyroiditis")
    items.append((sample, f"P{n}", X, sym, bc, None)); print(f.name, X.shape, flush=True)
# patient 4 has two files (4_1, 4_2) of one GEO sample: kept under one sample name, concatenated
X, sym, obs = by_sample(items, G, lambda n: meta[n][1], "Microwell-seq")
obs["sample"] = obs["sample"]
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame([dict(patient=v[0], sample=k, tissue=v[1], subtype="papillary thyroid carcinoma" if v[1] == "tumour" else "normal adjacent",
         hashimoto=v[2], treatment="not stated", response="", time_fields="none") for k, v in meta.items()])
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": list(MARKER_MAP.values())})
finish(G, pb, donors, cm)
