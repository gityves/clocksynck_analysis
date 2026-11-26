"""GSE191288 (PTC, 3 patients x left/right lobe + one non-tumour thyroid, 10x h5): marker-rule compartments."""
from geobuild import *
G = "GSE191288"
raw = OUT / G / "raw" / "x"
meta = {"T1L": ("P1", "tumour", "RET/FARP1 fusion", "N0"), "T1R": ("P1", "tumour", "BRAF V600E", "N0"),
        "T2L": ("P2", "tumour", "BRAF V600E", "N1b"), "T2R": ("P2", "tumour", "BRAF V600E", "N1b"),
        "T3L": ("P3", "tumour", "BRAF V600E", "N1a"), "T3R": ("P3", "tumour", "BRAF V600E", "N1a"),
        "NT": ("P4", "adjacent_normal", "wild type", "NA")}
items = []
for f in sorted(raw.glob("*.h5")):
    s = f.stem.split("_", 1)[1]
    X, sym, bc = read_h5_10x(f)
    items.append((s, meta[s][0], X, sym, bc, None)); print(s, X.shape, flush=True)
X, sym, obs = by_sample(items, G, lambda n: meta[n][1], "10x_cells")
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame([dict(patient=v[0], sample=k, tissue=v[1], subtype="papillary thyroid carcinoma" if v[1] == "tumour" else "non-tumour thyroid",
         genotype=v[2], lymph_node_stage=v[3], treatment="not stated", response="", time_fields="none") for k, v in meta.items()])
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": list(MARKER_MAP.values())})
finish(G, pb, donors, cm)
