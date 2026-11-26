"""GSE141445 (Chen 2021 prostate, 13 tumours pooled in one matrix; patient = barcode suffix -1..-13 = 10x aggregation index).
Raw count matrix used. No author cell labels deposited -> marker rule. Mapping of suffix to the paper's patient ids is not
documented in GEO."""
from geobuild import *
G = "GSE141445"
f = OUT / G / "raw" / "x" / "GSM4203181_data.raw.matrix.txt.gz"
X, sym, bc = read_dense_stream(f)
print(X.shape, X.data[:10], flush=True)
suf = pd.Series(bc).str.extract(r"-(\d+)$", expand=False).astype(int).to_numpy()
comp = np.array([MARKER_MAP.get(x, "drop") for x in marker_comp(X, sym)], dtype=object)
o = pd.DataFrame({"donor": [f"{G}|S{s}" for s in suf], "sample": [f"S{s}" for s in suf],
                  "compartment": ["tumour:" + c if c != "drop" else "drop" for c in comp], "stratum": "10x_cells"})
pb = pseudobulk(X, sym, o, "counts")
donors = pd.DataFrame({"patient": [f"S{i}" for i in range(1, 14)]}).assign(sample=lambda d: d.patient,
    subtype="prostate cancer (12 primary + 1 lymph-node metastasis; which one is not identified in GEO)",
    treatment="not stated", response="", time_fields="none", note="patient = aggregation suffix; mapping to paper ids unknown")
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": ["tumour:" + v for v in MARKER_MAP.values()]})
finish(G, pb, donors, cm)
