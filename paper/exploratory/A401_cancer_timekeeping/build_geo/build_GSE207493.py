"""GSE207493 (ccRCC, 19 unsorted tumour homogenates, 10x): marker-rule compartments (no author labels deposited)."""
import numpy as np, pandas as pd, scipy.sparse as sp
from geobuild import *
G = "GSE207493"
raw = OUT / G / "raw" / "x"
Xs, obs, syms = [], [], None
for d in sorted(p for p in raw.iterdir() if p.is_dir()):
    X, sym, bc = read_10x(d)
    comp = marker_comp(X, sym)
    c = np.array([MARKER_MAP.get(x, "drop") for x in comp])
    o = pd.DataFrame({"donor": f"{G}|{d.name}", "sample": d.name, "compartment": ["tumour:" + x if x != "drop" else x for x in c],
                      "stratum": "10x_v3_cells"})
    if syms is None: syms = sym
    assert (sym == syms).all()
    Xs.append(X); obs.append(o); print(d.name, X.shape, flush=True)
X = sp.vstack(Xs).tocsr(); obs = pd.concat(obs, ignore_index=True)
pb = pseudobulk(X, syms, obs, "counts")
donors = pd.DataFrame({"patient": sorted(obs["sample"].unique())}).assign(sample=lambda d: d.patient, subtype="ccRCC",
        treatment="not stated (GEO: total tumour homogenate)", response="", time_fields="none")
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": ["tumour:" + v for v in MARKER_MAP.values()]})
finish(G, pb, donors, cm)
