"""GSE148071 (Wu 2021 NSCLC, 42 patients, dense gene x cell tables): marker-rule compartments; histology from histology.csv."""
import re
from geobuild import *
G = "GSE148071"
raw = OUT / G / "raw" / "x"
hist = pd.read_csv(OUT / G / "histology.csv") if (OUT / G / "histology.csv").exists() else pd.DataFrame(columns=["patient", "histology"])
items = []
for f in sorted(raw.glob("*_exp.txt.gz"), key=lambda p: int(re.search(r"_P(\d+)_", p.name).group(1))):
    n = re.search(r"_(P\d+)_", f.name).group(1)
    X, sym, bc = read_dense(f)
    items.append((n, n, X, sym, bc, None)); print(n, X.shape, flush=True)
X, sym, obs = by_sample(items, G, lambda n: "tumour", "dense_10x")
pb = pseudobulk(X, sym, obs, "counts")
age_sex = pd.DataFrame({"patient": [f"P{i}" for i in range(1, 43)]})
donors = age_sex.merge(hist, how="left", on="patient").assign(sample=lambda d: d.patient, treatment="not stated (advanced NSCLC biopsies, GEO)",
        response="", time_fields="none", tissue="tumour (biopsy, tissue site per paper)")
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": ["tumour:" + v for v in MARKER_MAP.values()]})
finish(G, pb, donors, cm)
