"""GSE197461 (cervical cancer, 3 SCC + 5 ADC, 10x 5'): marker-rule compartments."""
from geobuild import *
G = "GSE197461"
raw = OUT / G / "raw" / "x"
hpv = {"SCC_1": "HPV58", "SCC_2": "HPV16", "SCC_3": "HPV16", "ADC_1": "HPV16", "ADC_2": "HPV18", "ADC_3": "HPV45", "ADC_4": "HPV-negative", "ADC_5": "HPV-negative"}
items = []
for s in hpv:
    pre = next(raw.glob(f"GSM*_{s}_barcodes.tsv.gz")).name.replace("barcodes.tsv.gz", "")
    X, sym, bc = read_10x(raw / pre)
    items.append((s, s, X, sym, bc, None)); print(s, X.shape, flush=True)
X, sym, obs = by_sample(items, G, lambda n: "tumour", "10x_5prime")
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame([dict(patient=k, sample=k, subtype="cervical " + k[:3], hpv=v, treatment="not stated", response="", time_fields="none") for k, v in hpv.items()])
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": ["tumour:" + v for v in MARKER_MAP.values()]})
finish(G, pb, donors, cm)
