"""GSE155698 (Peng 2019-type PDAC, 10x, 16 treatment-naive tumours + 3 adjacent normal; PBMC samples skipped): marker rule."""
import re
from geobuild import *
G = "GSE155698"
raw = OUT / G / "raw" / "x"
items, meta = [], {}
for d in sorted(p for p in raw.iterdir() if p.is_dir()):
    n = d.name
    adj = n.startswith("AdjNorm")
    pat = re.sub(r"[AB]$", "", n.split("TISSUE_")[1])
    pat = f"AdjNorm{pat}" if adj else f"P{pat}"  # adjacent-normal samples are not linked to tumour patients in GEO
    if (d / "filtered_feature_bc_matrix").is_dir():
        X, sym, bc = read_10x(d / "filtered_feature_bc_matrix")
    elif (d / "filtered_gene_bc_matrices").is_dir():  # TISSUE_16: v2 layout, uncompressed
        X, sym, bc = read_10x(d / "filtered_gene_bc_matrices", gz=False)
    else:  # TISSUE_14: h5 only
        X, sym, bc = read_h5_10x(d / "filtered_feature_bc_matrix.h5")
    meta[n] = (pat, "adjacent_normal" if adj else "tumour")
    items.append((n, pat, X, sym, bc, None)); print(n, X.shape, flush=True)
X, sym, obs = by_sample(items, G, lambda n: meta[n][1], lambda n: "10x_HiSeq4000(GPL20301)" if n in ("PDAC_TISSUE_13", "PDAC_TISSUE_16") else "10x_NovaSeq6000(GPL24676)")
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame([dict(patient=v[0], sample=k, tissue=v[1], subtype="PDAC", treatment="Treatment naive (GEO)", response="",
        time_fields="none", note="11A/11B are two samples of patient 11" if v[0] == "P11" else
        ("adjacent normal not linked to a tumour patient in GEO" if v[1] == "adjacent_normal" else "")) for k, v in meta.items()])
cm = pd.DataFrame({"source": "cancer_common marker rule", "label": list(MARKER_MAP), "compartment": list(MARKER_MAP.values())})
finish(G, pb, donors, cm)
