"""Per-sample metadata tables (no expression) of the downloaded series, plus whole-matrix scale diagnostics (global min/max/median over all probes,
no gene-level values). Run with python -I."""
import sys, gzip
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import parse_hdr as P
D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_gastric_lung"); OUT = Path(__file__).resolve().parents[1]
FILES = {g: [f"{g}_series_matrix.txt.gz"] for g in ["GSE62254", "GSE15459", "GSE84437", "GSE26253", "GSE30219", "GSE37745", "GSE41271", "GSE81089"]}
def hdr(path):
    p = OUT / "hdr" / (path.name.replace(".txt.gz", ".hdr")); return P.parse(p)
scale = []
for s, fs in FILES.items():
    tabs = []
    for f in fs:
        path = D / s / f; h = hdr(path); t = P.table(h); t.insert(0, "series", s); t["file"] = f; tabs.append(t)
        if s == "GSE81089":  # series matrix carries no expression table (FPKM only in a supplementary file, not downloaded)
            scale.append(dict(series=s, file=f, n_probes=0, n_samples=len(t), log2_like=None)); continue
        ex = pd.read_csv(path, sep="\t", comment="!", index_col=0, quotechar='"', na_values=["null", "NA", ""]).astype(np.float32)
        v = ex.values; scale.append(dict(series=s, file=f, n_probes=ex.shape[0], n_samples=ex.shape[1], global_min=float(np.nanmin(v)), global_median=float(np.nanmedian(v)),
                                         global_max=float(np.nanmax(v)), frac_nan=float(np.isnan(v).mean()), samples_match=list(ex.columns) == list(t.gsm),
                                         log2_like=bool(np.nanmax(v) < 100)))
    m = pd.concat(tabs, ignore_index=True)
    m.to_csv(OUT / "meta" / f"{s}_samples.csv", index=False)
pd.DataFrame(scale).to_csv(OUT / "scale_diagnostics.csv", index=False); print(pd.DataFrame(scale).to_string())
