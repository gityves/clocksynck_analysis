"""Parse GSE13093 (Vollmers et al. 2009, mouse liver, GPL1261) into clock-gene tables. Untrusted download: parse only.

    python -I parse_gse13093.py <data dir>
Writes <dir>/samples.csv (gsm, group, ct) and <dir>/expression_clock.csv (samples x GTEx clock genes, log2; per gene
the probe set with the highest mean over all samples).
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GENES = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"]
ALIAS = {"CIART": ["Ciart", "Gm129"]}

d = Path(sys.argv[1])
meta, rows = {}, []
with gzip.open(d / "GSE13093_series_matrix.txt.gz", "rt") as f:
    for line in f:
        if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
            meta[line.split("\t")[0]] = [x.strip().strip('"').strip() for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith('"ID_REF"'):
            cols = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            break
    for line in f:
        if line.startswith("!series_matrix_table_end"):
            break
        rows.append(line.rstrip("\n").split("\t"))
X = pd.DataFrame(rows, columns=cols).set_index("ID_REF")
X.index = X.index.str.strip('"')
X = X.apply(pd.to_numeric, errors="coerce")
titles = meta["!Sample_title"]
s = pd.DataFrame({"gsm": meta["!Sample_geo_accession"], "title": titles})
s["group"] = [re.sub(r"_CT\d+$", "", t).strip() for t in titles]
s["ct"] = [float(m.group(1)) if (m := re.search(r"CT(\d+)$", t)) else np.nan for t in titles]
s.to_csv(d / "samples.csv", index=False)

with gzip.open(d / "GPL1261.annot.gz", "rt", errors="replace") as f:
    lines = f.read().split("\n")
b = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
e = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
hdr = lines[b + 1].split("\t")
ann = pd.DataFrame([l.split("\t") for l in lines[b + 2:e]], columns=hdr).set_index("ID")
sym = ann["Gene symbol"].dropna()
print("value range:", float(np.nanmin(X.values)), float(np.nanmax(X.values)))
out = {}
for g in GENES:
    names = ALIAS.get(g, [g.capitalize()])
    probes = [p for p in sym.index if sym[p] in names and p in X.index]
    if not probes:
        print("missing", g)
        continue
    best = X.loc[probes].mean(1).idxmax()
    out[g] = X.loc[best]
E = pd.DataFrame(out)
if np.nanmax(X.values) > 100:  # MAS5 linear scale
    E = np.log2(E.clip(lower=1))
E.index.name = "gsm"
E.to_csv(d / "expression_clock.csv")
print(E.shape, s.group.value_counts().to_dict())
