"""Parse GSE180012 (mouse liver after partial hepatectomy PHx / sham, RNA-seq). Untrusted download: parse only.

    python -I parse_gse180012.py <data dir>
Writes into <data dir>:
  samples.csv             gsm, title, arm (sham|PHx), series (main|day|night), time_h (hours after surgery), rep, platform
  expression_main.parquet main series (48 samples, rows = titles) x mouse genes UPPER-CASED (human-style symbols); values are the
                          file-1 log2(RPKM + 1); symbols damaged by Excel (dates, NaN) dropped, duplicates: highest mean kept
  expression_dn_clock.csv day / night surgery series of file 2 (4 h, 24 h): log2(TPM + 1) of the 14 clock genes (sheet-wise sample names)
"""
import datetime
import gzip
import io
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CLOCK = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"]
d = Path(sys.argv[1])


def sheets(n):
    return pd.ExcelFile(io.BytesIO(gzip.open(d / f"GSE180012_processed_data_file_{n}.xls.gz").read()), engine="xlrd")


# ---- samples from the GEO series matrices -------------------------------------------------------------------------------------------
rows = []
for plat in ("GPL13112", "GPL24247"):
    meta = {}
    with gzip.open(d / f"GSE180012-{plat}_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith(("!Sample_title", "!Sample_geo_accession")):
                meta[line.split("\t")[0]] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "time point" in line:
                meta["tp"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
    for g, t, tp in zip(meta["!Sample_geo_accession"], meta["!Sample_title"], meta["tp"]):
        h = int(re.search(r"(\d+) h after", tp).group(1))
        rows.append(dict(gsm=g, title=t, arm="PHx" if "PHx" in tp else "sham", series="night" if "nighttime" in tp else "main", time_h=h,
                         rep=int(t[-1]), platform=plat))
S = pd.DataFrame(rows)
S.loc[S.series.eq("main") & S.time_h.isin([4, 24]) & S.title.str.startswith("N"), "series"] = "night"
S["title"] = S.title.str.replace("_Rep", "_rep")
S.to_csv(d / "samples.csv", index=False)

# ---- main series (file 1) ------------------------------------------------------------------------------------------------------------
cols, blocks = None, []
for arm in ("sham", "PHx"):
    x = sheets(1).parse(arm)
    ok = x.geneSymbol.map(lambda v: isinstance(v, str))
    x = x[ok].copy()
    x["sym"] = x.geneSymbol.str.upper()
    smp = [c for c in x.columns if re.match(r"\d+[SP]_rep\d$", str(c))]
    x["m"] = x[smp].mean(axis=1)
    x = x.sort_values("m", ascending=False).drop_duplicates("sym").set_index("sym")[smp]
    blocks.append(x.T)
E = pd.concat(blocks, axis=0, join="inner").astype("float32")
E.index.name = "title"
E.to_parquet(d / "expression_main.parquet")
miss = [g for g in CLOCK if g not in E.columns]
print("main", E.shape, "missing clock genes", miss)
assert not miss and set(E.index) <= set(S.title)

# ---- day / night series (file 2) -----------------------------------------------------------------------------------------------------
out = []
x2 = sheets(2)
for sh in x2.sheet_names:
    x = x2.parse(sh)
    sym = next(c for c in ("GeneSymbol", "gene_name") if c in x.columns)
    tpm = [c for c in x.columns if "TPM" in str(c)]
    x = x[x[sym].map(lambda v: isinstance(v, str))].copy()
    x["sym"] = x[sym].str.upper()
    x = x[x.sym.isin(CLOCK)].drop_duplicates("sym").set_index("sym")[tpm]
    x.columns = [re.sub(r"[_ ]?\(?TPM\)?", "", c).replace("_", "-") for c in x.columns]  # D4P-1, N4S-1 ...
    out.append(np.log2(x.T + 1).reindex(columns=CLOCK))
D = pd.concat(out)
D = D[~D.index.duplicated()]
D.index.name = "sample"
D.to_csv(d / "expression_dn_clock.csv")
print("dn", D.shape, list(D.index))
print(S.groupby(["series", "arm", "time_h"]).size().unstack("time_h"))
