"""A196 data preparation (blind to the clock statistic): GSE93204 (NeoPalAna) series matrix and GPL6480 annotation parsed with the
A168 prep functions (same conventions as GSE87411: log2 ratios to the universal reference, probes averaged per symbol, genes missing in
> 20 % of samples dropped, A97 scores over the series) into ../output/neopal_{meta,expr}.parquet. No coherence is computed here.
Run: python -I prep_a196.py   (reads downloaded text files only)"""
import re
import sys
from pathlib import Path

A168 = Path(__file__).resolve().parents[2] / "A168_endocrine" / "scripts"
sys.path.insert(0, str(A168))
import pandas as pd  # noqa: E402

import a168_lib as L  # noqa: E402
import prep_a168 as P  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(exist_ok=True)
L.OUT = OUT  # `finish` writes to L.OUT
P.PLATFORM["GSE93204"] = "GPL6480"


def neopal():
    gse = "GSE93204"; h = P.header(P.D / gse / f"{gse}_series_matrix.txt.gz")
    t = h["!Sample_title"][0]; gsm = h["!Sample_geo_accession"][0]
    rows = []
    for g, x in zip(gsm, t):
        p, tp = re.fullmatch(r"Breast Cancer Array_Patient ID_(\w+?)_(BL|C1D1|C1D15|C1D15_Redo|Surg|S)", x).groups()
        rows.append((g, p, "NeoPalAna", {"C1D15_Redo": "C1D15", "S": "Surg"}.get(tp, tp), x))  # one redo, one "S" title
    m = pd.DataFrame(rows, columns=["gsm", "patient", "arm", "tp", "title"]).set_index("gsm")
    assert not m.duplicated(["patient", "tp"]).any()
    P.finish(gse, "neopal", m, P.expression(gse))
    pd.DataFrame(P.LOG, columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)
    print(m.tp.value_counts().to_string())


if __name__ == "__main__":
    neopal()
