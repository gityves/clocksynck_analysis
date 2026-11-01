"""GSE154200 (lesional / non-lesional skin, Munich): persons with both biopsies per diagnosis, from sample metadata only.
No expression value is read. Run: python -I meta_gse154200.py"""
from pathlib import Path
import pandas as pd

SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
OUT = Path(__file__).resolve().parents[1] / "output"
rows, cur = [], None
for l in open(SCCC / "skin_gse154200/GSE154200_gsm_brief.txt", encoding="utf-8", errors="replace"):
    l = l.rstrip("\n")
    if l.startswith("^SAMPLE"):
        cur = {"gsm": l.split("=")[1].strip()}
        rows.append(cur)
    elif l.startswith("!Sample_characteristics_ch1"):
        k, _, v = l.split("=", 1)[1].partition(":")
        cur[k.strip()] = v.strip()
    elif l.startswith("!Sample_source_name_ch1"):
        cur["source"] = l.split("=", 1)[1].strip()
    elif l.startswith("!Sample_title"):
        cur["title"] = l.split("=", 1)[1].strip()
m = pd.DataFrame(rows)
print(m.shape, m.columns.tolist(), m.source.value_counts().to_dict(), m.diagnosis.value_counts().to_dict())
m["state"] = m.source.map(lambda s: "lesional" if s.startswith("lesional") else "nonlesional" if s.startswith("non") else s)
c = m.groupby(["diagnosis", "patient id", "state"]).size().unstack(fill_value=0)
for k in ("lesional", "nonlesional"):
    c[k] = c.get(k, 0)
c["both"] = (c.lesional > 0) & (c.nonlesional > 0)
c["one_one"] = (c.lesional == 1) & (c.nonlesional == 1)
t = c.groupby("diagnosis").agg(patients=("both", "size"), both=("both", "sum"), exactly_1_plus_1=("one_one", "sum"))
t.loc["pooled"] = t.sum()
t["n_samples"] = [int((m.diagnosis == d).sum()) for d in t.index[:-1]] + [len(m)]
t["patients_in_>1_diagnosis"] = int((m.groupby("patient id").diagnosis.nunique() > 1).sum())
t.to_csv(OUT / "gse154200_counts.csv")
m.drop(columns=["title"]).to_csv(OUT / "gse154200_samples.csv", index=False)
print(t.to_string()); print(m.source.value_counts().to_dict())
