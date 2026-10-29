"""A161 step 1-2 (before any clock readout): phase labels, endometriosis status, patient overlap by sample id.
python -I s01_metadata.py"""
import sys, re
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import numpy as np, pandas as pd
import a161_lib as L

m = L.header(L.EN / "GSE51981_series_matrix.txt.gz")
print(m.columns.tolist())
m["gsm"] = m.geo_accession
t = m["tissue"].str.lower()
m["phase_label"] = np.select([t.str.contains("early secretory"), t.str.contains("mid-secretory|mid secretory"),
                              t.str.contains("late secretory|late-secretory"), t.str.contains("proliferative")],
                             ["early_secretory", "mid_secretory", "late_secretory", "proliferative"], "unknown")
m["phase"] = m.phase_label.map({"proliferative": "P", "early_secretory": "S", "mid_secretory": "S", "late_secretory": "S"}).fillna("NA")
m["endo"] = (m["endometriosis/no endometriosis"].str.lower() == "endometriosis").astype(int)
m["pathology"] = m.source_name_ch1
keep = ["gsm", "title", "tissue", "phase_label", "phase", "endometriosis/no endometriosis", "endo", "endometriosis severity", "source_name_ch1"]
m[keep].to_csv(L.OUT / "gse51981_samples.csv", index=False)
print(m.phase_label.value_counts().to_dict())
print(pd.crosstab(m.phase_label, m.endo))
print(m.source_name_ch1.value_counts().to_string())
print("duplicate titles:", m.title.duplicated().sum())

# overlap of sample ids
g6 = L.header(L.EN / "GSE6364_series_matrix.txt.gz"); g6["gsm"] = g6.geo_accession
g4 = L.soft_samples(L.SC / "GSE4888_family.soft.gz")
def pid(s):  # patient id: last token of title (PE_D_26A -> 26A; 'endometrium 598' -> 598)
    return re.split(r"[_ ]", s)[-1].upper()
a = pd.DataFrame({"gsm": m.gsm, "id": m.title.map(pid)}); a["series"] = "GSE51981"
b = pd.DataFrame({"gsm": g6.gsm, "id": g6.title.map(pid)}); b["series"] = "GSE6364"
c = pd.DataFrame({"gsm": g4.index, "id": g4.title.map(pid).values}); c["series"] = "GSE4888"
allid = pd.concat([a, b, c])
tab = allid.groupby("id").series.agg(lambda s: ",".join(sorted(s))).rename("series")
ov = pd.Series({"GSE51981&GSE6364": len(set(a.id) & set(b.id)), "GSE51981&GSE4888": len(set(a.id) & set(c.id)),
                "GSE6364&GSE4888": len(set(b.id) & set(c.id)), "n_51981": a.id.nunique(), "n_6364": b.id.nunique(), "n_4888": c.id.nunique()})
print(ov)
ov.to_csv(L.OUT / "patient_overlap_summary.csv", header=["n"])
allid.to_csv(L.OUT / "patient_ids_all_series.csv", index=False)
print(sorted(set(a.id) & set(b.id)))
print(sorted(set(b.id)))
g6[["gsm", "title"]].assign(char=g6["char"] if "char" in g6 else "").to_csv(L.OUT / "gse6364_samples.csv", index=False)
