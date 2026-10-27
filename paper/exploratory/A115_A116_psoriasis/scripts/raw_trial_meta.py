"""Tidy per-sample metadata (patient, tissue, timepoint, arm) of the candidate trial series from the raw series matrices on
disk (metadata only, no expression). Output: output/trial_meta_<series>.csv and output/a116_trial_counts.csv (patients per
arm with lesional AND uninvolved biopsies at baseline and at each follow-up). No offset or response statistic is computed.

    python -I raw_trial_meta.py
"""
import gzip, re, sys, collections
from pathlib import Path
import pandas as pd
ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/skin_lesional")
OUT = Path(__file__).resolve().parents[1] / "output"


def read(s):
    f = next((ROOT / s).glob("*series_matrix.txt.gz"))
    F = collections.defaultdict(list)
    with gzip.open(f, "rt", errors="replace") as h:
        for line in h:
            if line.startswith("!series_matrix_table_begin"): break
            if line.startswith("!Sample_"):
                p = line.rstrip("\n").split("\t"); F[p[0]].append([x.strip('"') for x in p[1:]])
    n = len(F["!Sample_title"][0])
    ch = {}
    for col in F["!Sample_characteristics_ch1"]:
        k = col[0].split(":")[0].strip() if col and ":" in col[0] else None
        if k: ch[k] = [x.split(":", 1)[1].strip() if ":" in x else "" for x in col]
    return F["!Sample_title"][0], F["!Sample_geo_accession"][0], ch, n


def tidy(s):
    title, gsm, ch, n = read(s)
    d = pd.DataFrame({"gsm": gsm, "title": title})
    for k, v in ch.items(): d[k] = v
    return d


def norm(s, d):
    t = d.title
    if s == "GSE283313":
        m = t.str.extract(r"Pt (\d+) W(\d+) (NL|L) (\w+)"); return pd.DataFrame(dict(patient=m[0], week=m[1].astype(float), tissue=m[2].map({"NL": "NL", "L": "LS"}), arm=m[3]))
    if s == "GSE69967":
        m = t.str.extract(r"(\d+):(\w+),(NL|LS),Day (\d+)"); return pd.DataFrame(dict(patient=m[0], week=m[3].astype(float) / 7, tissue=m[2], arm=m[1]))
    if s == "GSE252029":
        return pd.DataFrame(dict(patient=d["subject id"], week=d["time point"].str.extract(r"(\d+)")[0].astype(float), tissue=d["skin"], arm=d["treatment"]))
    if s == "GSE136757":
        return pd.DataFrame(dict(patient=d["subjid"], week=d["time"].str.extract(r"(\d+)")[0].astype(float) / 7, tissue=d["tissue"].map({"lesional skin": "LS", "non-lesional skin": "NL"}), arm=d["treatment"]))
    if s == "GSE193128":
        return pd.DataFrame(dict(patient=d["patient_no"], week=d["week"].str.extract(r"(\d+)")[0].astype(float), tissue=d["disease state (nl/l)"].map({"NL": "NL", "L": "LS"}), arm="?"))
    if s == "GSE137218":
        return pd.DataFrame(dict(patient=t.str.extract(r"patient(\d+)")[0], week=d["day"].astype(float) / 7, tissue=d["disease state (ls/nl)"], arm="secukinumab"))
    if s == "GSE197056":
        return pd.DataFrame(dict(patient=d["individual"], week=d["timepoint"].str.extract(r"(\d+)")[0].astype(float).fillna(0), tissue=d["tissue"], arm=d["treatment"]))
    if s == "GSE268958":
        return pd.DataFrame(dict(patient=d["subjectid"], week=d["timepoint"].str.extract(r"(\d+)")[0].astype(float).fillna(0), tissue=d["tissue"].str.lower().map(lambda x: "NL" if "non" in x else ("LS" if "les" in x else x)), arm="brodalumab"))
    if s == "GSE85034":
        return pd.DataFrame(dict(patient=t.str.extract(r"Subject (\d+)")[0], week=d["timepoint"], tissue=d["timepoint"], arm=d["treatment"]))
    return None


rows = []
for s in ["GSE283313", "GSE69967", "GSE252029", "GSE136757", "GSE193128", "GSE137218", "GSE197056", "GSE268958", "GSE85034"]:
    d = tidy(s); n = norm(s, d)
    d.to_csv(OUT / f"trial_meta_raw_{s}.csv", index=False)
    if s == "GSE85034":
        print(s, n.groupby(["tissue", "arm"]).size()); continue
    n = n.dropna(subset=["patient"]); n["week"] = n.week.round(1)
    # baseline = week 0 (day 0); follow-up visits = other weeks
    for arm, g in n.groupby("arm"):
        base = set(g[(g.week == 0) & (g.tissue == "LS")].patient) & set(g[(g.week == 0) & (g.tissue == "NL")].patient)
        for wk in sorted(set(g.week.dropna()) - {0}):
            fu_ls = set(g[(g.week == wk) & (g.tissue == "LS")].patient); fu_nl = set(g[(g.week == wk) & (g.tissue == "NL")].patient)
            rows.append(dict(series=s, arm=arm, followup_week=wk, n_baseline_LS_and_NL=len(base), n_LS_both_visits=len(base & fu_ls),
                             n_LS_and_NL_both_visits=len(base & fu_ls & fu_nl), n_arm_patients=g.patient.nunique()))
r = pd.DataFrame(rows); r.to_csv(OUT / "a116_trial_counts.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500); print(r.to_string())
