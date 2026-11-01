"""C4 scout: per patient, malignant cells (copy-number call) in the legacy CoPhaser tables and how many are cycling (cycling_logit > 0);
cell counts only, no clock-gene column is read. Run with python -I."""
import pandas as pd, numpy as np
from pathlib import Path
OUT = Path(__file__).resolve().parent
D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian"); CC = Path("/home/maxine/Documents/paychere/scCircadianCancer/outputs/states/cc")
rows = []
for f in sorted(CC.glob("*.csv.gz")):
    st = f.name.split(".")[0]; cf = D / "cnv" / f"{st}.csv.gz"
    if not cf.exists(): continue
    t = pd.read_csv(f, usecols=["Unnamed: 0", "patient", "kind", "total_umi", "theta", "cycling_logit"], index_col=0)
    t = t[t.kind == "tumour"]
    c = pd.read_csv(cf, usecols=["sample", "barcode", "cnv_call"]); c.index = c["sample"].astype(str) + "_" + c["barcode"].astype(str)
    m = t.join(c.cnv_call, how="inner"); m = m[m.cnv_call == "malignant"]
    for p, g in m.groupby("patient"):
        rows.append(dict(study=st, patient=p, n_mal=len(g), n_cycling=int((g.cycling_logit > 0).sum()), median_umi=g.total_umi.median()))
R = pd.DataFrame(rows); R["frac_cycling"] = R.n_cycling / R.n_mal
R.to_csv(OUT / "per_patient_malignant_cycling_counts.csv", index=False)
S = R.groupby("study").agg(n_pat=("patient", "size"), n_ge300=("n_mal", lambda x: int((x >= 300).sum())),
                           n_ge300_cyc100=("n_cycling", lambda x: 0), med_frac_cycling=("frac_cycling", "median"))
S["n_pat_ge300_and_cycling_ge150"] = R[(R.n_mal >= 300) & (R.n_cycling >= 150)].groupby("study").size()
S["n_pat_ge300_and_cycling_ge100"] = R[(R.n_mal >= 300) & (R.n_cycling >= 100)].groupby("study").size()
S = S.drop(columns="n_ge300_cyc100").fillna(0); print(S.to_string()); print(S.sum(numeric_only=True))
S.to_csv(OUT / "per_study_cycling_counts.csv")
print(R.frac_cycling.describe())
