"""A127 feasibility (blind): data rule and stratum sizes. python -I feasibility.py"""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import numpy as np, pandas as pd
import a127_lib as L

m, X = L.load()
g2m_all = L.sh.hallmark_g2m()
rule = dict(n_clock_present=sum(g in X.columns for g in L.GENES), n_clock=len(L.GENES),
            n_g2m_present=sum(g in X.columns for g in g2m_all), n_g2m=len(g2m_all),
            n_purity_present=sum(g in X.columns for g in L.sh.PURITY), n_purity=len(L.sh.PURITY),
            clock_genes_missing=";".join(g for g in L.GENES if g not in X.columns),
            n_tumours=len(m), n_unique_patients=m.patient.nunique())
rule["g2m_fraction"] = rule["n_g2m_present"] / rule["n_g2m"]
pd.Series(rule).to_csv(L.OUT / "feasibility_data_rule.csv", header=False)
print(rule)
cl = X[L.GENES]
print("clock gene log2FPKM mean/sd:\n", pd.DataFrame(dict(mean=cl.mean(), sd=cl.std(), frac_floor=(cl < -3.3).mean())).round(2))
pd.DataFrame(dict(mean=cl.mean(), sd=cl.std(), frac_floor=(cl < -3.3).mean())).to_csv(L.OUT / "feasibility_clock_genes.csv")
m["pur_t"] = m.groupby("pam50").purity.transform(lambda s: L.tertile(s))
m["pur_t_cohort"] = L.tertile(m.purity)
t = m.groupby(["pam50", "pur_t"]).size().unstack().assign(total=lambda d: d.sum(1))
t["stratum_ge30"] = (m.groupby(["pam50", "pur_t"]).size() >= 30).groupby("pam50").sum()
t["tumours_per_prolif_tertile_min"] = (m.groupby(["pam50", "pur_t"]).size() // 3).groupby("pam50").min()
t.to_csv(L.OUT / "feasibility_pam50_purity_tertile.csv"); print(t)
e = m.groupby(["pur_t_cohort", "er"]).size().unstack(); e.to_csv(L.OUT / "feasibility_er_by_purity_tertile.csv"); print(e)
c = pd.crosstab(m.pam50, m.er); c.to_csv(L.OUT / "feasibility_pam50_by_er.csv"); print(c)
print(pd.crosstab(m.pam50, m.platform)); print(m[["age at diagnosis","tumor size","lymph node status","overall survival days","overall survival event"]].notna().sum())
print("events", m["overall survival event"].astype(float).sum())
