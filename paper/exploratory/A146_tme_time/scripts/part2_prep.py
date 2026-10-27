"""A146 part 2 blind preparation: inventory and power for GSE266919 (no program scores, no time or phase used).
Inputs: A140 Q2 sample table and cells-per-sample table (cell counts only). Outputs: output/part2_prep_*.csv."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

A140 = Path(__file__).resolve().parents[2] / "A140_dip_cyclops/output/q2"
OUT = Path(__file__).resolve().parents[1] / "output"
MIN = 50
s = pd.read_csv(A140 / "samples.csv").drop_duplicates("Sample").set_index("Sample")
c = pd.read_csv(A140 / "cells_per_sample.csv").set_index("Sample")
d = s.join(c, how="outer").fillna({k: 0 for k in c.columns})
d["T_all"] = d.CD4_T + d.CD8_T
comps = ["CD4_T", "CD8_T", "T_all", "NK", "myeloid", "B"]

# 1. inventory
inv = [("samples", len(d)), ("patients", d.Patient.nunique()),
       ("patients_with_2plus_samples", int((d.groupby("Patient").size() > 1).sum()))]
for col in ["Group", "Tissue"]:
    inv += [(f"samples_{col}={k}", int(v)) for k, v in d[col].value_counts().items()]
inv += [(f"patients_with_{k}", int(d[d.Group == k].Patient.nunique())) for k in d.Group.unique()]
inv += [("recorded_collection_time_fields", 0)]
pd.DataFrame(inv, columns=["item", "n"]).to_csv(OUT / "part2_prep_inventory.csv", index=False)

# compartments: patients with >=MIN cells, three patient definitions
rows = []
pre_breast = d[(d.Group.str.startswith("Pre")) & (d.Tissue == "breast")] if d.Tissue.isin(["breast"]).any() else d.iloc[:0]
sets = {"any_sample": d, "pre_treatment": d[d.Group.str.startswith("Pre")], "pre_breast": pre_breast}
for name, sub in sets.items():
    for k in comps:
        rows.append(dict(patient_set=name, compartment=k, n_samples=int((sub[k] >= MIN).sum()),
                         n_patients=int(sub[sub[k] >= MIN].Patient.nunique())))
for k in ["cancer_epithelium", "fibroblast", "endothelial", "mural"]:
    for name in sets:
        rows.append(dict(patient_set=name, compartment=k, n_samples=0, n_patients=0))
comp = pd.DataFrame(rows); comp.to_csv(OUT / "part2_prep_compartments.csv", index=False)

# 3. power (Fisher z, one-sided 0.05, 80 %)
def mde(n): return float(np.tanh((norm.ppf(.95) + norm.ppf(.8)) / np.sqrt(n - 3))) if n > 3 else np.nan
prog = [("P1", "cancer_epithelium"), ("P2", "cancer_epithelium"), ("P3", "cancer_epithelium"),
        ("P9", "cancer_epithelium"), ("P9", "myeloid"), ("P4", "T_all"), ("P10", "T_all"),
        ("P5", "myeloid"), ("P6", "myeloid"), ("P7", "all_cells_fraction"), ("P8", "all_cells_fraction")]
p = []
for pr, k in prog:
    for name in sets:
        n = int(comp[(comp.patient_set == name) & (comp.compartment == k)].n_patients.iloc[0]) \
            if k != "all_cells_fraction" else 0
        p.append(dict(program=pr, compartment=k, patient_set=name, n_patients_compartment=n,
                      n_patients_with_internal_time=0, mde_abs_r_if_time_existed=mde(n),
                      mde_abs_r_actual=np.nan,
                      note="no stromal host -> no internal time" if k != "all_cells_fraction"
                      else "only immune cells deposited; no all-cell denominator"))
pd.DataFrame(p).to_csv(OUT / "part2_prep_power.csv", index=False)
pd.DataFrame([dict(n=n, mde_abs_r=mde(n)) for n in [10, 20, 30, 40, 44, 60, 76]]).to_csv(OUT / "part2_prep_mde_table.csv", index=False)
print(pd.DataFrame(inv)); print(comp[comp.n_patients > 0].to_string()); print(pd.DataFrame(p).query("patient_set=='pre_breast'").to_string())
