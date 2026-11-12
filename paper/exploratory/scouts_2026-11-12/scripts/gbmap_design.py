"""Scout (design only): per study x broad class design table for the Extended GBmap. Counts and UMI totals only.

Study key = author + suspension (cell / nucleus). Patient = donor_id (no donor occurs in two authors of the atlas).
Cells are restricted to a chosen assay set: 'droplet' (10x 3'/5', Drop-seq) and 'umi_all' (droplet + microwell-seq,
BD Rhapsody, CEL-seq2, STRT-seq; everything except Smart-seq2).
Class mapping (published annotation_level_3, cell_type; malignant = annotation_level_1 Neoplastic):
  myeloid = TAM-BDM, TAM-MG, Mono, DC, Neutrophil, Mast; oligodendrocyte = Oligodendrocyte, OPC (non-malignant);
  T/NK = CD4/CD8, NK; endothelial; mural = Mural cell; neuron; astrocyte; other = B cell, Plasma B, RG.
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "output"
import h5py
from anndata.io import read_elem

SRC = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_gbmap/extended_gbmap_230c8701.h5ad"
cols = ["assay", "suspension_type", "author", "donor_id", "annotation_level_1", "annotation_level_3", "stage", "MGMT"]
with h5py.File(SRC, "r") as f:
    obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in cols}, index=read_elem(f["obs"]["_index"]))
cnt = pd.read_csv(OUT / "gbmap_cell_counts.csv.gz", index_col=0)
assert (cnt.index == obs.index).all()
obs = obs.join(cnt)
for c in cols:
    obs[c] = obs[c].astype(str)

MAP = {"TAM-BDM": "myeloid", "TAM-MG": "myeloid", "Mono": "myeloid", "DC": "myeloid", "Neutrophil": "myeloid",
       "Mast": "myeloid", "Oligodendrocyte": "oligodendrocyte", "OPC": "oligodendrocyte", "CD4/CD8": "T_NK",
       "NK": "T_NK", "Endothelial": "endothelial", "Mural cell": "mural", "Neuron": "neuron",
       "Astrocyte": "astrocyte", "B cell": "other", "Plasma B": "other", "RG": "other"}
obs["cls"] = np.where(obs.annotation_level_1 == "Neoplastic", "malignant", obs.annotation_level_3.map(MAP))
assert obs.cls.notna().all()
obs["study"] = obs.author + "|" + obs.suspension_type
DROPLET = lambda a: a.startswith("10x") or a == "Drop-seq"
obs["droplet"] = obs.assay.map(DROPLET)
obs["umi"] = obs.assay != "Smart-seq2"
CLASSES = ["malignant", "myeloid", "oligodendrocyte", "T_NK", "endothelial", "mural", "neuron", "astrocyte", "other"]
OTHERS = CLASSES[1:]
MIN = 100


def design(d, label):
    pb = d.groupby(["study", "donor_id", "cls"]).agg(n=("clock_umi", "size"), clock=("clock_umi", "sum"),
                                                     mgmt=("mgmt_umi", "sum")).reset_index()
    pb_all = pb.groupby(["donor_id", "cls"], as_index=False)[["n", "clock", "mgmt"]].sum()  # a patient counts once
    rows = []
    for study, g in list(pb.groupby("study")) + [("ALL", pb_all)]:
        mal = set(g[(g.cls == "malignant") & (g.n >= MIN)].donor_id)
        for c in CLASSES:
            x = g[g.cls == c]
            big = x[x.n >= MIN]
            r = dict(set=label, study=study, cls=c, n_patients_any=x.donor_id.nunique(), n_patients_ge100=len(big),
                     n_patients_ge100_with_malignant=len(set(big.donor_id) & mal) if c != "malignant" else len(big),
                     median_cells_per_patient=x.n.median() if len(x) else np.nan,
                     median_clock_umi_pseudobulk=big.clock.median() if len(big) else np.nan,
                     median_clock_umi_per_cell=(big.clock / big.n).median() if len(big) else np.nan)
            if c == "malignant":
                r["median_mgmt_umi_pseudobulk"] = big.mgmt.median() if len(big) else np.nan
                r["frac_patients_mgmt_ge10"] = (big.mgmt >= 10).mean() if len(big) else np.nan
            rows.append(r)
    return pd.DataFrame(rows)


res = []
for label, mask in [("droplet", obs.droplet), ("umi_all", obs.umi)]:
    res.append(design(obs[mask], label))
res = pd.concat(res)
res.to_csv(OUT / "gbmap_design_table.csv", index=False)

# patient-level table (public clinical metadata of the atlas; counts only)
d = obs[obs.droplet]
pt = d.groupby(["study", "donor_id"]).agg(cells=("cls", "size"), malignant=("cls", lambda s: (s == "malignant").sum()),
                                           stage=("stage", lambda s: s.mode().iat[0]),
                                           mgmt_labels=("MGMT", lambda s: ",".join(sorted(set(s))))).reset_index()
for c in OTHERS:
    pt[c] = d[d.cls == c].groupby(["study", "donor_id"]).size().reindex(pd.MultiIndex.from_frame(pt[["study", "donor_id"]])).fillna(0).astype(int).values
pt.to_csv(OUT / "gbmap_patient_table_droplet.csv", index=False)

# patient-level summary of metadata for patients with >=100 malignant droplet cells
ok = pt[pt.malignant >= MIN].copy()
norm = {"M": "methylated", "Meth": "methylated", "UM": "unmethylated", " UM": "unmethylated", "Un": "unmethylated",
        "NM": "unmethylated"}


def mg(s):
    v = {norm[x] for x in s.split(",") if x in norm}
    return v.pop() if len(v) == 1 else ("conflict" if v else "unknown")


ok["mgmt_status"] = ok.mgmt_labels.map(mg)
summ = pd.concat([ok.stage.value_counts().rename("n").rename_axis("value").reset_index().assign(variable="stage"),
                  ok.mgmt_status.value_counts().rename("n").rename_axis("value").reset_index().assign(variable="mgmt_status")])
summ.to_csv(OUT / "gbmap_metadata_summary_ge100malignant_droplet.csv", index=False)
print(res[res.study == "ALL"].to_string())
print(summ.to_string())

# breadth: malignant patients (>=100 droplet malignant cells, patient counted once) by number of other classes with >=100 cells
pc = d.groupby(["donor_id", "cls"]).size().unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0)
pc = pc[pc.malignant >= MIN]
brd = pd.Series({"patients_ge100_malignant": len(pc),
                 "myeloid": int((pc.myeloid >= MIN).sum()),
                 "myeloid_and_oligodendrocyte": int(((pc.myeloid >= MIN) & (pc.oligodendrocyte >= MIN)).sum()),
                 "myeloid_oligo_TNK": int(((pc.myeloid >= MIN) & (pc.oligodendrocyte >= MIN) & (pc.T_NK >= MIN)).sum()),
                 "n_other_classes_ge100_0": int(((pc[OTHERS] >= MIN).sum(axis=1) == 0).sum()),
                 "n_other_classes_ge100_1": int(((pc[OTHERS] >= MIN).sum(axis=1) == 1).sum()),
                 "n_other_classes_ge100_2": int(((pc[OTHERS] >= MIN).sum(axis=1) == 2).sum()),
                 "n_other_classes_ge100_3plus": int(((pc[OTHERS] >= MIN).sum(axis=1) >= 3).sum())}, name="n")
brd.to_csv(OUT / "gbmap_breadth_droplet.csv")
print(brd.to_string())
