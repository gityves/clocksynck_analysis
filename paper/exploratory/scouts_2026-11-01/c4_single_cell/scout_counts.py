"""C4 scout: per-study / per-patient malignant-cell counts, UMI depth and clock-gene detection (expression counts only; no phase or
coupling statistic). Reads the stored pooled table (cancer_sc) and the legacy CoPhaser cell tables + copy-number calls.
Run: python -I scout_counts.py (sccc env, <= 4 threads)."""
import sys, gzip, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "/home/maxine/Documents/paychere/scCircadianOverview/paper")
from common import data
OUT = Path(__file__).resolve().parent
D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
CC = Path("/home/maxine/Documents/paychere/scCircadianCancer/outputs/states/cc")
G = "ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR2D1 NR1D2 DBP RORA TEF CIART".replace(" NR2D1", "").split()

pb, don = data.load("cancer_sc")
mal = pb[pb.compartment == "tumour:malignant"].copy()
mal["patient"] = mal.donor.map(don.patient)
mal["cancer"] = mal.donor.map(don.cancer)
mal["umi_mean"] = mal.depth / mal.n_cells
mal.drop(columns=G).to_csv(OUT / "malignant_per_patient_pseudobulk_meta.csv", index=False)  # no clock counts kept
paired = set()
try:
    import importlib; sys.path.insert(0, str(Path("/home/maxine/Documents/paychere/scCircadianOverview/paper/results/figure_5/scripts/cancer")))
    import _shared as sh
    paired = set(sh.sc_builder().paired_epithelium(pb))
except Exception as e:
    print("paired failed", e)
mal["a87"] = mal.donor.isin(paired)
rows = []
for st, g in mal.groupby("study"):
    for lab, h in (("all", g), ("a87", g[g.a87])):
        if len(h) == 0: continue
        rows.append(dict(study=st, cancer=h.cancer.iloc[0], set=lab, n_patients=len(h), n_mal_total=int(h.n_cells.sum()),
                         med_cells=h.n_cells.median(), min_cells=h.n_cells.min(),
                         n_ge300=int((h.n_cells >= 300).sum()), n_ge500=int((h.n_cells >= 500).sum()),
                         n_ge1000=int((h.n_cells >= 1000).sum()), mean_umi_per_cell_median_patient=h.umi_mean.median()))
S = pd.DataFrame(rows); S.to_csv(OUT / "per_study_malignant_counts.csv", index=False)
print(S.to_string())

# cell level: detection and median UMI of malignant cells (copy-number call), CoPhaser cell tables
det = []
for f in sorted(CC.glob("*.csv.gz")):
    st = f.name.split(".")[0]
    cf = D / "cnv" / f"{st}.csv.gz"
    hdr = pd.read_csv(f, nrows=0).columns
    t = pd.read_csv(f, usecols=[hdr[0], "sample", "patient", "kind", "compartment", "total_umi"] + [f"clock:{g}" for g in G], index_col=0)
    t = t[t.kind == "tumour"]
    if cf.exists():
        c = pd.read_csv(cf, usecols=["sample", "barcode", "cnv_call"])
        c.index = c["sample"].astype(str) + "_" + c["barcode"].astype(str)
        t = t.join(c.cnv_call, how="inner"); m = t[t.cnv_call == "malignant"]; lab = "cnv malignant"
    else:
        m = t[t.compartment.astype(str).str.contains("epithel")]; lab = "epithelial (no cnv call)"
    if len(m) == 0:
        det.append(dict(study=st, label="none")); continue
    r = dict(study=st, label=lab, n_cells_in_cycle_table=len(m), n_patients=m.patient.nunique(), median_umi=m.total_umi.median())
    for g in G: r["det_" + g] = (m[f"clock:{g}"] >= 1).mean()
    r["det_mean14"] = np.mean([r["det_" + g] for g in G])
    r["frac_cells_any_clock_umi"] = (m[[f"clock:{g}" for g in G]].sum(1) >= 1).mean()
    r["median_clock_umi_per_cell"] = m[[f"clock:{g}" for g in G]].sum(1).median()
    pp = m.groupby("patient").size()
    r["n_pat_ge300_in_table"] = int((pp >= 300).sum()); r["median_cells_per_pat_table"] = pp.median(); r["min_cells_per_pat_table"] = pp.min()
    det.append(r)
    m.groupby("patient").size().rename("n_mal_cells_in_cycle_table").reset_index().assign(study=st).to_csv(
        OUT / f"_tmp_{st}.csv", index=False)
Dd = pd.DataFrame(det); Dd.to_csv(OUT / "per_study_cell_level_detection.csv", index=False)
pd.set_option("display.width", 250); print(Dd.round(3).to_string())
parts = [pd.read_csv(p) for p in OUT.glob("_tmp_*.csv")]
pd.concat(parts).to_csv(OUT / "per_patient_malignant_cells_in_cycle_tables.csv", index=False)
for p in OUT.glob("_tmp_*.csv"): p.unlink()
