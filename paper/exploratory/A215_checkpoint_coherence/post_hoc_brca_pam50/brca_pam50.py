"""POST HOC (2026-11-06, user question "are METABRIC and TCGA opposite because of cancer type?"), descriptive.
TCGA BRCA alone, the three genomic contrasts of A210 (MYC amplified 2 vs neutral 0), A213 (FGA top vs bottom tertile) and
A215 (TP53 altered vs intact), with two cell definitions: proliferation tertile only (the registered TCGA cells) and
PAM50 x proliferation tertile (the METABRIC cells). >= 10 tumours per arm; labels permuted within cell (2,000); one-sided p
for a loss. Coherence, expression and proliferation from the A210 cache; PAM50 from the cBioPortal PanCanAtlas BRCA study.
    python -I brca_pam50.py
Output: brca_pam50.csv"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "A168_endocrine" / "scripts")); sys.path.insert(0, str(EXP.parents[1] / "paper"))
import a168_lib as L  # noqa: E402
from common.paths import DATA  # noqa: E402

GENES, tertile = list(L.GENES), L.A166.tertile
d = pd.read_parquet(EXP / "A210_myc_coherence/output/tcga_cache/BRCA.parquet").reset_index(drop=True)
pam = pd.read_csv(DATA / "tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep="\t", comment="#")
d["pam50"] = d.patient.map(pam.set_index("PATIENT_ID").SUBTYPE)
d["myc"] = d["sample"].map(pd.read_csv(EXP / "A210_myc_coherence/output/myc_cna_tcga.csv", index_col=0).myc_cna)
d["fga"] = d["sample"].map(pd.read_csv(EXP / "A213_heterogeneity_coherence/output/fga_tcga.csv", index_col=0).fga)
d["tp53"] = d["sample"].map(pd.read_csv(EXP / "A215_checkpoint_coherence/output/status_tcga.csv", index_col=0).tp53)
Y = np.log((d[GENES].to_numpy(float) + 0.5) / d.depth.to_numpy(float)[:, None])


def arms(g, meas):
    if meas == "myc":
        return g[g.myc == 2].index.to_numpy(), g[g.myc == 0].index.to_numpy()
    if meas == "tp53":
        return g[g.tp53 == 1].index.to_numpy(), g[g.tp53 == 0].index.to_numpy()
    t = tertile(g.fga.to_numpy(float))
    return g.index.to_numpy()[t == 2], g.index.to_numpy()[t == 0]


rows, rng = [], np.random.default_rng(215)
for meas in ("myc", "fga", "tp53"):
    m = d[d[meas].notna()].copy()
    for strata in ("prolif", "pam50_prolif"):
        mm = m[m.pam50.notna()].copy() if strata == "pam50_prolif" else m.copy()
        key = ["pam50"] if strata == "pam50_prolif" else []
        mm["t"] = mm.groupby(key).prolif.transform(tertile) if key else tertile(mm.prolif.to_numpy(float))
        cells = []
        for _, g in mm.groupby(key + ["t"]):
            a, b = arms(g, meas)
            if min(len(a), len(b)) >= 10:
                cells.append((a, b))
        T = lambda cs: np.mean([L.coh(Y, a, L.WF) - L.coh(Y, b, L.WF) for a, b in cs])
        obs = T(cells)
        null = []
        for _ in range(2000):
            cs = []
            for a, b in cells:
                p = rng.permutation(np.r_[a, b]); cs.append((p[:len(a)], p[len(a):]))
            null.append(T(cs))
        null = np.array(null)
        rows.append(dict(measure=meas, cells_by=strata, n_cells=len(cells), n_tumours=int(sum(len(a) + len(b) for a, b in cells)), T=obs,
                         null_mean=null.mean(), null_sd=null.std(ddof=1), excess=obs - null.mean(), p_one_lower=(1 + np.sum(null <= obs)) / (1 + len(null))))
    print(rows[-2:], flush=True)
# PAM50 composition of the arms (pooled over proliferation tertiles)
comp = []
for meas in ("myc", "fga", "tp53"):
    m = d[d[meas].notna() & d.pam50.notna()]
    a, b = arms(m, meas)
    for arm, idx in (("altered_or_top", a), ("reference", b)):
        comp.append(dict(measure=meas, arm=arm, n=len(idx), **m.loc[idx].pam50.value_counts(normalize=True).round(3).to_dict()))
out = pd.DataFrame(rows); out.to_csv(HERE / "brca_pam50.csv", index=False); pd.DataFrame(comp).to_csv(HERE / "brca_pam50_composition.csv", index=False)
print(out.round(4).to_string()); print(pd.DataFrame(comp).to_string())
