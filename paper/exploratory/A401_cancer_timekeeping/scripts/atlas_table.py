"""Per-cancer atlas of the tumour clock readouts (descriptive assembly of scored outputs; no new statistic).

Sources (each value is copied from its scored output):
  A401 output/score/per_cancer.csv (14 genes): single-cell patients / cohorts, cancer-cell coupling c_mal (90 %), z_mal against
       matched sets, cancer-cell amplitude relative to fibroblasts A_mal, A138 bulk pairing, A83 clock-shaped share
  A403b output/a403b/d1.csv: immune-hot slope on tumour-stroma loudness per cancer (90 %)
  A304 output/score/a3_per_cancer.csv: TCGA pairs, tumour-normal correlation, coupling to body time c_T;
       output/score/per_cancer.csv: H304.1 gain high / low proliferation
  A304r output/a304r/per_type.csv: CPTAC-3 pairs and gain high / low
                                                     -> output/atlas/atlas.csv, atlas.md
    python atlas_table.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
EXP = HERE.parent
OUT = HERE / "output" / "atlas"
CPTAC = {"LUAD": "LUAD", "CCRCC": "KIRC", "LSCC": "LUSC", "SKCM": "SKCM", "HNSCC": "HNSC", "PDAC": "PAAD", "UCEC": "UCEC", "STAD": "STAD"}


def fmt(v, lo=None, hi=None, d=2):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return ""
    if not np.isfinite(v):
        return ""
    s = f"{v:.{d}f}"
    if lo is not None and np.isfinite(lo) and np.isfinite(hi):
        s += f" ({lo:.{d}f} to {hi:.{d}f})"
    return s


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(HERE / "output/score/per_cancer.csv")
    a = a[a.readout.astype(str) == "14"].set_index("cancer")
    d1 = pd.read_csv(HERE / "output/a403b/d1.csv").query("readout == 'immune_hot'").set_index("cancer")
    a3 = pd.read_csv(EXP / "A304_cancer_gain_departure/output/score/a3_per_cancer.csv").set_index("cancer")
    pc = pd.read_csv(EXP / "A304_cancer_gain_departure/output/score/per_cancer.csv").query("contrast == 'H304.1'").set_index("cancer")
    cp = pd.read_csv(EXP / "A304_cancer_gain_departure/output/a304r/per_type.csv")
    cp["cancer"] = cp.type.map(CPTAC)
    cp = cp.set_index("cancer")
    cancers = sorted((set(a.index) | set(a3.index) | set(cp.index)) - {"unspecified", "kidney_unspecified", "SCLC", "SARC"})
    rows = []
    for c in cancers:
        r = dict(cancer=c)
        if c in a.index:
            x = a.loc[c]
            r.update(sc_patients=int(x.patients), sc_cohorts=int(x.cohorts), c_mal=x.c_mal, c_mal_lo=x.c_lo90, c_mal_hi=x.c_hi90,
                     z_mal=x.z_mal, A_mal=x.A_mal, bulk_pairing_A138=x.bulk_pairing, a83_share=x.a83_share)
        if c in d1.index:
            x = d1.loc[c]
            r.update(colder_louder_slope=x.slope, cl_lo=x.lo90, cl_hi=x.hi90, cl_patients=int(x.patients))
        if c in a3.index:
            x = a3.loc[c]
            r.update(tcga_pairs=int(x.n), tcga_corr=x["corr"], c_body=x["c_T"])
        if c in pc.index:
            x = pc.loc[c]
            r.update(tcga_gain_ratio=x.G_case / x.G_ref if x.G_ref > 0 else np.nan)
        if c in cp.index:
            x = cp.loc[c]
            r.update(cptac_pairs=int(x.n), cptac_gain_ratio=x.G_high / x.G_low if x.G_low > 0 else np.nan)
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "atlas.csv", index=False)
    md = ["| cancer | single-cell patients (cohorts) | cancer-cell coupling c (90 %) | z vs matched | cancer-cell amplitude / fibroblast | colder-louder slope (90 %) | TCGA pairs | TCGA tumour-normal corr | coupling to body time | TCGA gain ratio high / low prolif | CPTAC pairs | CPTAC gain ratio |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in T.iterrows():
        g = r.get
        md.append("| " + " | ".join([
            r.cancer,
            f"{int(g('sc_patients'))} ({int(g('sc_cohorts'))})" if pd.notna(g("sc_patients")) else "",
            fmt(g("c_mal"), g("c_mal_lo"), g("c_mal_hi")), fmt(g("z_mal")), fmt(g("A_mal")),
            fmt(g("colder_louder_slope"), g("cl_lo"), g("cl_hi")),
            str(int(g("tcga_pairs"))) if pd.notna(g("tcga_pairs")) else "", fmt(g("tcga_corr")), fmt(g("c_body")),
            fmt(g("tcga_gain_ratio")), str(int(g("cptac_pairs"))) if pd.notna(g("cptac_pairs")) else "", fmt(g("cptac_gain_ratio"))]) + " |")
    (OUT / "atlas.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
