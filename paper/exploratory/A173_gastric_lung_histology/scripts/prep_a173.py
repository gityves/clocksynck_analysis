"""A173 data preparation (blind to the clock statistic): GEO series matrices of the gastric and lung-histology cohorts
(scout `scouts_2026-10-31/gastric_lung_histology`) parsed with the A171 parser into ../output/{key}_meta.parquet and
{key}_expr.parquet (log2, probes averaged per symbol; G2M / purity-proxy scores over the kept tumours). No coherence here.
Run: python -I prep_a173.py   (reads downloaded text files only; nothing downloaded is executed)"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A171_luad_lihc" / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
import a168_lib as L  # noqa: E402
import prep_a171 as Q  # noqa: E402

OUT = HERE.parent / "output"
OUT.mkdir(exist_ok=True)
L.OUT = OUT
Q.P.D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_gastric_lung")
Q.META = HERE.parents[1] / "scouts_2026-10-31/gastric_lung_histology/meta"
Q.SYMBOL_COLUMNS = (*Q.SYMBOL_COLUMNS, "Symbol")  # Illumina DASL full table (GPL8432)
LOG10_SERIES = {"GSE62254"}  # deposited log10 (0.57-4.41): x log2(10) (registered)


def series(gse, key, gpl, keep, cols):
    m = pd.read_csv(Q.META / f"{gse}_samples.csv", dtype=str)
    m = m[keep(m)].set_index("gsm")
    m = m[[c for c in cols if c in m]].rename(columns=cols)
    ex = Q.expression_file(gse, f"{gse}_series_matrix.txt.gz", gpl)
    if gse in LOG10_SERIES:
        ex = ex * np.log2(10)
        Q.P.LOG.append((gse, "deposited_scale", "log10 -> log2"))
    Q.P.finish(gse, key, m, ex)


if __name__ == "__main__":
    allrows = lambda m: m.gsm.notna()
    hist = lambda *v: (lambda m: m.histology.str.strip().isin(v))
    series("GSE62254", "gse62254", "GPL570", allrows, {"patient": "patient"})
    series("GSE84437", "gse84437", "GPL6947", allrows, {"ptstage": "pt", "pnstage": "pn"})
    series("GSE15459", "gse15459", "GPL570", allrows, {})
    series("GSE26253", "gse26253", "GPL8432", allrows, {"pathological stage": "stage"})
    series("GSE37745", "gse37745", "GPL570", hist("adeno", "squamous"), {"histology": "histology", "tumor stage": "stage"})
    series("GSE41271", "gse41271", "GPL6884", hist("Adenocarcinoma", "Squamous"), {"histology": "histology", "final patient stage": "stage"})
    series("GSE30219", "gse30219", "GPL570", hist("ADC", "SQC"), {"histology": "histology", "pt stage": "pt"})
    pd.DataFrame(Q.P.LOG, columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)
    for k in ("gse62254", "gse84437", "gse15459", "gse26253", "gse37745", "gse41271", "gse30219"):
        m = pd.read_parquet(OUT / f"{k}_meta.parquet")
        print(k, len(m), m.histology.value_counts().to_dict() if "histology" in m else "", flush=True)
