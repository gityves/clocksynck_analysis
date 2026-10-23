"""Inflamed vs non-inflamed gut mucosa of the same endoscopy (bulk): six GEO series of Crohn's disease and UC.

Source: `$SCCC_DATA/ibd_bulk/processed/` (GSE107499, GSE179285, GSE193677, GSE208303, GSE66407, GSE95095; parsed by
legacy `paired_damping/scripts/build_ibd.py`; PROVENANCE in `ibd_bulk/`).
Legacy loader replaced: `paired_damping/scripts/a55_ibd.py:load` (A55; `--disease CD|UC`, `--same-segment`).
Used by: A55, A79, A82, A84, A88.

Choices (as legacy): condition `nonlesional` = non-inflamed, `lesional` = inflamed (scout convention); unit = patient;
cohort = series (a disease filter is applied on `disease` before pairing; GSE179285 is a log2 ratio to a reference,
GSE193677 log2 adjusted counts, GSE208303 counts, the rest log2 intensity). `segment` and `same_segment` (inflamed and
non-inflamed biopsies from the same bowel segment) are kept. `gene_means()` returns the expression-matched control
input; GSE179285 has no absolute level and uses the GSE193677 non-inflamed means, as the legacy `match_means`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "ibd_bulk" / "processed"
LEVELS = ("nonlesional", "lesional")
EPI = ["EPCAM", "VIL1", "CDH17", "KRT20"]  # composition genes of A55
IMM = ["PTPRC", "CD3E", "CD68"]


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = T.lesional_series(SRC, by_disease=False)
    return samples, T.units_table(samples, ("disease", "segment", "platform"))


def gene_means() -> dict[str, pd.Series]:
    ref = pd.read_parquet(SRC / "GSE193677.parquet")
    ref = ref[ref.status == "nonlesional"]
    out = {}
    for f in sorted(SRC.glob("GSE*.parquet")):
        d = T.read_all_genes(f, T.META)
        non = d[d.status == "nonlesional"]
        non = non[non.patient.isin(set(d[d.status == "lesional"].patient))].set_index("patient").sort_index()
        genes = [c for c in non.columns if c not in T.META]
        if non.value_type.iloc[0] == "log2_ratio":
            rg = [c for c in ref.columns if c not in T.META]
            lm = ref[rg].astype(float).mean()
            out[f.stem] = lm[lm.index.isin(non.columns)]
        else:
            out[f.stem] = T.mean_log_expression(non, T.META, non.value_type.iloc[0] == "counts")
    return out
