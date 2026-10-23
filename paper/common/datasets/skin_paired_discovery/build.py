"""Lesional vs non-lesional skin of the same patient and visit: discovery series (bulk, microarray and RNA-seq).

Source: GEO series in `$SCCC_DATA/skin_lesional/processed/` (17 series: psoriasis and atopic dermatitis; sample
matrices parsed once by the legacy `paired_damping/scripts/build_skin.py`, all genes, one row per biopsy).
Legacy loaders replaced: `paired_damping/scripts/a53_skin.py:series` (A53, per series x disease) and
`a54_skin.py:load("discovery")` (A54: the nine psoriasis series `DISCOVERY_PSORIASIS`, at least 10 pairs).
Used by: A53, A54, A55-A57 baselines, A73, A74, A78-A84 (psoriasis cohorts), A88.

Choices (as legacy): unit = patient (both samples of a pair come from one visit), condition = nonlesional / lesional,
cohort = `series|disease` (GSE121212 holds psoriasis and AD), clock genes in panel symbols (BMAL1 -> ARNTL),
counts cohorts keep `depth` (genome-wide count), array cohorts (`value_type` log2_intensity) have depth NaN. Pairs are
formed by `clocksync.paired` (the legacy keeps the patients present in both statuses, at least 10 pairs).
`gene_means()` returns the expression-matched control input of the legacy (`match_means`: mean log fraction or mean
log2 intensity of the non-lesional samples, per gene).

    samples, units = build()
    cs.paired(samples, unit="unit", condition="condition", levels=LEVELS, cohort="cohort", counts="value_type",
              pseudo_means=gene_means())
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "skin_lesional" / "processed"
DISCOVERY_PSORIASIS = ["GSE30999", "GSE117468", "GSE201827", "GSE136757", "GSE85034", "GSE252029", "GSE197056",
                       "GSE268958", "GSE121212"]
LEVELS = ("nonlesional", "lesional")


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(samples, units): every sample of the 17 discovery series (all diseases)."""
    samples = T.lesional_series(SRC)
    return samples, T.units_table(samples)


def psoriasis_discovery(samples: pd.DataFrame) -> pd.DataFrame:
    """The A54 discovery set: psoriasis samples of the nine pre-registered series."""
    return samples[(samples.disease == "psoriasis") & samples.series.isin(DISCOVERY_PSORIASIS)]


def gene_means(samples: pd.DataFrame | None = None) -> dict[str, pd.Series]:
    """cohort -> mean log expression of all genes in the non-lesional samples (input of `paired(pseudo_means=)`)."""
    out = {}
    for f in sorted(SRC.glob("GSE*.parquet")):
        d = T.read_all_genes(f, T.META)
        for dis, g in d.groupby("disease"):
            non = g[g.status == "nonlesional"]
            les = set(g[g.status == "lesional"].patient)
            non = non[non.patient.isin(les)].set_index("patient").sort_index()
            out[f"{f.stem}|{dis}"] = T.mean_log_expression(non, T.META, non.value_type.iloc[0] == "counts")
    return out
