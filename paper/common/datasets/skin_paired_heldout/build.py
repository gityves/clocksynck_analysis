"""Lesional vs non-lesional psoriatic skin of the same visit: the seven held-out series of A54.

Source: GEO series in `$SCCC_DATA/skin_lesional/processed_heldout/` (GSE13355, GSE137218, GSE153007, GSE193128,
GSE283313, GSE66511, GSE69967; parsed by legacy `paired_damping/scripts/build_skin_heldout.py`).
Legacy loader replaced: `paired_damping/scripts/a54_skin.py:load("heldout")` (psoriasis, at least 10 pairs).
Used by: A54 (H54.1-H54.3 replication), A78-A84 (psoriasis cohorts).

Choices (as legacy): unit = patient, condition = nonlesional / lesional, cohort = `series|disease`; GSE13355 is
log2 centred and GSE153007 a ratio to a reference (no absolute level): their expression-matched controls use the mean
log2 intensity of the GSE30999 non-lesional samples (`GSE30999_nonlesional_means.csv`), as the legacy `match_means`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "skin_lesional" / "processed_heldout"
LEVELS = ("nonlesional", "lesional")
NO_LEVEL = ("log2_centred", "glog2_ratio_to_reference")


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = T.lesional_series(SRC)
    samples = samples[samples.disease == "psoriasis"].reset_index(drop=True)
    return samples, T.units_table(samples)


def gene_means() -> dict[str, pd.Series]:
    """cohort -> mean log expression of the non-lesional samples (reference means for the two no-level series)."""
    ref = pd.read_csv(SRC / "GSE30999_nonlesional_means.csv", index_col=0).iloc[:, 0]
    out = {}
    for f in sorted(SRC.glob("GSE*.parquet")):
        d = T.read_all_genes(f, T.META)
        d = d[d.disease == "psoriasis"]
        non = d[d.status == "nonlesional"]
        non = non[non.patient.isin(set(d[d.status == "lesional"].patient))].set_index("patient").sort_index()
        genes = [c for c in non.columns if c not in T.META]
        if non.value_type.iloc[0] in NO_LEVEL:
            out[f"{f.stem}|psoriasis"] = ref[ref.index.isin(genes)]
        else:
            out[f"{f.stem}|psoriasis"] = T.mean_log_expression(non, T.META, non.value_type.iloc[0] == "counts")
    return out
