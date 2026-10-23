"""Lesional vs non-lesional skin of other conditions (acne, alopecia areata, hidradenitis suppurativa, vitiligo,
discoid lupus, actinic keratosis, mycosis fungoides; morphea, sarcoidosis and melasma on a linear scale).

Source: GEO series in `$SCCC_DATA/skin_lesional/processed_other/` (14 series) and `processed_other_linear/`
(3 series, log2 of linear intensity; legacy: outside the A59 rule, descriptive only; column `source` tells them apart).
Parsed by legacy `paired_damping/scripts/build_skin_other.py`.
Legacy loader replaced: `paired_damping/scripts/a59_other.py:load` (per series x condition, at least 10 pairs).
Used by: A59, A79 (companion), A82, A84, A88.

Choices (as legacy): unit = patient, condition = nonlesional / lesional, cohort = `series|disease` (the legacy names it
`series_disease`). Mycosis fungoides is `disease == "mycosis_fungoides"` (A59 sensitivity without cutaneous lymphoma).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

ROOT = DATA / "skin_lesional"
LEVELS = ("nonlesional", "lesional")
FOLDERS = {"standard": "processed_other", "linear": "processed_other_linear"}


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    parts = []
    for tag, folder in FOLDERS.items():
        t = T.lesional_series(ROOT / folder)
        t.insert(3, "source", tag)
        parts.append(t)
    samples = pd.concat(parts, ignore_index=True)
    return samples, T.units_table(samples)


def gene_means() -> dict[str, pd.Series]:
    """cohort -> mean log expression of the non-lesional samples; series without an absolute level (log2 centred,
    ratio to a reference: GSE207744) use the GSE30999 non-lesional means, as the legacy `match_means`."""
    ref = pd.read_csv(ROOT / "processed_heldout" / "GSE30999_nonlesional_means.csv", index_col=0).iloc[:, 0]
    out = {}
    for folder in FOLDERS.values():
        for f in sorted((ROOT / folder).glob("GSE*.parquet")):
            d = T.read_all_genes(f, T.META)
            for dis, g in d.groupby("disease"):
                non = g[g.status == "nonlesional"]
                non = non[non.patient.isin(set(g[g.status == "lesional"].patient))].set_index("patient").sort_index()
                if non.value_type.iloc[0] in ("log2_centred", "glog2_ratio_to_reference"):
                    genes = [c for c in non.columns if c not in T.META]
                    out[f"{f.stem}|{dis}"] = ref[ref.index.isin(genes)]
                else:
                    out[f"{f.stem}|{dis}"] = T.mean_log_expression(non, T.META, non.value_type.iloc[0] == "counts")
    return out
