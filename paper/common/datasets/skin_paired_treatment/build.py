"""Psoriasis biopsies before and after treatment, lesional and non-lesional skin of the same patient.

Source: `$SCCC_DATA/skin_lesional/processed_treatment/` (GSE117468 week 12, GSE201827 week 12 and week 52, array;
parsed by legacy `paired_damping/scripts/build_skin_treatment.py`).
Legacy loader replaced: `paired_damping/scripts/a56_treatment.py:visits` (and its copies in a57, a73, a74, a79-a80).
Used by: A56 (H56.2, interaction post hoc), A57, A74, A80.

Choices (as legacy): four samples per patient (visit baseline / post x status). `unit` = patient; the lesional vs
non-lesional pair of one visit is the pair of `clocksync.paired` when the table is filtered to one visit; the
baseline vs post comparison of A56 swaps visit labels within patient. `cohort` = file (series, `GSE201827_w52` is the
week-52 follow-up) |disease. Columns `visit`, `treatment`, `arm` as parsed (placebo arms contain "placebo").
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "skin_lesional" / "processed_treatment"
LEVELS = ("nonlesional", "lesional")


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    parts = []
    for f in sorted(SRC.glob("GSE*.parquet")):
        t = T.lesional_series(SRC, [f.name])
        t["cohort"] = f.stem + "|" + t.disease.astype(str)
        parts.append(t)
    samples = pd.concat(parts, ignore_index=True)
    return samples, T.units_table(samples, ("disease", "platform", "treatment", "arm"))
