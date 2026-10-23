"""Timed mouse skin: whole epidermis and epidermal stem cells over the day (8 GEO series, arrays and RNA-seq).

Source: `$SCCC_DATA/mouse_skin_timed/timed_<GSE>.parquet` (GSE84510, GSE84511, GSE84521 stem cells, microarray;
GSE190032, GSE190034, GSE174155, GSE294336, GSE294339 whole epidermis; mouse symbols mapped to the human clock panel;
parsed by legacy `paired_damping/scripts/build_timed_epidermis.py`; PROVENANCE per series in `mouse_skin_timed/`).
Legacy loader replaced: `paired_damping/scripts/a72_mouse_timed.py` (`pd.read_parquet(timed_<GSE>.parquet)`).
Used by: A72 (PRIMARY and STRICT sets below), A78 (timed reference cue sets), descriptive per-gene acrophases.

Table: one row per sample, `unit` = sample, `condition` = population (stem / whole_epidermis), `cohort` = series,
`zt` (h, modulo 24 as the legacy), lighting (LD, DD), age, sex, `value_type` (counts for GSE294336 / GSE294339;
the others are arrays or the authors' log2 normalised values, depth NaN or the library size for counts).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "mouse_skin_timed"
PRIMARY = ["GSE84510", "GSE84511", "GSE84521", "GSE190032", "GSE294336", "GSE294339"]  # A72 eligibility: LD, non-aged
STRICT = ["GSE84510", "GSE84511", "GSE190032", "GSE294336"]  # young adults only
SERIES = ["GSE84510", "GSE84511", "GSE84521", "GSE174155", "GSE190032", "GSE190034", "GSE294336", "GSE294339"]
META = ["sample", "population", "rep_in_time", "zt", "lighting", "age", "sex", "platform", "value_type", "depth"]


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    parts = []
    for g in SERIES:
        d = T.read_clock(SRC / f"timed_{g}.parquet", META).reset_index(drop=True)
        d["zt"] = d.zt.astype(float) % 24
        d.insert(0, "cohort", g)
        d.insert(0, "unit", d["sample"].astype(str))
        d.insert(1, "condition", d.population.astype(str))
        parts.append(d)
    t = pd.concat(parts, ignore_index=True)
    lead = ["unit", "condition", "cohort", "zt", "value_type", "depth"]
    t = t[lead + [c for c in t.columns if c not in lead + T.GENES] + T.GENES]
    units = t.set_index(t.cohort + "|" + t.unit)[["zt", "lighting", "age", "sex", "platform"]].rename_axis("unit_key")
    return t, units
