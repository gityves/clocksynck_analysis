"""Lesional vs uninvolved psoriatic skin resolved by layer: Visium (GSE206391, GSE202011) and GeoMx (GSE310145).

Source: `$SCCC_DATA/skin_layers/processed/{GSE206391,GSE202011,GSE310145}.parquet` (marker-defined layers for the
Visium series, author segments for GeoMx; parsed by legacy `paired_damping/scripts/build_skin_layers.py`; counts).
Legacy loader replaced: `paired_damping/scripts/a62_layers.py:load` (A62), also read by `a77_geomx_layers.py` (A77).
Used by: A62, A77.

`build()` returns the unaggregated table (one row per spot group or GeoMx segment): unit = patient, condition = status,
cohort = series, `layer` (epidermis / dermis), `region` (GeoMx segment type), `site`. `layer_pairs()` applies the
A62 aggregation: rows of a layer (GeoMx: epidermis = epidermis + basal_layer segments, dermis = vim_high + vim_low;
derm_pap, lesional only, is left out) are summed per patient x status (counts and depth), giving the paired table
of that layer (one row per patient x status, ready for `clocksync.paired(cohort="cohort")`).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "skin_layers" / "processed"
SERIES = ["GSE206391", "GSE202011", "GSE310145"]
LEVELS = ("nonlesional", "lesional")
GEOMX = {"epidermis": ("epidermis", "basal_layer"), "dermis": ("vim_high", "vim_low")}


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = T.lesional_series(SRC, [f"{g}.parquet" for g in SERIES], by_disease=False)
    return samples, T.units_table(samples, ("disease", "platform"))


def layer_pairs(samples: pd.DataFrame, layer: str, series=SERIES, pso_only: bool = False, regions=None) -> pd.DataFrame:
    """A62 aggregation of one layer: counts and depth summed per cohort x unit x condition."""
    d = samples[samples.layer == layer]
    if pso_only:
        d = d[d.disease != "PsA"]
    g = d[d.series != "GSE310145"]
    geo = d[d.series == "GSE310145"]
    geo = geo[geo.region.isin(regions or GEOMX[layer])]
    d = pd.concat([g, geo])
    d = d[d.series.isin(series)]
    num = ["depth"] + T.GENES
    x = d.groupby(["cohort", "unit", "condition"], as_index=False)[num].sum()
    x["value_type"] = "counts"
    x["layer"] = layer
    return x
