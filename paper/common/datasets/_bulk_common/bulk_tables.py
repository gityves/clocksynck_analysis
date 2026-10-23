"""Helpers shared by the bulk and paired dataset builders (not a dataset): sample tables for `clocksync.paired`.

A paired table has one row per sample: `unit` (patient, resection or visit shared by the two samples of a pair),
`condition` (e.g. nonlesional / lesional), `cohort` (series or cancer: pairs are pooled over cohorts), `value_type`
(counts, or the array scale), `depth` (genome-wide count; NaN for arrays), the 14 clock genes in `clocksync`
symbols (BMAL1 -> ARNTL), and descriptive columns. Counts are raw counts, arrays are log2 intensities as in the
processed files; `paired(..., counts="value_type")` mixes them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import clocksync as cs
from clocksync.genes import canonical, gene_columns

GENES = list(cs.gtex().genes)


def clock_columns(names) -> list[str]:
    """Names (as stored) of the 14 clock genes in `names`, in panel order."""
    cols = gene_columns(list(names), GENES)
    if (cols < 0).any():
        raise KeyError([g for g, c in zip(GENES, cols) if c < 0])
    return [list(names)[i] for i in cols]


def read_clock(path, meta: list[str], extra: tuple[str, ...] = ()) -> pd.DataFrame:
    """Read the metadata columns present in `meta`+`extra` and the 14 clock genes of a processed parquet (all-gene
    columns are not read), clock genes renamed to the panel symbols."""
    names = pq.ParquetFile(path).schema.names
    gc = clock_columns(names)
    keep = [c for c in list(meta) + list(extra) if c in names]
    d = pd.read_parquet(path, columns=keep + gc)
    return d.rename(columns=dict(zip(gc, GENES)))


def read_all_genes(path, meta: list[str]) -> pd.DataFrame:
    """All gene columns of a processed parquet (for expression-matched controls), BMAL1 renamed to ARNTL."""
    d = pd.read_parquet(path)
    d = d.rename(columns={"BMAL1": "ARNTL"}) if "ARNTL" not in d else d
    return d


def mean_log_expression(non: pd.DataFrame, meta: list[str], counts: bool, depth=None) -> pd.Series:
    """Legacy `match_means`: mean log fraction (counts, genes with mean count >= 1) or mean log2 intensity of the
    reference samples, per gene (paired_damping/a53_skin.py, a54_skin.py)."""
    genes = [c for c in non.columns if c not in meta]
    x = non[genes].astype(float)
    if counts:
        d = non.depth.to_numpy(float)
        return np.log(x.div(d, axis=0).mean() + 1e-12)[x.mean() >= 1]
    return x.mean()


def pairs(samples: pd.DataFrame, levels, cohort="cohort", unit="unit", condition="condition") -> pd.DataFrame:
    """Number of complete pairs per cohort (a unit with both conditions)."""
    out = {}
    for k, g in samples.groupby(cohort, sort=False):
        a = set(g.loc[g[condition] == levels[0], unit])
        b = set(g.loc[g[condition] == levels[1], unit])
        out[k] = len(a & b)
    return pd.Series(out, name="n_pairs")


META = ["patient", "status", "sample", "dataset", "disease", "platform", "value_type", "depth", "n_replicates", "visit",
        "treatment", "arm", "layer", "region", "site", "segment", "same_segment", "n_units"]


def lesional_series(folder: Path, files=None, extra: tuple[str, ...] = (), by_disease: bool = True) -> pd.DataFrame:
    """One sample table from the processed parquets of `folder` (columns patient, status, ..., genes).

    unit = patient, condition = status (nonlesional / lesional; for the gut: non-inflamed / inflamed), series = the
    GEO accession, cohort = series|disease (series if by_disease is False). Rows keep the file order, files sorted.
    """
    files = sorted(Path(folder).glob("GSE*.parquet")) if files is None else [Path(folder) / f for f in files]
    out = []
    for f in files:
        d = read_clock(f, META, extra)
        d.insert(0, "series", f.stem.split("_")[0])
        d["cohort"] = d.series + "|" + d.disease.astype(str) if by_disease else d.series
        d = d.rename(columns={"patient": "unit", "status": "condition"})
        d["unit"] = d["unit"].astype(str)
        out.append(d)
    t = pd.concat(out, ignore_index=True)
    lead = ["unit", "condition", "cohort", "series", "disease", "value_type", "depth"]
    rest = [c for c in t.columns if c not in lead + GENES]
    return t[lead + rest + GENES]


def units_table(samples: pd.DataFrame, cols=("disease", "platform", "dataset")) -> pd.DataFrame:
    """One row per (cohort, unit), index `cohort|unit`: first value of the descriptive columns."""
    cols = [c for c in cols if c in samples]
    u = samples.drop_duplicates(["cohort", "unit"]).set_index(["cohort", "unit"])[cols]
    u.index = [f"{a}|{b}" for a, b in u.index]
    return u.rename_axis("unit_key")
