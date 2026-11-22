"""Loaders of the Jed deliveries with the derived tables of `derive.py` attached.

Unspliced counts: use kind="intron_direct" (METHOD.md section 9). kind="intron" (gene body minus exon) loses clock
genes that overlap another gene in unstranded libraries (CLOCK, PER2 zero in GSE326530); `bulk` warns when asked for it
on a delivery with unstranded runs.

Counts stay raw (clocksync reads counts with their depth); `norm="cpm"` / `"logcpm"` divides by the depth of the same
count type and is meant for QC and plots, not for clocksync input.

    from common.mapping import load
    counts, sheet = load.bulk("cardosomoreira2019_postnatal", "exon")            # library x symbol, one row per library
    cells = load.cells("kameneva2021_fetal_adrenal")                            # one row per delivered barcode
    x = load.adata("kameneva2021_fetal_adrenal", "GSM4446535")                  # AnnData, obs = the cells table
    counts, sheet = load.pseudobulk("kameneva2021_fetal_adrenal", "unspliced")  # donor|class x symbol, A244 cells
"""
import warnings

import anndata as ad
import numpy as np
import pandas as pd

from clocksync.genes import canonical
from common.paths import DATA

SRC, DEST = DATA / "sccc_mapping", DATA / "sccc_mapping_derived"
KINDS = ("exon", "intron", "intron_direct", "genebody")


def _norm(counts: pd.DataFrame, depth: pd.Series, norm: str | None) -> pd.DataFrame:
    if norm is None:
        return counts
    cpm = counts.div(depth.to_numpy(), axis=0) * 1e6
    return np.log1p(cpm) if norm == "logcpm" else cpm


def sheet(name: str) -> pd.DataFrame:
    return pd.read_csv(DEST / name / "sheet.csv")


def bulk(name: str, kind: str = "exon", norm: str | None = None, used_only: bool = True):
    """(counts, sheet): runs summed per library, genes as canonical symbols (duplicates summed); libraries the
    authors did not use are dropped when `used_only` and the sheet has a `used` column."""
    assert kind in KINDS
    if kind == "intron":
        q = pd.read_csv(SRC / name / "qc.tsv", sep="\t")
        n0 = int((q.strand_used.astype(str) == "0").sum())
        if n0:
            warnings.warn(f"{name}: {n0} unstranded runs; counts_intron loses genes overlapping another gene "
                          "(CLOCK, PER2 ...). Use kind='intron_direct' (METHOD.md section 9).", stacklevel=2)
    s = sheet(name)
    if used_only and "used" in s:
        s = s[s.used.astype(bool)]
    c = pd.read_csv(SRC / name / f"counts_{kind}.tsv.gz", sep="\t", index_col=0)
    sym = c.pop("symbol").astype(str).map(canonical).to_numpy() if "symbol" in c else c.index
    c = c[s.run].groupby(sym).sum().T.groupby(s.library.to_numpy()).sum()
    lib = s.groupby("library", sort=False).agg(
        {**{k: "first" for k in s.columns if not k.startswith("depth_") and k not in ("run", "library")},
         **{k: "sum" for k in s.columns if k.startswith("depth_")}, "run": " ".join})
    lib = lib.loc[c.index]
    return _norm(c, lib[f"depth_{kind}"], norm), lib


def cells(name: str) -> pd.DataFrame:
    return pd.read_parquet(DEST / name / "cells.parquet")


def adata(name: str, gsm: str) -> ad.AnnData:
    """AnnData of one sample (by its corrected GSM: the file is taken from the cells table)."""
    c = cells(name)
    c = c[c.gsm == gsm]
    x = ad.read_h5ad(SRC / c.file.iloc[0])
    x.obs = c.set_index("barcode").loc[x.obs_names]
    return x


def pseudobulk(name: str, layer: str = "spliced", group: str | None = "a244_class", min_cells: int = 200,
               min_umi: int = 500):
    """(counts, sheet) per donor and `group` over the cells whose `group` is non-empty, layer `X` / `spliced` /
    `unspliced` / `ambiguous`; groups with < `min_cells` cells dropped. `group=None`: one pseudobulk per donor over the
    cells that are not blood (whole organ). Cells below `min_umi` Gene UMIs are left out (A242 floor)."""
    c = cells(name)
    keep = (c[group].astype(str).ne("") if group else ~c.blood.astype(bool)) & (c.total_counts >= min_umi)
    key = c.donor + ("|" + c[group].astype(str) if group else "")
    acc, sym = {}, None
    for f in c.file.unique():
        x = ad.read_h5ad(SRC / f)
        m = (c.file == f).to_numpy() & keep.to_numpy()
        sub = c[m].set_index("barcode")
        X = (x.X if layer == "X" else x.layers[layer])[x.obs_names.get_indexer(sub.index)]
        if sym is None:
            sym = pd.Index([canonical(g) for g in x.var.symbols.astype(str)])
        for k in key[m].unique():
            rows = (key[m] == k).to_numpy()
            acc[k] = acc.get(k, 0) + np.asarray(X[rows].sum(0)).ravel()
    counts = pd.DataFrame(acc, index=sym).groupby(level=0).sum().T
    t = c[keep].assign(id=key[keep]).groupby("id").agg(donor=("donor", "first"), n_cells=("barcode", "size"),
                                                         depth=("total_counts", "sum"), age=("age", "first"))
    if group:
        t["compartment"] = t.index.str.split("|").str[1]
    t = t[t.n_cells >= min_cells]
    return counts.loc[t.index], t
