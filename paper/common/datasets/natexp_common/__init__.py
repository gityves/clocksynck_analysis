"""Helpers shared by the builders of the natural-experiment datasets (COVID-19, placenta, fetal atlases, mouse embryo).

Every builder returns `(pb, donors)`: the standard pseudobulk table of `clocksync` (donor, compartment, [stratum],
n_cells, depth, the 14 clock genes) and a donor table indexed by donor. Sources that `clocksync.pseudobulk` cannot
read directly (one csv per sample, 10x h5, loom, MatrixMarket streams, summed partials) are read in small chunks
into AnnData objects of the 14 clock genes with the genome-wide count of every cell in `obs["depth"]`; these feed the
same aggregation (`clocksync.pseudobulk(..., depth="depth")`), so all datasets are summed by one code path.

Builders keep every pseudobulk with at least one cell: thresholds (min cells per pseudobulk, min donors per
compartment) belong to the analyses that are pre-registered with them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from clocksync.genes import CLOCK_GENES, canonical, gene_columns

GENES = list(CLOCK_GENES)
PB_COLUMNS = ["donor", "compartment", "stratum", "n_cells", "depth"] + GENES


def standard(pb: pd.DataFrame, extra=()) -> pd.DataFrame:
    """Standard column order (optional columns only if present), sorted by donor and compartment."""
    cols = [c for c in PB_COLUMNS if c in pb] + [c for c in extra if c in pb]
    return pb[cols].sort_values(["donor", "compartment"]).reset_index(drop=True)


def clock_adata(X, symbols, depth, obs: pd.DataFrame, genes=GENES) -> ad.AnnData:
    """Cells x `genes` AnnData from a cells x all-genes count matrix (dense or sparse) and gene symbols; `depth` is
    each cell's genome-wide count. Raises if a gene is absent."""
    cols = gene_columns(symbols, genes)
    if (cols < 0).any():
        raise KeyError(f"genes not found: {[g for g, c in zip(genes, cols) if c < 0]}")
    X = sp.csr_matrix(X)[:, cols] if sp.issparse(X) else sp.csr_matrix(np.asarray(X)[:, cols])
    obs = obs.copy()
    obs["depth"] = np.asarray(depth, float)
    obs.index = obs.index.astype(str)
    return ad.AnnData(X.astype(np.float64), obs=obs, var=pd.DataFrame(index=[canonical(g) for g in genes]))


def read_legacy(path, donor="donor_id", compartment="cell_type", extra=()) -> pd.DataFrame:
    """Convert a legacy `pseudobulk.h5ad` (X = summed counts of a 68-gene panel; obs donor_id, cell_type, n_cells,
    depth) into the standard table. `extra`: further obs columns kept."""
    a = ad.read_h5ad(path)
    cols = gene_columns(a.var_names, GENES, missing="raise")
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    t = pd.DataFrame({"donor": a.obs[donor].astype(str).to_numpy(),
                      "compartment": a.obs[compartment].astype(str).to_numpy(),
                      "n_cells": a.obs["n_cells"].astype(int).to_numpy(), "depth": a.obs["depth"].to_numpy(float)})
    for c in extra:
        t[c] = a.obs[c].to_numpy()
    return pd.concat([t, pd.DataFrame(X[:, cols].astype(float), columns=GENES)], axis=1)


def majority_stratum(pb: pd.DataFrame, assay: str = "assay", vote=None) -> pd.DataFrame:
    """Legacy rule for sources with several assays: each donor keeps the rows of its majority assay (most cells over
    the rows in `vote`, default all); the assay becomes the stratum."""
    v = pb if vote is None else pb[np.asarray(vote)]
    major = (v.groupby(["donor", assay]).n_cells.sum().reset_index().sort_values("n_cells")
             .drop_duplicates("donor", keep="last").set_index("donor")[assay])
    out = pb[pb[assay].to_numpy() == major.reindex(pb.donor).to_numpy()].copy()
    out["stratum"] = out.pop(assay)
    return out


def sum_tables(tables, keys=("donor", "compartment")) -> pd.DataFrame:
    """Add up pseudobulk tables built on disjoint pieces of the same cells (counts, depth and cell numbers)."""
    t = pd.concat(tables, ignore_index=True)
    num = [c for c in t.columns if c not in keys and pd.api.types.is_numeric_dtype(t[c])]
    other = [c for c in t.columns if c not in keys and c not in num]
    out = t.groupby(list(keys), as_index=False)[num].sum()
    if other:
        out = out.merge(t.drop_duplicates(list(keys))[list(keys) + other], on=list(keys))
    return out


def compare(a: pd.DataFrame, b: pd.DataFrame, donors=None, keys=("donor", "compartment"), tol=0.5) -> dict:
    """Row-by-row agreement of two standard tables (restricted to `donors`): same rows, same counts and depths."""
    if donors is not None:
        a, b = a[a.donor.isin(donors)], b[b.donor.isin(donors)]
    m = a.merge(b, on=list(keys), how="outer", suffixes=("_a", "_b"), indicator=True)
    both = m[m._merge == "both"]
    cols = ["n_cells", "depth"] + GENES
    dmax = max(float(np.abs(both[f"{c}_a"] - both[f"{c}_b"]).max()) for c in cols) if len(both) else np.nan
    rel = float(max((np.abs(both[f"{c}_a"] - both[f"{c}_b"]) / both[f"{c}_b"].clip(lower=1)).max() for c in cols)) if len(both) else np.nan
    return dict(rows_a=len(a), rows_b=len(b), only_a=int((m._merge == "left_only").sum()),
                only_b=int((m._merge == "right_only").sum()), max_abs_diff=dmax, max_rel_diff=rel,
                ok=bool(len(both) and not (m._merge != "both").any() and dmax <= tol))
