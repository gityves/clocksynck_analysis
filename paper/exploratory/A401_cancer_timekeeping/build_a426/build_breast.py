"""A426 R1: Wu 2021 (GSE176078) and Pal 2021 (GSE161529) breast tumours with the immune compartments.

The A388 genome tables of the two atlases hold the malignant and stromal classes only. Here every cell class is built from the
counts (same files, same parsing as `a388.build_wu` / `build_pal`):
  Wu   author labels (`celltype_major`): Cancer Epithelial -> malignant, CAFs -> fibroblast, Endothelial, PVL -> mural,
       T-cells -> T_NK, B-cells and Plasmablasts -> B_plasma, Myeloid; Normal Epithelial -> epithelial
  Pal  authors' tumour clusters -> malignant (as A388); every other cell by the `cancer_common.compartments` marker rule (the
       Seurat objects carry no cell-type labels)
Sample median UMI = median nCount (Wu: per patient) or median total UMI of the sample group's cells (Pal). Counts only.
Output `$SCCC_DATA/a426_r1/{GSE176078,GSE161529}/pseudobulk.parquet`.

    build_breast.py wu | pal
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import a426_common as C

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER / "exploratory/A170_sc_breast_er/scripts"))
sys.path.insert(0, str(PAPER / "common/datasets"))
import prep_a170 as P  # noqa: E402
from cancer_common.compartments import marker_compartments  # noqa: E402

WU_MAP = {"Cancer Epithelial": "malignant", "CAFs": "fibroblast", "Endothelial": "endothelial", "PVL": "mural", "T-cells": "T_NK",
          "B-cells": "B_plasma", "Plasmablasts": "B_plasma", "Myeloid": "myeloid", "Normal Epithelial": "epithelial"}


def frame(rows, source, stratum):
    """rows: [(patient, class, n_cells, sample_median_umi, Series of counts)] -> A426 pseudobulk frame."""
    X = pd.DataFrame({i: r[4] for i, r in enumerate(rows)}).T.fillna(0.0)
    X.columns = [str(g) for g in X.columns]
    X = X.T.groupby(level=0).sum().T.reset_index(drop=True).astype(np.float32)
    meta = pd.DataFrame({"donor": [f"{source}|{r[0]}" for r in rows], "sample": [r[0] for r in rows],
                         "compartment": [f"tumour:{r[1]}" for r in rows], "stratum": stratum, "n_cells": [r[2] for r in rows],
                         "depth": X.sum(axis=1).to_numpy(), "median_umi": np.nan, "units": "counts", "sample_median_umi": [r[3] for r in rows]})
    out = pd.concat([meta, X], axis=1)
    return out[out.n_cells >= 20].reset_index(drop=True)


def wu():
    meta = pd.read_csv(P.WU / "metadata.csv", index_col=0)
    bc = pd.read_csv(P.WU / "count_matrix_barcodes.tsv", header=None)[0].to_numpy()
    genes = pd.read_csv(P.WU / "count_matrix_genes.tsv", header=None)[0].to_numpy()
    meta = meta.reindex(bc)
    pat = meta["orig.ident"].to_numpy()
    cls = meta.celltype_major.map(WU_MAP).to_numpy()
    classes = sorted(set(WU_MAP.values()))
    units = [(p, c) for p in sorted(set(pat)) for c in classes]
    ix = {u: j for j, u in enumerate(units)}
    u_of = np.array([ix.get((p, c), -1) for p, c in zip(pat, cls)])
    U, G = len(units), len(genes)
    acc, n = np.zeros(G * U), 0
    for ch in pd.read_csv(P.WU / "count_matrix_sparse.mtx", sep=" ", skiprows=2, header=None, chunksize=20_000_000,
                          dtype={0: np.int64, 1: np.int64, 2: np.float64}):
        g, c, v = ch[0].to_numpy() - 1, ch[1].to_numpy() - 1, ch[2].to_numpy()
        k = u_of[c] >= 0
        acc += np.bincount(g[k] * U + u_of[c][k], weights=v[k], minlength=G * U)
        n += len(ch)
        C.log(f"Wu: {n:,} entries")
    acc = acc.reshape(G, U)
    ncell = pd.Series(u_of[u_of >= 0]).value_counts()
    med = meta.groupby("orig.ident").nCount_RNA.median()
    rows = [(p, c, int(ncell.get(j, 0)), float(med[p]), pd.Series(acc[:, j], index=genes)) for j, (p, c) in enumerate(units)]
    return frame(rows, "GSE176078", "GSE176078|Wu")


def pal():
    z = zipfile.ZipFile(P.PAL_ZIP)
    rows = []
    for obj, (er, tum) in P.PAL_OBJ.items():
        md = P.R.read_rds(P.D / f"GSE161529/figshare17058077/SeuratObject_{obj}.rds")["meta.data"]
        vals = dict(zip(list(md["names"].data), md.data))
        cells = pd.DataFrame({"group": P.fac(vals["group"]), "cluster": P.fac(vals["seurat_clusters"]).astype(int)}, index=list(md["row.names"].data))
        for g, d in cells.groupby("group"):
            bc, sym, gi, ci, v = P.read10x(z, g.replace("_", "-"))
            pos = pd.Series(np.arange(len(bc)), index=[f"{g}_{b}" for b in bc])[d.index].to_numpy()
            X = sp.csr_matrix((v, (ci, gi)), shape=(len(bc), len(sym)))[pos]
            umi = np.asarray(X.sum(1)).ravel()
            comp = np.asarray(marker_compartments(X, sym, umi)).astype(object)
            comp[d.cluster.isin(tum).to_numpy()] = "malignant"
            for k in sorted(set(comp) - {"unassigned"}):
                i = np.flatnonzero(comp == k)
                rows.append((g, k, len(i), float(np.median(umi)), pd.Series(np.asarray(X[i].sum(0)).ravel(), index=sym)))
            C.log(f"Pal {obj} {g}: {len(d)} cells")
    return frame(rows, "GSE161529", "GSE161529|Pal")


if __name__ == "__main__":
    which = sys.argv[1]
    pb = wu() if which == "wu" else pal()
    src = "GSE176078" if which == "wu" else "GSE161529"
    d = C.DATA / "a426_r1" / src
    d.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(d / "pseudobulk.parquet", index=False)
    pd.DataFrame().to_csv(d / "samples.csv", index=False)
    C.log(f"wrote {len(pb)} rows, {pb.donor.nunique()} patients to {d}")
