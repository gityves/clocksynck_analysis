"""A170 data preparation (registered in paper/docs/registered_experiments.md, A170): per-patient pseudobulks (summed raw
counts, all genes) of Wu 2021 (GSE176078) and Pal 2021 (GSE161529) breast tumours, by cell class: malignant, whole
tumour (all cells), and for Wu also CAF and endothelial. No clock statistic is computed here.

    python -I prep_a170.py
Reads downloaded files only (mtx, tsv, csv, .rds parsed by the pure-python `rds_parse`); nothing downloaded is executed.
Writes ../output/: units.csv (unit, cohort, patient, er, cls, n_cells, depth), counts_<cohort>.parquet (units x genes).
"""
import gzip
import io
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rds_parse as R  # noqa: E402

D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sc_breast_er")
SCOUT = HERE.parents[2] / "exploratory/scouts_2026-10-31/sc_breast_er"
OUT = HERE.parent / "output"
OUT.mkdir(exist_ok=True)
WU = D / "GSE176078/extract/Wu_etal_2021_BRCA_scRNASeq"
PAL_ZIP = D / "GSE161529/figshare17058077/HumanBreast10X.zip"
PAL_OBJ = {"ERTotal": ("ER+", {0, 4}), "TNBC": ("ER-", {0, 2})}  # authors' Tum clusters (0-based seurat_clusters)
WU_CLASS = {"Cancer Epithelial": "malignant", "CAFs": "caf", "Endothelial": "endothelial"}
CHUNK = 20_000_000


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def wu():
    meta = pd.read_csv(WU / "metadata.csv", index_col=0)
    er = pd.read_csv(SCOUT / "wu_reconciled.csv").set_index("orig.ident").final_ER
    bc = pd.read_csv(WU / "count_matrix_barcodes.tsv", header=None)[0].to_numpy()
    genes = pd.read_csv(WU / "count_matrix_genes.tsv", header=None)[0].to_numpy()
    meta = meta.reindex(bc)
    assert meta["orig.ident"].notna().all()
    pat = meta["orig.ident"].to_numpy()
    cls = meta.celltype_major.map(WU_CLASS).to_numpy()
    units = [f"Wu|{p}|{c}" for p in sorted(set(pat)) for c in ("whole", *WU_CLASS.values())]
    ix = {u: j for j, u in enumerate(units)}
    u_whole = np.array([ix[f"Wu|{p}|whole"] for p in pat])  # each cell is in its patient's 'whole' unit
    u_cls = np.array([ix.get(f"Wu|{p}|{c}", -1) for p, c in zip(pat, cls)])  # and in its class unit, if any
    U, G = len(units), len(genes)
    acc, n = np.zeros(G * U), 0
    for ch in pd.read_csv(WU / "count_matrix_sparse.mtx", sep=" ", skiprows=2, header=None, chunksize=CHUNK,
                          dtype={0: np.int64, 1: np.int64, 2: np.float64}):
        g, c, v = ch[0].to_numpy() - 1, ch[1].to_numpy() - 1, ch[2].to_numpy()
        acc += np.bincount(g * U + u_whole[c], weights=v, minlength=G * U)
        k = u_cls[c] >= 0
        acc += np.bincount(g[k] * U + u_cls[c][k], weights=v[k], minlength=G * U)
        n += len(ch)
        say(f"Wu: {n:,} entries")
    acc = acc.reshape(G, U)
    ncell = pd.Series(u_whole).value_counts().add(pd.Series(u_cls[u_cls >= 0]).value_counts(), fill_value=0)
    X = pd.DataFrame(acc.T, index=units, columns=genes)
    X = X.T.groupby(level=0).sum().T  # duplicated symbols summed
    rows = [dict(unit=u, cohort="Wu", patient=u.split("|")[1], er=er[u.split("|")[1]], cls=u.split("|")[2],
                 n_cells=int(ncell.get(j, 0)), depth=float(X.loc[u].sum())) for j, u in enumerate(units)]
    return X, pd.DataFrame(rows)


def read10x(z: zipfile.ZipFile, sample: str):
    base = f"Data/{sample}/"
    bc = pd.read_csv(io.BytesIO(gzip.decompress(z.read(base + "barcodes.tsv.gz"))), header=None)[0].to_numpy()
    ft = pd.read_csv(io.BytesIO(gzip.decompress(z.read(base + "features.tsv.gz"))), header=None, sep="\t")
    mt = pd.read_csv(io.BytesIO(gzip.decompress(z.read(base + "matrix.mtx.gz"))), sep=" ", comment="%", header=None)
    shape, mt = mt.iloc[0], mt.iloc[1:]
    assert int(shape[0]) == len(ft) and int(shape[1]) == len(bc)
    return bc, ft[1].to_numpy(), mt[0].to_numpy(int) - 1, mt[1].to_numpy(int) - 1, mt[2].to_numpy(float)


def fac(c):
    """An R factor (or character vector) as strings."""
    if "levels" not in c:
        return np.array([str(x) for x in c.data], dtype=object)
    lv = list(c["levels"].data)
    return np.array([lv[i - 1] for i in np.asarray(c.data)], dtype=object)


def pal():
    z = zipfile.ZipFile(PAL_ZIP)
    rows, mats = [], {}
    for obj, (er, tum) in PAL_OBJ.items():
        md = R.read_rds(D / f"GSE161529/figshare17058077/SeuratObject_{obj}.rds")["meta.data"]
        vals = dict(zip(list(md["names"].data), md.data))
        cl, grp = fac(vals["seurat_clusters"]).astype(int), fac(vals["group"])
        cells = pd.DataFrame({"group": grp, "cluster": cl}, index=list(md["row.names"].data))
        for g, d in cells.groupby("group"):
            sample = g.replace("_", "-")
            bc, sym, gi, ci, v = read10x(z, sample)
            names = np.array([f"{g}_{b}" for b in bc])
            pos = pd.Series(np.arange(len(names)), index=names)
            miss = (~d.index.isin(pos.index)).sum()
            assert miss == 0, (g, miss)
            for c, sel in (("whole", d.index), ("malignant", d.index[d.cluster.isin(tum)])):
                keep = np.zeros(len(bc), bool)
                keep[pos[sel].to_numpy()] = True
                s = np.bincount(gi[keep[ci]], weights=v[keep[ci]], minlength=len(sym))
                u = f"Pal|{g}|{c}"
                mats[u] = pd.Series(s, index=sym).groupby(level=0).sum()
                rows.append(dict(unit=u, cohort="Pal", patient=g, er=er, cls=c, n_cells=len(sel), depth=float(s.sum())))
            say(f"Pal {obj} {g}: {len(d)} cells, malignant {int(d.cluster.isin(tum).sum())}")
    return pd.DataFrame(mats).T, pd.DataFrame(rows)


if __name__ == "__main__":
    Xp, up = pal()
    Xp.astype(np.float32).to_parquet(OUT / "counts_Pal.parquet")
    Xw, uw = wu()
    Xw.astype(np.float32).to_parquet(OUT / "counts_Wu.parquet")
    u = pd.concat([uw, up], ignore_index=True)
    u.to_csv(OUT / "units.csv", index=False)
    say(u.groupby(["cohort", "cls", "er"]).size().to_string())
