"""Scout (design only): per-cell UMI of the 14 clock genes (summed), MGMT and total UMI from the Extended GBmap raw counts.

Reads raw.X in row chunks, keeps only the columns of interest. No gene-level clock statistic is computed.
Output: output/gbmap_cell_counts.csv.gz (cell id, clock UMI sum, MGMT UMI, total UMI).
"""
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from clocksync.config import limit_threads
from clocksync.genes import CLOCK_GENES, gene_columns

limit_threads()
SRC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_gbmap/extended_gbmap_230c8701.h5ad")
OUT = Path(__file__).resolve().parents[1] / "output" / "gbmap_cell_counts.csv.gz"
CHUNK = 40_000

with h5py.File(SRC, "r") as f:
    var = read_elem(f["raw/var"])
    names = var["feature_name"].astype(str)
    clock = gene_columns(names, CLOCK_GENES, missing="warn")
    mg = gene_columns(names, ["MGMT"], missing="raise")[0]
    print("clock genes found:", int((clock >= 0).sum()), "of", len(CLOCK_GENES), "MGMT col", mg, flush=True)
    clock = clock[clock >= 0]
    lut = np.zeros(len(names), dtype=np.int8)
    lut[clock] = 1
    lut[mg] = 2
    obs_ids = read_elem(f["obs"]["_index"]) if "_index" in f["obs"] else None
    X = f["raw/X"]
    indptr = X["indptr"][:]
    n = len(indptr) - 1
    res = np.zeros((n, 3), dtype=np.float64)  # clock sum, MGMT, total
    t0 = time.time()
    for a in range(0, n, CHUNK):
        b = min(a + CHUNK, n)
        s, e = indptr[a], indptr[b]
        d = X["data"][s:e].astype(np.float64)
        ix = X["indices"][s:e]
        rows = np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))
        res[a:b, 2] = np.bincount(rows, weights=d, minlength=b - a)
        k = lut[ix]
        m = k == 1
        res[a:b, 0] = np.bincount(rows[m], weights=d[m], minlength=b - a)
        m = k == 2
        res[a:b, 1] = np.bincount(rows[m], weights=d[m], minlength=b - a)
        print(f"{b}/{n} cells, {time.time() - t0:.0f}s", flush=True)

df = pd.DataFrame(res, columns=["clock_umi", "mgmt_umi", "total_umi"], index=obs_ids)
OUT.parent.mkdir(exist_ok=True)
df.to_csv(OUT)
print("written", OUT, df.shape)
