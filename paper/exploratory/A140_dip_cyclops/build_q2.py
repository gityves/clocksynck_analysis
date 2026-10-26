"""A140 Q2: pseudobulks of GSE266919 (Zhang et al., Cancer Cell 2025), per sample x compartment.

Reads the five deposited SingleCellExperiment .rds.gz files with the pure-python reader `rds_parse.py` (nothing is
executed, R is not used; run with `python -I`). One object at a time, in two passes: pass 1 reads the cell metadata
(Sample, Patient, Group, Tissue), the gene names and the matrix shape (the big vectors are only summarised); pass 2
streams the dense double count matrix (genes x cells, column major) in chunks of cells and sums the raw UMIs per sample.

Outputs (output/q2/, git-ignored):
  pb_<compartment>.csv.gz   sample x gene raw UMI sums (all genes of the object), plus n_cells
  cells_per_sample.csv      sample x compartment cell counts (0 where the object has no cell of the sample)
  samples.csv               Sample, Patient, Group, Tissue (aggregates only; no per-cell rows)

    python -I build_q2.py [--subsample N] [--jobs J]      (objects in parallel; --subsample: first N samples, check only)
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))  # paper/

import numpy as np
import pandas as pd
import scipy.sparse as sp

import rds_parse as R

RAW = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian")) / "tnbc_zhang2025" / "raw"
OBJECTS = {"CD4_T": "CD4Tcell", "CD8_T": "CD8Tcell", "myeloid": "Myeloid", "B": "Bcell", "NK": "NKcell"}
CHUNK_CELLS = 1500  # 1500 cells x ~21k genes x 8 B = 250 MB per chunk
META = ["Sample", "Patient", "Group", "Tissue"]


def named(x):
    return dict(zip(x["names"].data, x.data))


def pass1(path: Path):
    """Cell metadata (Sample, Patient, Group, Tissue), genes, shape of the counts."""
    o = R.read_rds(str(path))
    cols = named(o["colData"]["listData"])
    df = {}
    for k in META:
        v = cols[k]
        lv = v["levels"].data if "levels" in v else None
        df[k] = [lv[i - 1] if i > 0 else None for i in v.data] if lv is not None else list(v.data)
    a = named(o["assays"]["data"]["listData"])["counts"]
    dim = [int(x) for x in a["dim"].data]
    genes = list(a["dimnames"].data[0].data)
    return pd.DataFrame(df), genes, dim


def pass2(path: Path, ng: int, nc: int, codes: np.ndarray, nsamp: int, limit: int | None):
    """Stream the counts and sum them per sample (codes: sample index of each cell, -1 = skip)."""
    sums = np.zeros((nsamp, ng))
    state = {"done": 0}
    orig = R.Reader.num

    def num(self, t, n):
        if t != 14 or n != ng * nc:
            return orig(self, t, n)
        done = 0
        while done < nc:
            k = min(CHUNK_CELLS, nc - done)
            X = np.frombuffer(self.rd(8 * ng * k), dtype=">f8").reshape(k, ng)  # cells x genes
            c = codes[done:done + k]
            ok = np.flatnonzero(c >= 0)
            if len(ok):
                ind = sp.csr_matrix((np.ones(len(ok)), (c[ok], ok)), shape=(nsamp, k))
                sums[...] += ind @ X
            done += k
            if (done // CHUNK_CELLS) % 20 == 0:
                print(f"  {path.name}: {done}/{nc} cells", flush=True)
        state["done"] = done
        return R.BigVec("real", n, np.zeros(1), True, 0.0, 0.0)

    R.Reader.num = num
    try:
        R.read_rds(str(path))
    finally:
        R.Reader.num = orig
    assert state["done"] == nc, "count matrix not streamed"
    return sums


def build_one(args):
    comp, fname, out, subsample = args
    path = RAW / f"GSE266919_{fname}.rds.gz"
    meta, genes, (ng, nc) = pass1(path)
    assert len(meta) == nc and len(genes) == ng
    samples = sorted(meta.Sample.dropna().unique())
    if subsample:
        samples = samples[:subsample]
    idx = {s: i for i, s in enumerate(samples)}
    codes = meta.Sample.map(idx).fillna(-1).astype(int).to_numpy()
    n_cells = np.bincount(codes[codes >= 0], minlength=len(samples))
    print(f"{comp}: {nc} cells, {ng} genes, {len(samples)} samples", flush=True)
    sums = pass2(path, ng, nc, codes, len(samples), subsample)
    pb = pd.DataFrame(sums, index=pd.Index(samples, name="Sample"), columns=genes)
    pb.insert(0, "n_cells", n_cells)
    assert np.allclose(sums, np.round(sums))
    pb.to_csv(out / f"pb_{comp}.csv.gz")
    sm = meta.dropna(subset=["Sample"]).drop_duplicates("Sample").set_index("Sample")[META[1:]]
    print(f"{comp}: done", flush=True)
    return comp, pd.Series(n_cells, index=samples), sm


def main(argv=None):
    from common import cli
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    cli.setup(args)
    out = HERE / "output" / ("q2_subsample" if args.subsample else "q2")
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(c, f, out, args.subsample) for c, f in OBJECTS.items()]
    if args.jobs > 1:
        with ProcessPoolExecutor(min(args.jobs, len(jobs))) as ex:
            res = list(ex.map(build_one, jobs))
    else:
        res = [build_one(j) for j in jobs]
    counts = pd.DataFrame({c: n for c, n, _ in res}).fillna(0).astype(int)
    counts.index.name = "Sample"
    counts.to_csv(out / "cells_per_sample.csv")
    sm = pd.concat([s for _, _, s in res]).groupby(level=0).first().loc[counts.index]
    sm.to_csv(out / "samples.csv")
    print(f"samples {len(sm)}; cells per compartment {counts.sum().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
