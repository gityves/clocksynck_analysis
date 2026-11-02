"""A189: genome-wide Cao 2020 pseudobulks (organ and class level) for the expression-matched pseudo-clock companions.

Same rows and rules as the `fetal_cao` clock-gene tables (unit = Fetus_id, compartment = Organ or Organ:class with
`build.CLASSMAP`, depth = all-gene count of the cell). The GEO loom is decompressed once (21 GB) into
`$SCCC_DATA/natexp/cao2020/loom/` and deleted at the end (data cap). Slabs of 6,400 cells (all genes) are summed per
group with a sparse indicator product in worker processes and accumulated in the parent. Duplicate gene symbols: the
first row is kept (as the clock-gene builder). Check: the clock-gene columns must equal the stored `fetal_cao` tables.

Writes `$SCCC_DATA/natexp/cao2020/genome/pseudobulk_{organ,class}.h5ad` (rows = fetus x compartment, X = counts,
obs: donor, compartment, n_cells, depth) and PROVENANCE lines.

    python build_genome.py [--jobs 4] [--max-cells N]
"""
import argparse
import hashlib
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[3]), str(HERE.parents[3] / "common" / "datasets")]
import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import data  # noqa: E402
from fetal_cao import build as B  # noqa: E402

OUT = B.ROOT / "genome"
LEVELS = ("organ", "class")
CODES = None  # per-process: {level: group code of every cell}, number of groups, kept gene rows


def init(codes, n_groups, keep):
    global CODES
    limit_threads(1)
    CODES = (codes, n_groups, keep)


def slab(args):
    loom, c0, c1 = args
    codes, n_groups, keep = CODES
    with h5py.File(loom, "r") as f:
        X = f["matrix"][:, c0:c1]
    if (X < 0).any() or (X != np.round(X)).any():
        raise ValueError(f"slab {c0}: counts are not non-negative integers")
    depth = X.sum(0, dtype=np.float64)
    Xt = sp.csr_matrix(X[keep].T.astype(np.float64))
    out = {}
    for lev in LEVELS:
        c = codes[lev][c0:c1]
        M = sp.csr_matrix((np.ones(len(c)), (c, np.arange(len(c)))), shape=(n_groups[lev], len(c)))
        rows = np.unique(c)
        out[lev] = (rows, (M @ Xt)[rows].toarray(), np.bincount(c, minlength=n_groups[lev])[rows],
                    np.bincount(c, weights=depth, minlength=n_groups[lev])[rows])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--max-cells", type=int, default=None)
    a = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    gz = B.ROOT / "raw" / "GSE156793_S3_gene_count.loom.gz"
    loom = B.decompressed(gz, B.ROOT / "loom")
    print(f"decompressed {time.time() - t0:.0f}s", flush=True)
    try:
        obs, sym, n_genes = B.cell_attributes(loom, B.ROOT / "raw" / "GSE156793_S1_metadata_cells.txt.gz")
        canon = pd.Series([canonical(s) for s in sym])
        keep = np.flatnonzero(~canon.duplicated().to_numpy())
        genes = canon.iloc[keep].tolist()
        comp = {"organ": obs.Organ, "class": obs.Organ + ":" + obs.Main_cluster_name.map(B.classify)}
        keys, codes, n_groups = {}, {}, {}
        for lev in LEVELS:
            k = pd.MultiIndex.from_arrays([obs.Fetus_id, comp[lev]])
            codes[lev], keys[lev] = pd.factorize(k)
            n_groups[lev] = len(keys[lev])
        n = len(obs) if a.max_cells is None else min(a.max_cells, len(obs))
        acc = {lev: (np.zeros((n_groups[lev], len(keep))), np.zeros(n_groups[lev]), np.zeros(n_groups[lev]))
               for lev in LEVELS}
        tasks = [(loom, c0, min(c0 + B.SLAB, n)) for c0 in range(0, n, B.SLAB)]
        print(f"{len(obs)} cells, {len(keep)} genes, {len(tasks)} slabs", flush=True)
        with Pool(a.jobs, initializer=init, initargs=(codes, n_groups, keep)) as pool:
            for i, r in enumerate(pool.imap_unordered(slab, tasks)):
                for lev, (rows, X, nc, dp) in r.items():
                    acc[lev][0][rows] += X
                    acc[lev][1][rows] += nc
                    acc[lev][2][rows] += dp
                if (i + 1) % 50 == 0:
                    print(f"{i + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    finally:
        if a.max_cells is None:
            Path(loom).unlink(missing_ok=True)
            print("loom deleted", flush=True)
    for lev in LEVELS:
        X, nc, dp = acc[lev]
        ok = nc > 0
        o = pd.DataFrame({"donor": keys[lev].get_level_values(0), "compartment": keys[lev].get_level_values(1),
                          "n_cells": nc, "depth": dp})[ok].reset_index(drop=True)
        g = ad.AnnData(sp.csr_matrix(X[ok]), obs=o.set_axis(o.index.astype(str)), var=pd.DataFrame(index=genes))
        # check against the stored clock-gene table
        ref, _ = data.load("fetal_cao", level=lev)
        m = o.reset_index().merge(ref, on=["donor", "compartment"], suffixes=("", "_ref"))
        cg = [c for c in ref.columns if c in set(genes) and c not in ("donor", "compartment")]
        if a.max_cells is None:
            col = g.var_names.get_indexer(cg)
            new = g.X[m["index"].to_numpy()][:, col].toarray()
            diff = np.abs(new - m[cg].to_numpy(float)).max()
            ddep = np.abs(m.depth - m.depth_ref).max()
            print(f"{lev}: {len(o)} rows ({len(ref)} stored), {len(cg)} clock genes max |diff| {diff:.3g}, "
                  f"depth max |diff| {ddep:.3g}", flush=True)
        sfx = "" if a.max_cells is None else "_check"
        path = OUT / f"pseudobulk_{lev}{sfx}.h5ad"
        g.write_h5ad(path, compression="gzip")
        h = hashlib.sha256(path.read_bytes()).hexdigest()
        with open(OUT / "PROVENANCE.md", "a") as f:
            f.write(f"- {time.strftime('%Y-%m-%d')} `{path.name}`: {path.stat().st_size} bytes, sha256 {h}; built by "
                    f"`paper/exploratory/A189_fetal_within/scripts/build_genome.py` from `raw/{gz.name}`\n")
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
