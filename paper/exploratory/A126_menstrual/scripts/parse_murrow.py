"""A126: Murrow GSE198732 RM raw counts (dgCMatrix .rds.gz) -> donor x cell-type pseudobulks over the panel genes. Python only (repo's
`stream_rds.py` R-serialization reader, copied unchanged), nothing executed. Run with python -I.
    gzip -dc RM.raw.counts.matrix.rds.gz | python -I parse_murrow.py SPILL_DIR            (step 1: spill i/x/p, genes, cells)
    python -I parse_murrow.py SPILL_DIR --pseudobulk META_RDS_GZ                           (step 2: pseudobulks)
"""
import argparse, gzip, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import stream_rds as S
OUT = Path(__file__).resolve().parents[1] / "output"


class MR(S.Reader):
    def spill_target(self, t):
        return {("i", 13): "counts_i.int32", ("x", 14): "counts_x.float32"}.get((self.path[-1] if self.path else "", t))


def step1(spill):
    spill = Path(spill); spill.mkdir(parents=True, exist_ok=True)
    rd = MR(sys.stdin.buffer, spill, 0); rd.header(); obj = rd.item()
    (spill / "structure.txt").write_text("\n".join(S.tree(obj)) + "\n")
    dim = obj.a("Dim").val; dn = obj.a("Dimnames").val
    genes, cells = (list(S.compact_or_vector(d)) for d in dn)
    np.save(spill / "counts_p.npy", np.asarray(obj.a("p").val, dtype=np.int64))
    (spill / "genes.txt").write_text("\n".join(genes) + "\n"); (spill / "cells.txt").write_text("\n".join(cells) + "\n")
    (spill / "x_check.json").write_text(json.dumps(dict(rd.check, dim=[int(x) for x in dim], spilled=rd.spilled), indent=1))
    print(dim, rd.check, flush=True)


def meta(path):
    f = gzip.open(path, "rb")
    if f.peek(2)[:2] == b"\x1f\x8b": f = gzip.open(f, "rb")  # GEO file is gzipped twice
    rd = S.Reader(f, Path("."), 0); rd.header(); o = rd.item()
    cols = {k: S.column(v) for k, v in zip(S.names_of(o), o.val)}
    rn = S.compact_or_vector(o.a("row.names"))
    n = len(next(iter(cols.values())))
    return pd.DataFrame(cols, index=pd.Index(rn if len(rn) == n else range(n), name="cell"))


def step2(spill, meta_path):
    spill = Path(spill)
    m = meta(meta_path)
    genes = (spill / "genes.txt").read_text().split("\n")[:-1]; cells = (spill / "cells.txt").read_text().split("\n")[:-1]
    p = np.load(spill / "counts_p.npy"); I = np.memmap(spill / "counts_i.int32", dtype="<i4", mode="r"); X = np.memmap(spill / "counts_x.float32", dtype="<f4", mode="r")
    print("meta", m.shape, "cells", len(cells), "same names", bool((m.index.astype(str) == np.array(cells)).all()) if len(m) == len(cells) else "n/a", flush=True)
    assert len(m) == len(cells)
    panel = (OUT / "panel_genes.txt").read_text().split()
    sys.path.insert(0, str(Path(__file__).parent)); from a126_lib import canonical
    pos = {}
    for i, g in enumerate(genes): pos.setdefault(canonical(g), i)
    pcols = {g: pos[canonical(g)] for g in panel if canonical(g) in pos}
    print("panel genes found", len(pcols), "of", len(panel), flush=True)
    glut = np.full(len(genes), -1); glut[list(pcols.values())] = np.arange(len(pcols))
    grp = (m.Sample.astype(str) + "|" + m.Type.astype(str)).astype("category"); gi = grp.cat.codes.to_numpy()
    M = np.zeros((len(grp.cat.categories), len(pcols))); depth = np.zeros(len(grp.cat.categories))
    cell_depth = np.zeros(len(cells)); nc = len(cells)
    for s in range(0, nc, 5000):
        a, b = p[s], p[min(s + 5000, nc)]
        ii, xx = np.asarray(I[a:b]), np.asarray(X[a:b], np.float64)
        col = np.repeat(np.arange(s, min(s + 5000, nc)), np.diff(p[s:min(s + 5000, nc) + 1]))
        np.add.at(cell_depth, col, xx)
        k = glut[ii]; ok = k >= 0
        np.add.at(M, (gi[col[ok]], k[ok]), xx[ok])
        if (s // 5000) % 5 == 0: print("cells", s, flush=True)
    np.add.at(depth, gi, cell_depth)
    info = m.assign(grp=grp.to_numpy()).groupby("grp", observed=True).agg(donor=("Sample", "first"), type=("Type", "first"), n_cells=("Type", "size"),
        sort=("Sort", lambda s: "|".join(sorted(set(s)))), age=("Age", "first"), bmi=("BMI", "first"), obesity=("Obesity", "first"), parity=("Parity", "first"),
        hc=("HC_use", "first"), batch=("Batch", lambda s: "|".join(sorted(set(s))))).reindex(grp.cat.categories)
    info["depth"] = depth; info["depth_nCount"] = m.groupby(grp.to_numpy(), observed=True).nCount_RNA.sum().reindex(grp.cat.categories).to_numpy()
    info.to_csv(OUT / "murrow_pb_groups.csv")
    pd.DataFrame(M, index=grp.cat.categories, columns=list(pcols)).astype(np.float32).to_parquet(OUT / "murrow_pb_counts.parquet")
    print("depth vs nCount_RNA ratio", float((info.depth / info.depth_nCount).median()), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("spill"); ap.add_argument("--pseudobulk"); a = ap.parse_args()
    step2(a.spill, a.pseudobulk) if a.pseudobulk else step1(a.spill)
