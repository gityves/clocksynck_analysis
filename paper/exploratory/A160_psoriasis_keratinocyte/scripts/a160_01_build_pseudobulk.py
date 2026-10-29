"""A160 step 1: all-gene pseudobulks of GSE173706 (PN, PP, healthy NS libraries) per library x group, from the processed
dense count matrices on disk (parsed only; run with `python -I`).

Cell sets are the archived annotations, reused unchanged (`legacy/.../paired_damping/output/`):
- `psoriasis_sc_cells_GSE173706.csv.gz` (QC-passing cells of PN/PP libraries, class; whole-skin = all of them),
- `healthy_bins_cells_nl_GSE173706.csv.gz` (A63 bin of every PN/PP keratinocyte),
- `healthy_bins_cells_GSE173706_NS.csv.gz` (class and A63 bin of the healthy NS libraries).
Groups per library: whole (all QC cells), kc (all keratinocytes), kc_cycling, kc_d1..kc_d4 (A63 bins). Raw count sums by
gene symbol (GENCODE v46). Output: $SCCC_DATA/psoriasis_sc/a160_pseudobulk.h5ad and output/a160_cells_per_library.csv.

    python -I a160_01_build_pseudobulk.py [--procs 4]
"""
import argparse, gzip, re, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
ROOT = DATA / "psoriasis_sc" / "GSE173706"
GTF = DATA / "reference" / "gencode.v46.basic.annotation.gtf.gz"
ARCH = Path("/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/validation/paired_damping/output")
OUT = Path(__file__).resolve().parents[1] / "output"
GROUPS = ["whole", "kc", "kc_cycling", "kc_d1", "kc_d2", "kc_d3", "kc_d4"]


def ensg_map():
    d = {}
    with gzip.open(GTF, "rt") as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            f = ln.split("\t")
            if f[2] == "gene":
                d[re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]] = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
    return d


def libraries():
    rows = []
    for f in sorted(ROOT.glob("*.csv.gz")):
        m = re.match(r"GSM\d+_(PN|PP|NS)-(\w+)\.csv\.gz", f.name)
        if not m:
            continue
        pid = m.group(2)
        visit = "V3" if (pid.endswith("V3") and m.group(1) != "NS") else "V1"
        pid = pid[:-2] if visit == "V3" else pid
        rows.append(dict(path=f, patient=pid, visit=visit, status={"PP": "lesional", "PN": "nonlesional", "NS": "healthy"}[m.group(1)]))
    t = pd.DataFrame(rows)
    ps = t[t.status != "healthy"]
    both = ps.groupby(["patient", "visit"]).status.nunique()
    ps = ps[[(p, v) in set(both[both == 2].index) for p, v in zip(ps.patient, ps.visit)]].reset_index(drop=True)
    ps["lib"] = ps.index.astype(str)  # as the legacy build (index of the paired PN/PP table)
    h = t[t.status == "healthy"].copy()
    h["lib"] = h.patient
    return pd.concat([ps, h], ignore_index=True)


def cells_table():
    c = pd.read_csv(ARCH / "psoriasis_sc_cells_GSE173706.csv.gz", index_col=0)
    b = pd.read_csv(ARCH / "healthy_bins_cells_nl_GSE173706.csv.gz", index_col=0)
    c["group_bin"] = b["final"].reindex(c.index)
    h = pd.read_csv(ARCH / "healthy_bins_cells_GSE173706_NS.csv.gz", index_col=0)
    h["group_bin"] = h["final"].where(h["final"].str.startswith("kc_"))
    h["status"] = "healthy"
    h["visit"] = "V1"
    return c, h


def one(args):
    r, ids_bin = args
    emap = ENSG
    ids, comp, gbin = ids_bin
    pos = {}
    out = {}
    first = True
    sums = None
    genes = []
    cols = None
    ncell = len(ids)
    Gm = np.zeros((ncell, len(GROUPS)), np.float32)
    Gm[:, 0] = 1
    Gm[:, 1] = (comp == "keratinocyte")
    for j, g in enumerate(GROUPS[2:], 2):
        Gm[:, j] = (gbin == g)
    parts = []
    for ch in pd.read_csv(r["path"], index_col=0, chunksize=3000):
        if cols is None:
            bc = [f"{r['lib']}_{x}" for x in ch.columns]
            ix = pd.Index(ids)
            sel = np.flatnonzero(np.isin(bc, ids))
            order = pd.Index(bc)[sel]
            perm = ix.get_indexer(order)
            assert len(sel) == len(ids), (r["path"].name, len(sel), len(ids))
            cols = sel
        x = ch.values[:, cols].astype(np.float32)
        parts.append((np.array([emap.get(i.split(".")[0], i) for i in ch.index]), x @ Gm[perm]))
    sym = np.concatenate([p[0] for p in parts])
    M = np.vstack([p[1] for p in parts])
    u, inv = np.unique(sym, return_inverse=True)
    S = np.zeros((len(u), len(GROUPS)), np.float64)
    np.add.at(S, inv, M)
    return r["lib"], u, S.T.astype(np.float32), Gm.sum(0)


ENSG = None


def main():
    global ENSG
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    ENSG = ensg_map()
    libs = libraries()
    c, h = cells_table()
    tasks = []
    for r in libs.to_dict("records"):
        d = h[h.lib.astype(str) == r["lib"]] if r["status"] == "healthy" else c[c.lib.astype(str) == r["lib"]]
        tasks.append((r, (d.index.to_numpy(), d.compartment.to_numpy(), d.group_bin.to_numpy())))
        print(r["status"], r["patient"], r["visit"], r["lib"], len(d), flush=True)
    with Pool(a.procs) as p:
        res = p.map(one, tasks, chunksize=1)
    genes = sorted(set().union(*[set(x[1]) for x in res]))
    gi = {g: i for i, g in enumerate(genes)}
    rows, X, meta = [], [], []
    for r, (lib, u, S, n) in zip(libs.to_dict("records"), res):
        idx = np.array([gi[g] for g in u])
        for k, g in enumerate(GROUPS):
            v = np.zeros(len(genes), np.float32)
            v[idx] = S[k]
            X.append(v)
            meta.append(dict(lib=lib, patient=r["patient"], visit=r["visit"], status=r["status"], group=g, n_cells=int(n[k]), depth=float(v.sum())))
    meta = pd.DataFrame(meta)
    meta.index = meta.lib + "|" + meta.group
    A = ad.AnnData(sp.csr_matrix(np.vstack(X)), obs=meta, var=pd.DataFrame(index=genes))
    A.write_h5ad(DATA / "psoriasis_sc" / "a160_pseudobulk.h5ad")
    meta.drop(columns=["depth"]).to_csv(OUT / "a160_cells_per_library.csv", index=False)
    print(meta[meta.group.isin(["kc", "whole"])].to_string(), flush=True)


if __name__ == "__main__":
    main()
