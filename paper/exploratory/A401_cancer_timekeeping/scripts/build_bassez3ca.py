"""A401: pseudobulk of Bassez 2021 breast anti-PD-1 (BIOKEY, 3CA copy), streamed from the MatrixMarket file.

Groups = sample x timepoint x compartment (Mast dropped), kept at >= 20 cells. No clock statistic.
Run: python -I build_bassez3ca.py   (reads the mtx extracted next to Cells.csv; delete it afterwards).
"""
import sys, hashlib, time
import numpy as np, pandas as pd

ROOT = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/bassez2021_3ca"
RAW = f"{ROOT}/raw/Data_Bassez2021_Breast"
COMP = {"Malignant": "tumour:malignant", "Fibroblast": "tumour:fibroblast", "Endothelial": "tumour:endothelial",
        "Myeloid": "tumour:myeloid", "Dendritic": "tumour:myeloid", "T_cell": "tumour:T_NK", "B_cell": "tumour:B_plasma"}
MIN_CELLS = 20

genes = pd.read_csv(f"{RAW}/Genes.txt", header=None, names=["g"])["g"].str.strip('"').to_numpy()
cells = pd.read_csv(f"{RAW}/Cells.csv", usecols=["cell_name", "sample", "cell_type", "timepoint", "nCount_RNA", "complexity",
                                                 "expansion", "cohort", "BC_type"], keep_default_na=False)
with open(f"{RAW}/Exp_data_UMIcounts.mtx") as f:
    line = f.readline()
    while line.startswith("%"):
        line = f.readline()
    ng, nc, nnz = map(int, line.split())
    header_lines = 1
    f.seek(0)
    header_lines = 0
    for l in f:
        header_lines += 1
        if not l.startswith("%"):
            break
assert ng == len(genes) and nc == len(cells), (ng, nc, len(genes), len(cells))
print(f"header ok: {ng} genes x {nc} cells, {nnz} entries; skip {header_lines} lines", flush=True)

cells["compartment"] = cells["cell_type"].map(COMP)
keep = cells["compartment"].notna()
cells["grp"] = (cells["sample"] + "|" + cells["timepoint"] + "|" + cells["compartment"]).where(keep)
codes, uniq = pd.factorize(cells["grp"])  # -1 for dropped
ngrp = len(uniq)
cidx = codes.astype(np.int64)
acc = np.zeros(ng * ngrp, dtype=np.float64)
cell_umi = np.zeros(nc)
cell_nf = np.zeros(nc)
t0, done = time.time(), 0
for ch in pd.read_csv(f"{RAW}/Exp_data_UMIcounts.mtx", sep=" ", header=None, skiprows=header_lines,
                      dtype=np.int32, chunksize=20_000_000, engine="c"):
    g = ch[0].to_numpy(np.int64) - 1
    c = ch[1].to_numpy(np.int64) - 1
    v = ch[2].to_numpy()
    cell_umi += np.bincount(c, weights=v, minlength=nc)
    cell_nf += np.bincount(c, minlength=nc)
    gi = cidx[c]
    m = gi >= 0
    acc += np.bincount(g[m] * ngrp + gi[m], weights=v[m], minlength=ng * ngrp)
    done += len(ch)
    print(f"{done/1e6:.0f}M / {nnz/1e6:.0f}M entries, {time.time()-t0:.0f}s", flush=True)
assert done == nnz
# Cells.csv nCount_RNA is NOT the column sum of the matrix (rank corr ~0 with it); `complexity` equals the number of
# detected genes per matrix column in every cell, which confirms the column order. The UMI total is taken from the matrix.
cells["umi"] = cell_umi
print(f"complexity == detected genes per column: {(cells['complexity'].to_numpy() == cell_nf).mean():.4f}; "
      f"nCount_RNA == column sum: {(cells['nCount_RNA'].to_numpy() == cell_umi).mean():.4f} "
      f"(Spearman {pd.Series(cell_umi).corr(cells['nCount_RNA'], method='spearman'):.3f})", flush=True)
mat = acc.reshape(ng, ngrp).T  # groups x genes

info = cells[keep].groupby("grp", sort=False).agg(n_cells=("cell_name", "size"), ncount=("umi", "sum"),
                                                  expansion=("expansion", "first"), cohort=("cohort", "first"),
                                                  BC_type=("BC_type", "first")).reindex(uniq)
info["depth"] = mat.sum(1)
# check: group depth vs sum of per-cell matrix column sums
agree = (info["depth"] == info["ncount"])
print(f"depth == sum of per-cell matrix UMI in {agree.sum()}/{len(info)} groups; max |rel diff| = "
      f"{(abs(info['depth']-info['ncount'])/info['ncount']).max():.3g}; total depth {info['depth'].sum():.0f} "
      f"vs total per-cell UMI {info['ncount'].sum():.0f}", flush=True)
med = cells.groupby(["sample", "timepoint"])["umi"].median()  # all cells of the sample x timepoint
parts = pd.Series(uniq).str.split("|", expand=True)
parts.columns = ["smp", "timepoint", "compartment"]
out = pd.DataFrame({"donor": "bassez2021|" + parts["smp"], "sample": parts["smp"] + "|" + parts["timepoint"],
                    "timepoint": parts["timepoint"], "compartment": parts["compartment"], "stratum": "bassez2021|10x",
                    "n_cells": info["n_cells"].to_numpy(), "depth": info["depth"].to_numpy().astype(np.int64)})
out["median_umi"] = [med[(s, t)] for s, t in zip(parts["smp"], parts["timepoint"])]
for k in ["expansion", "cohort", "BC_type"]:
    out[k] = info[k].to_numpy()
G = pd.DataFrame(mat.astype(np.int64), columns=genes).T.groupby(level=0).sum().T  # sum duplicate symbols
print(f"genes {ng} -> {G.shape[1]} unique symbols", flush=True)
keepg = (out["n_cells"] >= MIN_CELLS).to_numpy()
res = pd.concat([out[keepg].reset_index(drop=True), G[keepg].reset_index(drop=True)], axis=1)
res.to_parquet(f"{ROOT}/pseudobulk.parquet", index=False)
print(f"groups {ngrp} -> {len(res)} with >= {MIN_CELLS} cells", flush=True)

# donors
sm = pd.read_csv(f"{RAW}/Samples.csv")
pt = cells.groupby("sample").agg(timepoints=("timepoint", lambda s: "+".join(sorted(set(s)))),
                                 expansion=("expansion", lambda s: "/".join(sorted(set(s)))),
                                 cohort=("cohort", "first"), BC_type=("BC_type", "first")).reset_index()
pt = pt.rename(columns={"sample": "patient"})
don = pt.merge(sm.drop(columns=["patient"]).rename(columns={"sample": "patient"}), on="patient", how="left")
don.to_csv(f"{ROOT}/donors.csv", index=False)
print(f"donors.csv: {len(don)} patients", flush=True)

# counts
meta = res[["donor", "timepoint", "compartment", "n_cells", "expansion", "cohort"]]
rows = []
for (d, tp), g in meta.groupby(["donor", "timepoint"]):
    n = g.set_index("compartment")["n_cells"]
    mal = n.get("tumour:malignant", 0) >= MIN_CELLS
    f, e = n.get("tumour:fibroblast", 0) >= 50, n.get("tumour:endothelial", 0) >= 50
    rows.append(dict(donor=d, timepoint=tp, ok_any=mal and (f or e), ok_both=mal and f and e,
                     expansion=g["expansion"].iloc[0], cohort=g["cohort"].iloc[0]))
q = pd.DataFrame(rows)
print("\nmalignant (>=20) + fibroblast and/or endothelial (>=50), per timepoint:")
print(q.groupby("timepoint")[["ok_any"]].sum().to_string())
print("malignant + both fibroblast and endothelial (>=50), per timepoint:")
print(q.groupby("timepoint")[["ok_both"]].sum().to_string())
for col in ["ok_any", "ok_both"]:
    w = q.pivot(index="donor", columns="timepoint", values=col).fillna(False)
    both = w.index[w.all(axis=1) & (w.shape[1] == 2)]
    qq = q[q.donor.isin(both)].drop_duplicates("donor")
    print(f"\npatients with {col} at both time points, per expansion x cohort (n={len(both)}):")
    print(qq.groupby(["expansion", "cohort"]).size().to_string())
print("\nsha256", hashlib.sha256(open(f"{ROOT}/pseudobulk.parquet", "rb").read()).hexdigest())
