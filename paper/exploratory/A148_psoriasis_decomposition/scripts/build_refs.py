"""A148 reference step: all-gene pseudobulk of the healthy keratinocyte bins (kc_cycling, kc_d1..d4) and the non-keratinocyte classes.

The A63 bin pseudobulk csv holds the 14 clock genes only, so the all-gene bin sums are rebuilt from the raw healthy single-cell
matrices with the cell-to-bin assignments saved by the A63 build (archive output/healthy_bins_cells_<study>.csv.gz, reused
unchanged). Units: healthy donors of the A63 healthy table (GSE173706 healthy, Census Reynolds, palm-hip, sole-Boldo; GSE162183
excluded as in A63), a unit-bin counted if it has >= 20 cells. Counts are summed over units per bin.
Non-keratinocyte references (fibroblast, endothelial, immune = T + myeloid, melanocyte): nonlesional rows of
psoriasis_sc/pseudobulk_all.h5ad (A58 all-gene pseudobulks, same single-cell sources as the psoriasis bins).

    python -I build_refs.py        Output: output/ref_counts.csv.gz (genes x references, summed counts), output/ref_cells.csv
"""
import gzip, re, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
ARCH = Path("/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/validation/paired_damping/output")
OUT = Path(__file__).resolve().parents[1] / "output"
BINS = ["kc_cycling", "kc_d1", "kc_d2", "kc_d3", "kc_d4"]
MINC = 20


def ensg_map():
    d = {}
    with gzip.open(DATA / "reference/gencode.v46.basic.annotation.gtf.gz", "rt") as fh:
        for ln in fh:
            if ln[0] == "#": continue
            f = ln.split("\t")
            if f[2] != "gene": continue
            d[re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]] = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
    return d


def collapse(m, sym):
    sym = np.asarray(sym); u, inv = np.unique(sym, return_inverse=True)
    if len(u) == len(sym): return m, sym
    g = sp.csr_matrix((np.ones(len(sym), np.float32), (np.arange(len(sym)), inv)), shape=(len(sym), len(u)))
    return (m @ g).tocsr(), u


def sums(tag, cells, get):
    """get(lib) -> (csr cells x genes, symbols, barcodes); returns {bin: Series(gene sums)} and a cell-count table."""
    cells = cells[cells.final.isin(BINS)]
    acc, nc = {}, []
    for lib, sub in cells.groupby("lib"):
        m, sym, bc = get(lib)
        ids = pd.Index([f"{lib}_{b}" for b in bc])
        pos = ids.get_indexer(sub.index)
        assert (pos >= 0).all(), (tag, lib)
        for b, s in sub.groupby("final"):
            if len(s) < MINC: continue
            v = pd.Series(np.asarray(m[ids.get_indexer(s.index)].sum(0)).ravel(), index=sym)
            acc[b] = acc[b].add(v, fill_value=0) if b in acc else v
            nc.append(dict(study=tag, unit=lib, bin=b, n_cells=len(s)))
        del m
    return acc, nc


def main():
    res, ncs = {}, []
    emap = ensg_map()
    # GSE173706 healthy (NS) libraries
    cells = pd.read_csv(ARCH / "healthy_bins_cells_GSE173706_NS.csv.gz", index_col=0)
    files = {f.name.split("_NS-")[1].replace(".csv.gz", ""): f for f in (DATA / "psoriasis_sc/GSE173706").glob("GSM*_NS-*.csv.gz")}
    def get(lib):
        df = pd.read_csv(files[lib], index_col=0)
        sym = np.array([emap.get(i.split(".")[0], i) for i in df.index])
        m = sp.csr_matrix(df.values.T.astype(np.float32)); bc = df.columns.values
        m, sym = collapse(m, sym); return m, sym, bc
    a, n = sums("GSE173706_NS", cells, get); res["GSE173706_NS"] = a; ncs += n; print("GSE173706_NS", flush=True)
    for tag in ("palmhip", "soleboldo", "reynolds"):
        cells = pd.read_csv(ARCH / f"healthy_bins_cells_census_{tag}.csv.gz", index_col=0)
        x = ad.read_h5ad(DATA / f"healthy_skin/census_{tag}/{tag}.h5ad")
        X = x.X.tocsr().astype(np.float32); donors = x.obs.donor_id.astype(str).values; bcs = x.obs_names.values; var = x.var_names.values
        del x
        def get(lib, X=X, donors=donors, bcs=bcs, var=var):
            i = np.flatnonzero(donors == lib); m, sym = collapse(X[i], var); return m, sym, bcs[i]
        a, n = sums("census_" + tag, cells, get); res["census_" + tag] = a; ncs += n; print(tag, flush=True)
        del X
    nc = pd.DataFrame(ncs); nc.to_csv(OUT / "ref_cells.csv", index=False)
    tot = {b: sum((res[s][b] for s in res if b in res[s]), start=pd.Series(dtype=float)) for b in BINS}
    # sum properly with fill_value
    ref = {}
    for b in BINS:
        v = None
        for s in res:
            if b in res[s]: v = res[s][b] if v is None else v.add(res[s][b], fill_value=0)
        ref[b] = v
    R = pd.DataFrame(ref).fillna(0)
    # non-keratinocyte classes
    pa = ad.read_h5ad(DATA / "psoriasis_sc/pseudobulk_all.h5ad")
    o = pa.obs; X = pa.X.tocsr() if sp.issparse(pa.X) else sp.csr_matrix(pa.X)
    grp = {"fibroblast": ["fibroblast"], "endothelial": ["endothelial"], "immune": ["T", "myeloid"], "melanocyte": ["melanocyte"]}
    for k, cs_ in grp.items():
        i = np.flatnonzero((o.status == "nonlesional").values & o.compartment.isin(cs_).values & (o.n_cells >= MINC).values)
        R = R.join(pd.Series(np.asarray(X[i].sum(0)).ravel(), index=pa.var_names, name=k), how="outer")
        nc = pd.concat([nc, pd.DataFrame([dict(study="pseudobulk_all", unit="nonlesional rows", bin=k, n_cells=int(o.n_cells.values[i].sum()))])])
    R = R.fillna(0)
    R = R.T.groupby(level=0).sum().T
    R.to_csv(OUT / "ref_counts.csv.gz"); nc.to_csv(OUT / "ref_cells.csv", index=False)
    print(R.shape, R.sum().round(0).to_dict())
    print(nc.groupby("bin").agg(units=("unit", "nunique"), cells=("n_cells", "sum")))


if __name__ == "__main__":
    main()
