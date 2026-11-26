"""Shared helpers for the A401 GEO tumour pseudobulk builds (counts -> donor x compartment tables). No clock statistic."""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(ROOT / "paper" / "common" / "datasets"), str(ROOT / "clocksync" / "src")]
from cancer_common.compartments import marker_compartments  # noqa: E402

OUT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo_tumour_a401")
MIN_CELLS = 20
MARKER_MAP = {"epithelial": "epithelial", "fibroblast": "fibroblast", "endothelial": "endothelial", "mural": "mural",
              "myeloid": "myeloid", "T_NK": "T_NK", "B_plasma": "B_plasma"}  # mast / unassigned dropped
STROMA = ["fibroblast", "endothelial", "mural"]


def read_10x(prefix: Path | str, gz: bool = True):
    """Read `<prefix>matrix.mtx(.gz)`, features/genes, barcodes. Returns (cells x genes csr, symbols, barcodes)."""
    prefix = str(prefix)
    if Path(prefix).is_dir():
        prefix += "/"
    sfx = ".gz" if gz else ""
    X = scipy.io.mmread(f"{prefix}matrix.mtx{sfx}").tocsr().T.tocsr()
    for name in ("features.tsv", "genes.tsv"):
        p = Path(f"{prefix}{name}{sfx}")
        if p.exists():
            f = pd.read_csv(p, sep="\t", header=None)
            break
    sym = f.iloc[:, 1 if f.shape[1] > 1 else 0].astype(str).to_numpy()
    bc = pd.read_csv(f"{prefix}barcodes.tsv{sfx}", sep="\t", header=None)[0].astype(str).to_numpy()
    return X, sym, bc


def marker_comp(X, sym, umi=None):
    umi = np.asarray(X.sum(1)).ravel() if umi is None else umi
    return marker_compartments(X, sym, umi)


def label(comp, tissue):
    """Compartment strings `tumour:<c>` or `adjacent_normal:<c>`; tissue is 'tumour' / 'adjacent_normal' per cell."""
    return np.array([f"{t}:{c}" for t, c in zip(tissue, comp)])


def pseudobulk(X, sym, obs: pd.DataFrame, units: str, min_cells: int = MIN_CELLS) -> pd.DataFrame:
    """obs columns: donor, sample, compartment, stratum. X cells x genes (counts or normalised); duplicates summed.
    Returns one row per (donor, sample, compartment) with >= min_cells cells."""
    X = sp.csr_matrix(X)
    obs = obs.reset_index(drop=True)
    key = obs[["donor", "sample", "compartment", "stratum"]].astype(str).agg("\x1f".join, axis=1)
    ok = (obs.compartment.notna() & (obs.compartment != "drop")).to_numpy()
    codes, uniq = pd.factorize(key[ok])
    idx = np.flatnonzero(ok)
    G = sp.csr_matrix((np.ones(len(idx)), (codes, idx)), shape=(len(uniq), X.shape[0]))
    S = (G @ X).tocsc()
    ncell = np.asarray(G.sum(1)).ravel().astype(int)
    depth = np.asarray(S.sum(1)).ravel()
    # sum duplicate symbols
    gcodes, guniq = pd.factorize(pd.Series(sym).astype(str))
    M = sp.csr_matrix((np.ones(len(gcodes)), (np.arange(len(gcodes)), gcodes)), shape=(len(gcodes), len(guniq)))
    S = (S @ M).toarray()
    # median per-cell total of each group
    tot = np.asarray(X.sum(1)).ravel()
    med = pd.Series(tot[idx]).groupby(codes).median().reindex(range(len(uniq))).to_numpy()
    meta = pd.DataFrame([u.split("\x1f") for u in uniq], columns=["donor", "sample", "compartment", "stratum"])
    meta["n_cells"], meta["depth"], meta["median_umi"], meta["units"] = ncell, depth, med, units
    df = pd.concat([meta, pd.DataFrame(S, columns=guniq)], axis=1)
    return df[df.n_cells >= min_cells].reset_index(drop=True)


def finish(gse: str, pb: pd.DataFrame, donors: pd.DataFrame, cellmap: pd.DataFrame, extra_notes: str = ""):
    d = OUT / gse
    d.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(d / "pseudobulk.parquet", index=False)
    donors.to_csv(d / "donors.csv", index=False)
    cellmap.to_csv(d / "cell_type_mapping.csv", index=False)
    s = summary(pb)
    s.to_csv(d / "build_summary.csv", index=False)
    print(gse, s.to_string(index=False))
    return s


def summary(pb: pd.DataFrame) -> pd.DataFrame:
    """Per dataset counts: patients, samples, tumour compartments (>=20 cells), patients with tumour + >=2 stroma at >=50."""
    t = pb[pb.compartment.str.startswith("tumour:")]
    st = t[t.compartment.str.split(":").str[1].isin(STROMA) & (t.n_cells >= 50)]
    nst = st.groupby("donor").compartment.nunique()
    has_t = set(t.donor)
    return pd.DataFrame([dict(
        patients=pb.donor.nunique(), tumour_patients=len(has_t),
        tumour_samples=t[["donor", "sample"]].drop_duplicates().shape[0],
        pseudobulks=len(pb), cells=int(pb.n_cells.sum()),
        patients_ge1_stroma=int((nst >= 1).sum()), patients_ge2_stroma=int((nst >= 2).sum()),
        patients_tumour_epi_ge2_stroma=int(len(set(nst[nst >= 2].index) & set(
            t[t.compartment.isin(["tumour:malignant", "tumour:epithelial"]) & (t.n_cells >= 50)].donor))))])


def read_dense(path, sym_col=0, skip_cols=1):
    """Genes x cells text table (gz) -> (cells x genes csr, symbols, barcodes). First `skip_cols` columns are gene ids."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if isinstance(sym_col, int) and not isinstance(df.index, pd.RangeIndex):
        df = df.reset_index()  # header one field shorter than the rows: gene ids were read as the index
        skip_cols = 1
    sym = df.iloc[:, sym_col].astype(str).to_numpy() if isinstance(sym_col, int) else df[sym_col].astype(str).to_numpy()
    M = df.iloc[:, skip_cols:]
    return sp.csr_matrix(M.to_numpy(dtype=np.float32).T), sym, np.array(M.columns, dtype=str)


def read_h5_10x(path):
    import h5py
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = m["shape"][:]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(shape)).T.tocsr()
        ft = m["features"] if "features" in m else m
        sym = np.array([s.decode() for s in (ft["name"][:] if "name" in ft else m["gene_names"][:])])
        bc = np.array([s.decode() for s in m["barcodes"][:]])
    return X, sym, bc


def by_sample(items, G, tissue_of, stratum):
    """items: list of (sample_name, patient, X, sym, barcodes, author compartment array or None). Genes are aligned on the
    union of symbols (duplicates summed). Compartments: author's (already `drop` or a name) or the marker rule.
    Returns X (cells x union genes), symbols, obs."""
    union = pd.Index(pd.unique(np.concatenate([it[3] for it in items])))
    Xs, obs = [], []
    for name, pat, X, sym, bc, auth in items:
        comp = auth if auth is not None else np.array([MARKER_MAP.get(x, "drop") for x in marker_comp(X, sym)], dtype=object)
        P = sp.csr_matrix((np.ones(len(sym)), (np.arange(len(sym)), union.get_indexer(sym))), shape=(len(sym), len(union)))
        Xs.append((sp.csr_matrix(X) @ P).tocsr())
        tis = tissue_of(name)
        lab = np.array([f"{tis}:{c}" if c != "drop" else "drop" for c in comp])
        obs.append(pd.DataFrame({"donor": f"{G}|{pat}", "sample": name, "compartment": lab,
                                 "stratum": stratum(name) if callable(stratum) else stratum}))
    return sp.vstack(Xs).tocsr(), union.to_numpy(), pd.concat(obs, ignore_index=True)


def read_dense_stream(path, sep="\t"):
    """Big genes x cells text table (header = barcodes, first field of each row = gene) streamed into a sparse
    cells x genes csr (float32), row by row, without holding the dense table."""
    genes, rows = [], []
    with gzip.open(path, "rt") as f:
        bc = np.array(f.readline().rstrip("\n").split(sep))
        if bc[0] == "":
            bc = bc[1:]
        for line in f:
            g, rest = line.split(sep, 1)
            v = np.fromstring(rest, dtype=np.float32, sep=sep)
            nz = np.flatnonzero(v)
            rows.append((nz, v[nz]))
            genes.append(g)
    n = len(bc)
    ptr = np.cumsum([0] + [len(r[0]) for r in rows])
    G = sp.csr_matrix((np.concatenate([r[1] for r in rows]), np.concatenate([r[0] for r in rows]), ptr), shape=(len(rows), n))
    return G.T.tocsr(), np.array(genes), bc
