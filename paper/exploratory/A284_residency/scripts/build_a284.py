"""A284 tables: all-gene donor x compartment pseudobulks with resident and circulating T classes (labels registered in
paper/docs/registered_experiments.md, A284).

- breast (Kumar) and Kong ileum: the genome tables' own builder (`_atlas/genome.py`) with the registered Census
  `cell_type` labels placed first in the class list (first match wins), so their cells form `T_resident` /
  `T_circulating`; every other class is built as before (non-immune rows checked identical to the genome table).
- IBDverse ileum and rectum (healthy): the T author labels of the local raw file, summed per donor x tissue x class
  (`ibd_tsub`, below), added to the healthy rows of `ibd_ileum_genome` / `ibd_rectum_genome`.

    python -I build_a284.py ibd RAW_H5AD OUT_DIR     # raw IBDverse -> group pseudobulk of the two T classes
    (the score script calls `table(name)` for every atlas)
"""

from __future__ import annotations

import datetime
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common/datasets/_atlas"))

T_LABELS = {
    "organ_breast": ({"CD8-positive, alpha-beta memory T cell"},
                     {"effector memory CD8-positive, alpha-beta T cell", "effector memory CD4-positive, alpha-beta T cell",
                      "CD4-positive, alpha-beta T cell"}),
    "kong_ileum": ({"intraepithelial lymphocyte"}, {"naive T cell"}),
    "ibd": ({"CD8 TRM TRGC2", "Tissue resident GD T"}, {"CD4 CM/naive", "CD4 CM", "CD8 naive"}),
}
GENOME = {"organ_breast": "organ_breast_genome", "kong_ileum": "kong_ileum_genome",
          "ibd_ileum": "ibd_ileum_genome", "ibd_rectum": "ibd_rectum_genome"}
IBD_T = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/atlas_wg/ibdverse_tsub")
CHUNK = 100_000


def _exact(labels):
    return "^(" + "|".join(sorted(map(_re_escape, labels))) + ")$"


def _re_escape(s):
    import re
    return re.escape(s)


def _patched(name):
    """genome.build of `name` with the registered T labels first in the class list."""
    import genome as G
    import organ as O
    res, circ = T_LABELS[name]
    base = O.COARSE
    O.COARSE = (("T_resident", _exact(res)), ("T_circulating", _exact(circ))) + tuple(base)
    try:
        return G.build(GENOME[name])
    finally:
        O.COARSE = base


def table(name):
    """(pb, check) of an A284 atlas: all genes, the source genome table's rows plus T_resident / T_circulating."""
    import clocksync as cs
    clock = list(cs.gtex().genes)
    src, _ = importlib.import_module(f"common.datasets.{GENOME[name]}.build").build()
    if name in ("organ_breast", "kong_ileum"):
        pb, _ = _patched(name)
    else:
        tissue = name.removeprefix("ibd_")
        t = ibd_rows(tissue, set(src.donor), clock + [c for c in src.columns if c not in clock])
        pb = pd.concat([src, t], ignore_index=True)
    k = ["donor", "compartment"]
    new = pb.compartment.str.endswith(("T_resident", "T_circulating"))
    a = src.set_index(k)
    b = pb[~new].set_index(k)
    common = a.index.intersection(b.index)
    cols = clock + ["n_cells", "depth"]
    same = (a.loc[common, cols].astype(float).to_numpy() == b.loc[common, cols].astype(float).to_numpy()).all(1)
    check = dict(atlas=name, rows_source=len(a), rows_new_table_old_classes=len(b), common=len(common),
                 identical_rows=int(same.sum()),
                 resident_donors=int(pb[pb.compartment.str.endswith("T_resident")].donor.nunique()),
                 circulating_donors=int(pb[pb.compartment.str.endswith("T_circulating")].donor.nunique()))
    return pb, check


def ibd_rows(tissue, donors, genes):
    """Healthy donor x tissue x T-class rows (>= 20 cells) from the `ibd_tsub` group pseudobulk, in the genome table's
    format (donor `ibdverse:<id>`, compartment `<tissue>:<class>`, stratum = assay)."""
    import anndata as ad
    a = ad.read_h5ad(IBD_T / "pseudobulk.h5ad")
    o = a.obs
    m = (o.tissue == f"{tissue}|healthy").to_numpy() & (o.n_cells >= 20).to_numpy()
    X = pd.DataFrame(a.X[m].toarray(), columns=a.var.feature_name.astype(str))
    X = X.loc[:, ~X.columns.duplicated()]
    rows = pd.DataFrame(dict(donor="ibdverse:" + o.donor_id[m].astype(str).to_numpy(),
                             compartment=f"{tissue}:" + o.cell_type[m].astype(str).to_numpy(),
                             stratum="10x 3' v3", n_cells=o.n_cells[m].to_numpy(), depth=o.depth[m].to_numpy()))
    rows = pd.concat([rows, X.reindex(columns=[g for g in genes if g in X.columns]).fillna(0.0)], axis=1)
    return rows[rows.donor.isin(donors)].reset_index(drop=True)


def build_ibd(h5, out, classes=None, title="IBDverse T subclasses (A284)"):
    """Raw IBDverse -> group pseudobulk (donor x tissue|condition x class) of registered author labels, all genes.
    `classes`: {class: set of author labels} (default: the A284 T classes). Run with `python -I` (downloaded data)."""
    import anndata as ad
    import h5py
    import scipy.sparse as sp
    out = Path(out)
    if classes is None:
        res, circ = T_LABELS["ibd"]
        classes = {"T_resident": res, "T_circulating": circ}
    with h5py.File(h5, "r") as h:
        def cat(col):
            g = h["obs"][col]
            cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][()]], dtype=object)
            codes = g["codes"][()]
            return np.where(codes < 0, "NA", cats[np.clip(codes, 0, None)]).astype(str)
        sym = h["raw/var/feature_name"]
        sym = np.array([x.decode() if isinstance(x, bytes) else x for x in
                        (sym[()] if isinstance(sym, h5py.Dataset) else
                         np.array(sym["categories"][()])[sym["codes"][()]])], dtype=object)
        lab = cat("author_cell_type")
        cls = np.full(len(lab), "", dtype=object)
        for k, v in classes.items():
            cls[np.isin(lab, list(v))] = k
        dis, infl = cat("disease"), cat("author_inflammation_status")
        cond = np.where(dis == "normal", "healthy", np.where(infl == "inflamed", "crohn_inflamed", "crohn_uninflamed"))
        key = pd.Series(cat("donor_id")) + "|" + pd.Series(cat("tissue")) + "|" + pd.Series(cond) + "|" + pd.Series(cls)
        codes, uniq = pd.factorize(key.where((cls != "") & (cond == "healthy"), None))
        n, G = len(codes), len(uniq)
        X = h["raw/X"]
        indptr = X["indptr"][()]
        sums, depth, ncell = np.zeros((G, len(sym))), np.zeros(G), np.zeros(G, int)
        for a in range(0, n, CHUNK):
            b = min(a + CHUNK, n)
            c = codes[a:b]
            keep = c >= 0
            if not keep.any():
                continue
            lo, hi = indptr[a], indptr[b]
            M = sp.csr_matrix((X["data"][lo:hi], X["indices"][lo:hi], indptr[a:b + 1] - lo), shape=(b - a, len(sym)))
            A = sp.csr_matrix((np.ones(keep.sum()), (c[keep], np.flatnonzero(keep))), shape=(G, b - a))
            R = (A @ M).tocoo()
            np.add.at(sums, (R.row, R.col), R.data)
            depth += A @ np.asarray(M.sum(1)).ravel()
            ncell += np.bincount(c[keep], minlength=G)
    parts = pd.Series(uniq).str.split("|", expand=True)
    obs = pd.DataFrame(dict(donor_id=parts[0], tissue=parts[1] + "|" + parts[2], cell_type=parts[3], n_cells=ncell,
                            depth=depth))
    A = ad.AnnData(X=sp.csr_matrix(sums).astype(np.float32), obs=obs,
                   var=pd.DataFrame(dict(feature_name=sym.astype(str)), index=pd.Index(np.arange(len(sym)).astype(str))))
    out.mkdir(parents=True, exist_ok=True)
    A.write_h5ad(out / "pseudobulk.h5ad")
    (out / "PROVENANCE.md").write_text(
        f"# {title}\n- Built {datetime.date.today()} by "
        f"`paper/exploratory/A284_residency/scripts/build_a284.py ibd` from {h5} (raw/X): healthy cells, "
        f"author_cell_type {({k: sorted(v) for k, v in classes.items()})}; donor x tissue|condition x "
        f"class sums, all {len(sym)} genes.\n")
    print(obs.groupby(["tissue", "cell_type"]).size(), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "ibd":
        build_ibd(sys.argv[2], sys.argv[3])
