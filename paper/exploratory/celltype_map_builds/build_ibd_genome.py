"""All-gene donor x tissue x condition x class pseudobulks of IBDverse (cell-type program A3: genome-wide twin of the
`ibd_ileum` and `ibd_rectum` clock-gene tables). Self-contained (run with `python -I`; the input is downloaded data).

Grouping, classes and conditions exactly as the panel build (legacy desync/16_ibd/scripts/build_ibd.py, A14):
`author_cell_type_category` (T/ILC -> T, B, Plasma, Myeloid, Enterocyte, Colonocyte, Secretory, Stem); Mesenchymal split
by `cell_type` into endothelial (labels containing 'endothelial') and stromal; condition 'healthy' (disease normal),
'crohn_inflamed', 'crohn_uninflamed' (`author_inflammation_status`). Streams `raw/X` (integer UMI, CSR) in row chunks.
Every gene of raw/var instead of the panel; X stored sparse float32 (integer sums, exact).

Output (OUT_DIR): pseudobulk.h5ad (obs as the panel build; var feature_name), donors.csv, PROVENANCE.md.

    python -I build_ibd_genome.py H5AD OUT_DIR
"""

import datetime
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

CAT = {"T/ILC": "T", "B": "B", "Plasma": "plasma", "Myeloid": "myeloid", "Enterocyte": "enterocyte",
       "Colonocyte": "colonocyte", "Secretory": "secretory", "Stem": "stem"}
CHUNK = 100_000


def cat(h, col):
    g = h["obs"][col]
    cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][()]], dtype=object)
    codes = g["codes"][()]
    return np.where(codes < 0, "NA", cats[np.clip(codes, 0, None)]).astype(str)


def main(h5, out):
    out = Path(out)
    with h5py.File(h5, "r") as h:
        sym = h["raw/var/feature_name"]
        sym = np.array([x.decode() if isinstance(x, bytes) else x for x in
                        (sym[()] if isinstance(sym, h5py.Dataset) else
                         np.array(sym["categories"][()])[sym["codes"][()]])], dtype=object)
        c1, ct = cat(h, "author_cell_type_category"), cat(h, "cell_type")
        cls = pd.Series(c1).map(CAT).to_numpy(dtype=object)
        mes = c1 == "Mesenchymal"
        endo = np.char.find(ct.astype(str), "endothelial") >= 0
        cls[mes & endo] = "endothelial"
        cls[mes & ~endo] = "stromal"
        dis, infl = cat(h, "disease"), cat(h, "author_inflammation_status")
        cond = np.where(dis == "normal", "healthy", np.where(infl == "inflamed", "crohn_inflamed", "crohn_uninflamed"))
        key = pd.Series(cat(h, "donor_id")) + "|" + pd.Series(cat(h, "tissue")) + "|" + pd.Series(cond) + "|" + \
            pd.Series(cls).astype(str)
        ok = pd.notna(cls)
        codes, uniq = pd.factorize(key.where(ok, None))
        n = len(codes)
        X = h["raw/X"]
        indptr = X["indptr"][()]
        G = len(uniq)
        sums, depth, ncell = np.zeros((G, len(sym))), np.zeros(G), np.zeros(G, int)
        for a in range(0, n, CHUNK):
            b = min(a + CHUNK, n)
            lo, hi = indptr[a], indptr[b]
            M = sp.csr_matrix((X["data"][lo:hi], X["indices"][lo:hi], indptr[a:b + 1] - lo), shape=(b - a, len(sym)))
            c = codes[a:b]
            keep = c >= 0
            A = sp.csr_matrix((np.ones(keep.sum()), (c[keep], np.flatnonzero(keep))), shape=(G, b - a))
            R = (A @ M).tocoo()  # in place: no dense groups x genes temporary per chunk
            np.add.at(sums, (R.row, R.col), R.data)
            depth += A @ np.asarray(M.sum(1)).ravel()
            ncell += np.bincount(c[keep], minlength=G)
            print(a, n, flush=True)
    parts = pd.Series(uniq).str.split("|", expand=True)
    obs = pd.DataFrame(dict(donor_id=parts[0], tissue=parts[1] + "|" + parts[2], tissue_general=parts[1],
                            cell_type=parts[3], assay="10x 3' v3", n_cells=ncell, depth=depth))
    S = sp.csr_matrix(sums)
    del sums
    A = ad.AnnData(X=S.astype(np.float32), obs=obs,
                   var=pd.DataFrame(dict(feature_name=sym.astype(str)), index=pd.Index(np.arange(len(sym)).astype(str))))
    out.mkdir(parents=True, exist_ok=True)
    A.write_h5ad(out / "pseudobulk.h5ad")
    dd = pd.DataFrame(dict(donor_id=parts[0], cond=parts[2])).drop_duplicates()
    d = dd.groupby("donor_id").cond.agg(lambda s: "normal" if set(s) == {"healthy"} else "Crohn disease")
    d.rename("disease").rename_axis("donor_id").reset_index().to_csv(out / "donors.csv", index=False)
    (out / "PROVENANCE.md").write_text(
        f"# All-gene pseudobulk of IBDverse (CELLxGENE 0e5c75da)\n- Built {datetime.date.today()} by "
        f"`paper/exploratory/celltype_map_builds/build_ibd_genome.py` from {h5} (raw/X), grouping as the panel build "
        f"($SCCC_DATA/geo/IBDverse/pseudobulk, legacy build_ibd.py).\n- {n} cells, {len(obs)} donor x tissue x condition x "
        f"class rows, {d.size} donors, {len(sym)} genes (raw/var feature_name).\n")
    print(A)


if __name__ == "__main__":
    main(*sys.argv[1:3])
