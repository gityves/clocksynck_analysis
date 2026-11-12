"""A247 preparation: genome-wide patient x class pseudobulks of the Extended GBmap (droplet libraries only).

Registered in paper/docs/registered_experiments.md (A247). Computes no clock statistic: it sums raw counts and
recodes the patient metadata.

Source: `$SCCC_DATA/cellxgene_gbmap/extended_gbmap_230c8701.h5ad`, read backed (`raw/X` in row chunks).
Raw counts: `raw/X` holds integer UMI counts (checked on every chunk read: all values equal their rounding, minimum 1);
`X` holds normalised values (non-integer) and is not used. `raw/var.feature_name` holds the symbols (the index is Ensembl).

Cells: droplet assays (10x 3' v2 / v3, 10x 5' v1, Drop-seq); microwell-seq, BD Rhapsody, CEL-seq2, STRT-seq and
Smart-seq2 are out. Classes (published annotation): `malignant` = annotation_level_1 Neoplastic; `myeloid` = TAM-BDM,
TAM-MG, Mono, DC; `oligodendrocyte` = annotation_level_3 Oligodendrocyte (OPC out); `T_NK` = CD4/CD8, NK (descriptive).
One pseudobulk per patient (`donor_id`) x class with >= 100 cells. Stratum = study (`author`) | suspension (cells /
nuclei).

MGMT label, recoded from the metadata alone (the per-cell `MGMT` column, all cells of the donor): M, Meth ->
methylated; UM, Un, NM -> unmethylated; unknown, NT -> no call. A donor with both calls: `conflict`; with no call:
`unknown`.

Choices not fixed by the registration:
- `is_primary_data` False cells (Ruiz2022, Darmanis2017; 42,888 cells of 15 donors) are not duplicated inside the
  extended file (every cell index and `observation_joinid` is unique): they are kept, so each cell counts once.
  Dropping them would remove those donors.
- Patient = `donor_id`. A donor with cells in both suspension types (one donor) takes the suspension of the majority
  of its kept cells for all its classes; its counts are summed over suspensions.
- MGMT label and stage use all cells of the donor in the file (any assay), a patient property; `NT` is read as "not
  tested" (no call). Stage: the single non-unknown label of the donor, `mixed` if several, `unknown` if none;
  `recurrent` groups `Recurrent` and `Recurrent neo-aPD1`.
- Genes: symbols resolved with `clocksync.genes.canonical`, duplicated symbols summed.
- Depth = total UMI of the cells of the pseudobulk (all genes).

Outputs (`../output/`): pseudobulks.parquet (rows: donor x class; meta columns donor, compartment, stratum, study,
suspension, n_cells, depth, then genes), patients.csv, prep_summary.csv.

    python prep_a247.py
"""
from __future__ import annotations

import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

from clocksync.config import limit_threads
from clocksync.genes import canonical

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
SRC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_gbmap/extended_gbmap_230c8701.h5ad")
MIN_CELLS, CHUNK = 100, 50_000
CLASS_OF = {"TAM-BDM": "myeloid", "TAM-MG": "myeloid", "Mono": "myeloid", "DC": "myeloid",
            "Oligodendrocyte": "oligodendrocyte", "CD4/CD8": "T_NK", "NK": "T_NK"}
METHYLATED, UNMETHYLATED = {"M", "Meth"}, {"UM", "Un", "NM"}
META = ["donor", "compartment", "stratum", "study", "suspension", "n_cells", "depth"]


def read_obs() -> pd.DataFrame:
    cols = ["assay", "suspension_type", "author", "donor_id", "annotation_level_1", "annotation_level_3", "stage",
            "MGMT", "is_primary_data"]
    with h5py.File(SRC, "r") as f:
        o = pd.DataFrame({c: read_elem(f["obs"][c]) for c in cols}, index=read_elem(f["obs"]["_index"]))
    return o.astype({c: str for c in cols})


def mgmt_call(labels: pd.Series) -> str:
    s = {x.strip() for x in labels}
    m, u = bool(s & METHYLATED), bool(s & UNMETHYLATED)
    return "conflict" if m and u else "methylated" if m else "unmethylated" if u else "unknown"


def stage_of(labels: pd.Series) -> str:
    s = sorted(set(labels) - {"unknown"})
    return "unknown" if not s else s[0] if len(s) == 1 else "mixed"


def main() -> None:
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    o = read_obs()
    n_all = len(o)
    pat = o.groupby("donor_id").agg(MGMT_labels=("MGMT", lambda s: ",".join(sorted(set(s)))), mgmt_call=("MGMT", mgmt_call),
                                    stage=("stage", stage_of))
    pat["stage_group"] = pat.stage.replace({"Recurrent neo-aPD1": "recurrent", "Recurrent": "recurrent", "Primary": "primary"})
    droplet = o.assay.str.startswith("10x") | (o.assay == "Drop-seq")
    cls = pd.Series(np.where(o.annotation_level_1 == "Neoplastic", "malignant", o.annotation_level_3.map(CLASS_OF)), index=o.index)
    keep = (droplet & cls.notna()).to_numpy()
    # suspension of a patient: majority of its kept cells
    k = o[keep].assign(cls=cls[keep])
    susp = k.groupby("donor_id").suspension_type.agg(lambda s: s.value_counts().index[0])
    study = k.groupby("donor_id").author.agg(lambda s: s.value_counts().index[0])
    assert (k.groupby("donor_id").author.nunique() == 1).all(), "a donor occurs in two studies"
    n = k.groupby(["donor_id", "cls"]).size()
    groups = n[n >= MIN_CELLS].rename("n_cells").reset_index()
    gid = {(d, c): i for i, (d, c) in enumerate(zip(groups.donor_id, groups.cls))}
    cell_group = np.full(n_all, -1, dtype=np.int64)
    sel = np.flatnonzero(keep)
    cell_group[sel] = [gid.get((d, c), -1) for d, c in zip(o.donor_id.to_numpy()[sel], cls.to_numpy()[sel])]
    print(f"{n_all} cells, {int(droplet.sum())} droplet, {int((cell_group >= 0).sum())} in {len(groups)} pseudobulks", flush=True)

    with h5py.File(SRC, "r") as f:
        names = read_elem(f["raw/var"])["feature_name"].astype(str).map(canonical).to_numpy()
        genes, inv = np.unique(names, return_inverse=True)
        collapse = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(genes)))
        X = f["raw/X"]
        indptr = X["indptr"][:]
        assert len(indptr) - 1 == n_all
        S = np.zeros((len(groups), len(inv)))
        depth = np.zeros(len(groups))
        t0 = time.time()
        for a in range(0, n_all, CHUNK):
            b = min(a + CHUNK, n_all)
            g = cell_group[a:b]
            if (g >= 0).any():
                s, e = indptr[a], indptr[b]
                d = X["data"][s:e].astype(np.float64)
                assert d.min() >= 1 and np.all(d == np.round(d)), "raw/X is not integer counts"
                m = sp.csr_matrix((d, X["indices"][s:e], indptr[a:b + 1] - s), shape=(b - a, len(inv)))
                r = np.flatnonzero(g >= 0)
                ind = sp.csr_matrix((np.ones(len(r)), (g[r], r)), shape=(len(groups), b - a))
                S += (ind @ m).toarray()
                depth += np.asarray(ind @ np.asarray(m.sum(1)).ravel()).ravel()
            print(f"{b}/{n_all} cells, {time.time() - t0:.0f}s", flush=True)
    G = pd.DataFrame((collapse.T @ S.T).T, columns=genes).round().astype(np.int64)
    meta = pd.DataFrame({"donor": groups.donor_id.to_numpy(), "compartment": groups.cls.to_numpy(),
                         "stratum": (study[groups.donor_id].to_numpy() + "|" + susp[groups.donor_id].to_numpy()),
                         "study": study[groups.donor_id].to_numpy(), "suspension": susp[groups.donor_id].to_numpy(),
                         "n_cells": groups.n_cells.to_numpy(), "depth": depth})
    pb = pd.concat([meta, G], axis=1)
    assert np.allclose(pb[genes].sum(1), pb.depth), "gene sums differ from total UMI"
    pb.to_parquet(OUT / "pseudobulks.parquet", index=False)
    p = pat.loc[sorted(set(pb.donor))].join(pd.concat([study.rename("study"), susp.rename("suspension")], axis=1))
    p["stratum"] = p.study + "|" + p.suspension
    for c in ("malignant", "myeloid", "oligodendrocyte", "T_NK"):
        p[f"n_{c}"] = pb[pb.compartment == c].set_index("donor").n_cells.reindex(p.index).fillna(0).astype(int)
    p.index.name = "donor"
    p.reset_index().to_csv(OUT / "patients.csv", index=False)
    summ = pb.groupby("compartment").agg(pseudobulks=("donor", "size"), strata=("stratum", "nunique"),
                                         median_cells=("n_cells", "median"), median_depth=("depth", "median")).reset_index()
    summ.to_csv(OUT / "prep_summary.csv", index=False)
    print(summ.to_string(index=False), p.mgmt_call.value_counts().to_string(), p.stage_group.value_counts().to_string(), sep="\n")


if __name__ == "__main__":
    main()
