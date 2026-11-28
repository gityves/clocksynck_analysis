"""A426 R2: the 48 primary, untreated, unsorted CRC donors of the CELLxGENE single-cell CRC atlas (core atlas, dataset
05a8c945, 30.9 GB h5ad) that are not held on disk. The dataset is not in the Census, so its h5ad is read by HTTP range:
the obs columns (0.3 GB), then only the `raw/X` rows (CSR) of the selected cells, summed on the fly into donor x sample x
class pseudobulks (no cell-level matrix kept; about 4 GB transferred, never the full file). COUNTS ONLY.

Donors (counts only, as the scout `a426_step1`): studies not held on disk (Joanito 2022, Li 2023, Qi 2022, MUI Innsbruck,
Guo 2022, Chen 2021, Yang 2023); tumour samples with enrichment `naive` and treatment `naive`; >= 2 stromal classes
(Fibroblast, Endothelial cell, Pericyte of `cell_type_middle_crc_atlas`) >= 50 cells, >= 100 malignant cells, >= 50 immune
cells (the registered 48). Classes of the pseudobulks: the A401 CELLxGENE rule on the ontology `cell_type` (build_census_tumour
`compartment`: fibroblast / endothelial (not lymphatic) / mural / T_NK / B_plasma / myeloid / malignant; mast and other: dropped).
Stratum = cohort x platform. Output `$SCCC_DATA/crc_atlas_cxg/`.

    build_crc_atlas.py [--limit-blocks N]
"""
from __future__ import annotations

import argparse
import re
import time

import fsspec
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import a426_common as C

SRC = C.DATA / "crc_atlas_cxg"
URL = "https://datasets.cellxgene.cziscience.com/4a8b9568-965e-46b8-a427-baab6bf018e5.h5ad"
UA = {"User-Agent": "python-requests/2.31"}
NOT_HELD = ("Joanito_2022_Nat_Genet", "Li_2023_Cancer_Cell", "Qi_2022_Nat_Commun", "MUI_Innsbruck", "Guo_2022_JCI_Insight",
            "Chen_2021_Cell", "Yang_2023_Front_Oncol")
OBS = ["donor_id", "study_id", "sample_id", "sample_type", "enrichment_cell_types", "treatment_status_before_resection",
       "cell_type_middle_crc_atlas", "cell_type", "n_counts", "platform", "matrix_type"]
IMM = {"CD8", "CD4", "Plasma cell", "B cell", "Macrophage", "Monocyte", "Treg", "Neutrophil", "NK", "gamma-delta", "Dendritic cell",
       "Mast cell", "ILC", "NKT", "Eosinophil"}
MYELOID = re.compile(r"macrophage|monocyte|dendritic|myeloid|microglia|kupffer", re.I)
TNK = re.compile(r"\bT cell|T-regulatory|regulatory T|T helper|gamma-delta|T-cell|T lymphocyte|natural killer|\bNK\b|\bNKT\b|\bMAIT|mucosal invariant", re.I)
BPLASMA = re.compile(r"\bB cell|B-cell|plasma cell|plasmablast|plasmacyte", re.I)


def ontology_class(ct: str):
    """build_census_tumour.compartment for tumour cells of a dataset with a malignant label."""
    if re.search(r"pericyte|smooth muscle|mural", ct, re.I):
        return "mural"
    if re.search(r"lymphatic", ct, re.I):
        return None
    if re.search(r"endothelial", ct, re.I):
        return "endothelial"
    if re.search(r"fibroblast|stromal cell|stromal", ct, re.I) and not re.search(r"malignant|neoplastic|tumou?r cell|cancer cell", ct, re.I):
        return "fibroblast"
    if MYELOID.search(ct):
        return "myeloid"
    if BPLASMA.search(ct):
        return "B_plasma"
    if TNK.search(ct):
        return "T_NK"
    return "malignant" if re.search(r"malignant|neoplastic|tumou?r cell|cancer cell", ct, re.I) else None


def read_obs(h):
    out = {}
    for k in OBS:
        x = h["obs"][k]
        if isinstance(x, h5py.Group):
            cat = np.array([c.decode() if isinstance(c, bytes) else c for c in x["categories"][:]], dtype=object)
            codes = x["codes"][:]
            out[k] = np.where(codes < 0, "NA", cat[np.clip(codes, 0, None)]).astype(str)
        else:
            out[k] = x[:]
        C.log("obs", k)
    return pd.DataFrame(out)


def select_donors(d: pd.DataFrame):
    """Counts-only donor selection of the registration (the 48)."""
    w = d[(d.sample_type == "tumor") & (d.enrichment_cell_types == "naive") & (d.treatment_status_before_resection == "naive")
          & d.study_id.isin(NOT_HELD)].copy()
    w["key"] = w.study_id + "|" + w.donor_id
    mid = w.cell_type_middle_crc_atlas
    w["cls"] = np.where(mid.isin(IMM), "immune", np.where(mid == "Fibroblast", "fibroblast", np.where(mid == "Endothelial cell", "endothelial",
               np.where(mid == "Pericyte", "mural", np.where(w.cell_type == "malignant cell", "tumour", "other")))))
    t = pd.crosstab(w.key, w.cls)
    for c in ("immune", "fibroblast", "endothelial", "mural", "tumour"):
        if c not in t:
            t[c] = 0
    t["n_stroma_classes50"] = (t[["fibroblast", "endothelial", "mural"]] >= 50).sum(1)
    t["selected"] = (t.n_stroma_classes50 >= 2) & (t.tumour >= 100) & (t.immune >= 50)
    return t, w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-blocks", type=int, default=0)
    a = ap.parse_args()
    SRC.mkdir(parents=True, exist_ok=True)
    fs = fsspec.filesystem("http", client_kwargs={"headers": UA})
    f = h5py.File(fs.open(URL, block_size=16 << 20, cache_type="readahead"), "r")
    d = read_obs(f)
    t, w = select_donors(d)
    t.to_csv(SRC / "donor_selection_counts.csv")
    sel = set(t.index[t.selected])
    C.log(f"{len(sel)} donors selected of {len(t)} tumour donors of the not-held studies")
    assert len(sel) == 48, len(sel)
    d["key"] = d.study_id + "|" + d.donor_id
    m = (d.key.isin(sel) & (d.sample_type == "tumor") & (d.enrichment_cell_types == "naive") & (d.treatment_status_before_resection == "naive")).to_numpy()
    rows = np.flatnonzero(m)
    assert (d.matrix_type.to_numpy()[rows] == "raw counts").all()
    cls = np.array([ontology_class(c) or "drop" for c in d.cell_type.to_numpy()[rows]], dtype=object)
    smp = d.sample_id.to_numpy()[rows]
    key = d.key.to_numpy()[rows]
    plat = d.platform.to_numpy()[rows]
    gkey = pd.Series(key + "\x1f" + smp + "\x1f" + cls + "\x1f" + plat)
    ok = cls != "drop"
    codes, uniq = pd.factorize(gkey[ok])
    gidx = np.full(len(rows), -1)
    gidx[np.flatnonzero(ok)] = codes
    umi_obs = d.n_counts.to_numpy(float)[rows]
    pos = {r: i for i, r in enumerate(rows)}
    var = f["raw/var/feature_name"]
    cat = np.array([c.decode() if isinstance(c, bytes) else c for c in var["categories"][:]], dtype=object)
    sym = cat[var["codes"][:]].astype(str)
    ng = len(sym)
    acc = np.zeros((len(uniq), ng))
    rx = f["raw/X"]
    ip = rx["indptr"]
    gap = 200
    br = np.flatnonzero(np.diff(rows) > gap)
    starts, ends = np.r_[rows[0], rows[br + 1]], np.r_[rows[br], rows[-1]] + 1
    if a.limit_blocks:
        starts, ends = starts[: a.limit_blocks], ends[: a.limit_blocks]
    umi_chk = []
    t0 = time.time()
    for k, (s, e) in enumerate(zip(starts, ends)):
        for p0 in range(s, e, 20_000):
            p1 = min(e, p0 + 20_000)
            indptr = ip[p0:p1 + 1]
            s0, s1 = indptr[0], indptr[-1]
            X = sp.csr_matrix((rx["data"][s0:s1], rx["indices"][s0:s1], indptr - s0), shape=(p1 - p0, ng))
            need = np.flatnonzero(m[p0:p1])
            if not len(need):
                continue
            ii = np.array([pos[p0 + j] for j in need])
            g = gidx[ii]
            use = g >= 0
            Gm = sp.csr_matrix((np.ones(use.sum()), (g[use], need[use])), shape=(len(uniq), p1 - p0))
            acc += (Gm @ X).toarray()
            umi_chk.append((np.asarray(X[need].sum(1)).ravel(), umi_obs[ii]))
        C.log(f"block {k + 1}/{len(starts)} rows {e - s}, {time.time() - t0:.0f}s")
    chk = np.concatenate([a_ for a_, _ in umi_chk]), np.concatenate([b_ for _, b_ in umi_chk])
    C.log(f"UMI from counts vs obs n_counts: max |diff| {np.abs(chk[0] - chk[1]).max():.1f}, median ratio {np.median(chk[0] / np.maximum(chk[1], 1)):.3f}")
    # pseudobulk frame (symbols summed), sample median UMI over all selected tumour cells of the sample
    meta = pd.DataFrame([u.split("\x1f") for u in uniq], columns=["donor", "sample", "cls", "platform"])
    meta["compartment"] = "tumour:" + meta.cls
    meta["donor"] = "CRCatlas|" + meta.donor
    meta["stratum"] = "CRCatlas|" + meta.platform
    cnt = pd.Series(gkey[ok]).value_counts()
    meta["n_cells"] = [cnt[u] for u in uniq]
    meta["depth"] = acc.sum(1)
    sm = pd.Series(umi_obs).groupby(smp).median()
    meta["sample_median_umi"] = meta["sample"].map(sm).to_numpy()
    meta["median_umi"] = np.nan
    meta["units"] = "counts"
    gcodes, guniq = pd.factorize(pd.Series(sym))
    M = sp.csr_matrix((np.ones(len(gcodes)), (np.arange(len(gcodes)), gcodes)), shape=(len(gcodes), len(guniq)))
    S = np.asarray((sp.csr_matrix(acc) @ M).todense(), dtype=np.float32)
    pb = pd.concat([meta[["donor", "sample", "compartment", "stratum", "n_cells", "depth", "median_umi", "units", "sample_median_umi"]],
                    pd.DataFrame(S, columns=guniq)], axis=1)
    pb = pb[pb.n_cells >= 20].reset_index(drop=True)
    samples = (pd.DataFrame({"sample": smp, "donor": key, "platform": plat}).drop_duplicates("sample")
               .assign(median_umi=lambda x: x["sample"].map(sm).to_numpy(), passes_median_umi=lambda x: x.median_umi >= C.MIN_UMI))
    donors = C.donor_counts(pb)
    prov = ("# CELLxGENE single-cell CRC atlas, core atlas (collection 3a844375-..., dataset 05a8c945-bc12-414f-960d-a31943bbcdd1)\n\n"
            "- URL: https://datasets.cellxgene.cziscience.com/4a8b9568-965e-46b8-a427-baab6bf018e5.h5ad (H5AD, 30,875,155,333 bytes; NOT downloaded in full)\n"
            "- Read on 2026-11-28 by HTTP range requests (generic client agent): obs columns (counts-only donor selection) and the `raw/X` CSR rows of the\n"
            "  tumour cells of the 48 selected donors; summed into pseudobulks on the fly. The dataset is not in the Census (2025-11-08 / 2025-11-17).\n"
            "- No local copy of the h5ad; therefore no file sha256. Remote object: the dataset version id 4a8b9568-965e-46b8-a427-baab6bf018e5 in the URL.\n"
            "- Parsed only (h5py), nothing downloaded was executed. Built by `build_a426/build_crc_atlas.py`.\n")
    C.write_source(SRC, pb, samples, donors, prov)


if __name__ == "__main__":
    main()
