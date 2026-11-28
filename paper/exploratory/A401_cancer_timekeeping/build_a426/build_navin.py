"""A426 R3: Navin TNBC (CELLxGENE ceef2841-5333-46ac-92ef-ccbe0c20fe55, dataset 6f9de485; MD Anderson, 2026), raw counts.

Treatment-naive primary TNBC (collection statement); one sample per donor (the h5ad has no sample column). Classes = author
labels (Tumor -> malignant, Fibro, Endo, Peri -> mural, T -> T_NK, B -> B_plasma, Mye -> myeloid). Stratum = cohort x assay
(10x 3' v2 / v3). Cells are read donor block by donor block from `raw/X` (CSR). Output `$SCCC_DATA/navin_tnbc_cxg/`.

    build_navin.py [--limit N]
"""
from __future__ import annotations

import argparse

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import a426_common as C

SRC = C.DATA / "navin_tnbc_cxg"
H5 = SRC / "raw" / "e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad"
AUTHOR = {"Tumor": "malignant", "Fibro": "fibroblast", "Endo": "endothelial", "Peri": "mural", "T": "T_NK", "B": "B_plasma", "Mye": "myeloid"}


def cat(o, k):
    x = o[k]
    c = np.array([s.decode() if isinstance(s, bytes) else s for s in x["categories"][:]], dtype=object)
    return c[x["codes"][:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    parts, rows = [], []
    with h5py.File(H5, "r") as f:
        o = f["obs"]
        donor, ctype, assay = cat(o, "donor_id"), cat(o, "author_cell_type"), cat(o, "assay")
        sym = np.array([s.decode() for s in f["raw/var/gene_symbols"][:]])
        rx = f["raw/X"]
        indptr = rx["indptr"][:]
        order = pd.unique(donor)
        for k, d in enumerate(order[: a.limit or None]):
            i = np.flatnonzero(donor == d)
            assert (np.diff(i) == 1).all(), "donor block not contiguous"
            r0, r1 = i[0], i[-1] + 1
            s0, s1 = indptr[r0], indptr[r1]
            X = sp.csr_matrix((rx["data"][s0:s1], rx["indices"][s0:s1], indptr[r0:r1 + 1] - s0), shape=(r1 - r0, len(sym)))
            umi = np.asarray(X.sum(1)).ravel()
            comp = np.array([AUTHOR.get(x, "unassigned") for x in ctype[i]], dtype=object)
            chem = str(assay[i[0]])
            pb = C.sample_pseudobulk(X, sym, umi, comp, f"Navin|{d}", d, f"Navin|{chem}")
            parts.append(pb)
            rows.append(dict(sample=d, patient=d, cells=len(i), assay=chem, median_umi=float(np.median(umi)),
                             passes_median_umi=bool(np.median(umi) >= C.MIN_UMI)))
            if k % 10 == 0:
                C.log(f"{k + 1}/{len(order)} {d} cells {len(i)} median UMI {np.median(umi):.0f}")
    pb = C.concat_frames(parts)
    samples = pd.DataFrame(rows)
    donors = C.donor_counts(pb).merge(samples[["patient", "assay", "median_umi", "passes_median_umi"]].assign(donor=lambda d: "Navin|" + d.patient).drop(columns="patient"), on="donor", how="left")
    sums = (SRC / "raw" / "SHA256SUMS").read_text()
    prov = ("# Navin TNBC (CELLxGENE collection ceef2841-5333-46ac-92ef-ccbe0c20fe55, dataset 6f9de485-58cd-4342-bfc4-b3d3dd223aa8)\n\n"
            "- URL: https://datasets.cellxgene.cziscience.com/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad (H5AD, 4,522,982,499 bytes)\n"
            "- Downloaded 2026-11-28 (generic client agent; a first download was corrupted by two concurrent writers and was deleted and\n"
            "  redone; the obs columns of the second file were read completely without error).\n"
            "- Treatment: 'untreated samples of treatment-naive patients' (collection description; no per-donor field).\n"
            "- Parsed only (h5py), nothing downloaded was executed. Counts = `raw/X`. Built by `build_a426/build_navin.py`.\n\n## sha256\n\n```\n" + sums + "```\n")
    C.write_source(SRC, pb, samples, donors, prov)


if __name__ == "__main__":
    main()
