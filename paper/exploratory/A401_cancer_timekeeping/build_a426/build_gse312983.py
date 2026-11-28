"""A426 R4: GSE312983 (primary CRC snRNA, 10x 3', untreated resections), CellBender-filtered count matrices.

One sample per patient; classes by the `cancer_common.compartments` marker rule on the nuclei of the filtered matrix
(counts only, no clustering); sample median UMI = median total UMI of the filtered nuclei. Untreated samples only (GEO
`treatment status` = No). Output `$SCCC_DATA/gse312983/` (see `a426_common`).

    build_gse312983.py [--limit N]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import a426_common as C

SRC = C.DATA / "gse312983"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    meta = pd.read_csv(Path(__file__).resolve().parents[1] / "output/a426/build/GSE312983_samples.csv", index_col=0)
    meta = meta[meta["treatment status"] == "No"]
    files = sorted((SRC / "raw").glob("*_cellbender_filtered.h5"))
    if a.limit:
        files = files[: a.limit]
    parts, rows = [], []
    for k, f in enumerate(files):
        gsm = f.name.split("_")[0]
        pat = f.name.split("_")[1]
        X, sym, bc = C.G.read_h5_10x(f)
        X = sp.csr_matrix(X)
        umi = np.asarray(X.sum(1)).ravel()
        comp = C.G.marker_comp(X, sym, umi)
        pb = C.sample_pseudobulk(X, sym, umi, comp, f"GSE312983|{pat}", f.name[:-len("_cellbender_filtered.h5")], "GSE312983|snRNA 10x 3'")
        parts.append(pb)
        rows.append(dict(sample=f.name, gsm=gsm, patient=pat, nuclei=X.shape[0], median_umi=float(np.median(umi)),
                         stage=meta.loc[gsm, "Stage"], msi=meta.loc[gsm, "microsatellite status"], tissue=meta.loc[gsm, "tissue"],
                         passes_median_umi=bool(np.median(umi) >= C.MIN_UMI)))
        if k % 10 == 0:
            C.log(f"{k + 1}/{len(files)} {pat} nuclei {X.shape[0]} median UMI {np.median(umi):.0f}")
    pb = C.concat_frames(parts)
    samples = pd.DataFrame(rows)
    donors = C.donor_counts(pb).merge(samples[["patient", "median_umi", "passes_median_umi"]].assign(donor=lambda d: "GSE312983|" + d.patient).drop(columns="patient"), on="donor", how="left")
    sums = (SRC / "raw" / "SHA256SUMS").read_text()
    prov = ("# GSE312983: primary colorectal cancer snRNA-seq (GEO, 2026-09/10)\n\n"
            "- URL: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312983/suppl/GSE312983_RAW.tar (33 GB tar; members read by HTTP range,\n"
            "  only the CellBender-filtered h5 of the 164 untreated samples, 21.0 GB, not the full tar)\n"
            "- Downloaded 2026-11-28 (generic client agent). Tar offsets from the GEO file list (`output/a426/build/GSE312983_tar_index.csv` in the analysis repo) (header name and HDF5 signature verified per member).\n"
            "- Sample metadata: GEO series matrix (`output/a426/build/GSE312983_samples.csv`, GEO series matrix: treatment status = No).\n"
            "- Parsed only (h5py, pandas); nothing downloaded was executed. Built by `build_a426/build_gse312983.py`.\n\n## sha256 of the downloaded files\n\n```\n" + sums + "```\n")
    C.write_source(SRC, pb, samples, donors, prov)


if __name__ == "__main__":
    main()
