"""The single-cell cancer cohorts: eligible samples, per-sample loaders, compartments, copy-number calls.

Two sources, both processed files that are already on disk (no download, no mapping):
- scBaseCount (Arc Institute, 2026 release): one h5ad per sample (SRX) in `scbasecount_full/<study>/` (the lab's earlier
  downloads, `$SCCC_SCBC_FULL`, then `$SCCC_DATA/scbc_full/<study>/`); sample annotation (study, patient, sample type,
  treatment, median UMI) from the curated inventory `$SCCC_DATA/inventory/inventory_samples.csv`; 10x chemistry from the
  Arc study table (`ARC`);
- GEO cohorts: one h5ad per cohort, `$SCCC_DATA/geo/<GSE>/<GSE>.h5ad` (raw UMI counts; obs sample, patient, sample_type,
  sorting, treatment, total_umi, [chemistry]).
(The scripts that produced these files, legacy/scripts/inventory/*, fetch_scbc.py and geo/build_geo.py, are data
acquisition and are not ported.)

Ported from legacy/scripts/run_study.py (`select`, `scbc_loader`, `units`, `majority_chemistry`) and
legacy/scripts/cnv_calls.py (`calls`), with no change of behaviour.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path[:0] = [str(Path(__file__).resolve().parents[1]), str(Path(__file__).resolve().parents[3])]
from clocksync.genes import read_counts  # noqa: E402
from common.paths import DATA  # noqa: E402

from cancer_common import cnv  # noqa: E402
from cancer_common.compartments import marker_compartments  # noqa: E402

SCBC_FULL = (Path(os.environ.get("SCCC_SCBC_FULL", DATA.parent / "scbasecount_full")), DATA / "scbc_full")
GEO = DATA / "geo"
CNV = DATA / "cnv"
INVENTORY = DATA / "inventory" / "inventory_samples.csv"
ARC = Path(os.environ.get("SCCC_ARC_STUDIES",
                          "/home/maxine/Documents/paychere/scRNAseq_database/output/datasets/scbasecount_human_studies.csv"))
GEO_CANCER = {"GSE178341": "colorectal", "GSE166555": "colorectal", "GSE131907": "lung", "GSE181919": "head and neck",
              "GSE183904": "gastric", "GSE149614": "liver", "GSE200997": "colorectal", "GSE184362": "thyroid",
              "GSE132465": "colorectal", "GSE144735": "colorectal", "GSE159115": "kidney", "GSE181294": "prostate",
              "GSE326201": "liver", "GSE198315": "head and neck", "GSE212966": "pancreatic", "GSE270680": "gastric",
              "GSE294300": "colorectal", "GSE236696": "colorectal", "GSE315534": "colorectal", "GSE281736": "thyroid",
              "GSE271690": "colorectal", "GSE222315": "bladder", "GSE335070": "ampullary", "GSE203067": "esophageal"}
MIN_UMI = 1500  # median UMI per sample
KEEP = ("tumour", "adjacent_normal", "healthy", "metastasis")


def select(study: str) -> pd.DataFrame:
    """Inventory rows of the scBaseCount study's eligible samples (biopsy or resection, not treated, patient known,
    sample type in `KEEP`, median UMI >= `MIN_UMI`, gene-expression library), with the path of the file (`file`)."""
    inv = pd.read_csv(INVENTORY, low_memory=False)
    s = inv[(inv.study == study) & inv.source.isin(["biopsy", "resection"]) & (inv.treatment != "treated")
            & inv.patient.notna() & inv.sample_type.isin(KEEP) & (inv.median_umi >= MIN_UMI)]
    if "library" in s:
        s = s[s.library.fillna("gex") == "gex"]
    s = s.assign(file=[next((d / study / f"{x}.h5ad" for d in SCBC_FULL if (d / study / f"{x}.h5ad").exists()), None)
                       for x in s.srx])
    return s[s.file.notna()]


def scbc_loader(path):
    def load():
        a = read_counts(path)
        sym = a.var["gene_symbols"].astype(str).to_numpy()
        umi = a.obs.umi_count_Unique.to_numpy(float)
        return a.X, a.var_names, umi, marker_compartments(a.X, sym, umi), a.obs.index, sym
    return load


def units(study: str, where: str | None = None) -> tuple[list, str]:
    """(sample, patient, sample_type, load, chemistry, n_cells) per eligible sample, and the cancer type. `load()`
    returns (X, genes, total_umi, compartment, barcodes, symbols). `where` = 'obs_column=value' keeps a subset
    (GEO only). GEO cohorts: not treated, not CD45+-sorted, patient known, sample median UMI >= `MIN_UMI`."""
    if study.startswith("GSE"):
        a = read_counts(GEO / study / f"{study}.h5ad")
        if where:
            col, val = where.split("=", 1)
            a = a[(a.obs[col].astype(str) == val).to_numpy()].copy()
        o = a.obs
        ok = (o.sample_type.isin(KEEP) & (o.treatment != "treated") & ~o.sorting.str.contains(r"CD45\+|CD45p", regex=True)
              & o.patient.astype(str).ne(""))
        a = a[ok.to_numpy()].copy()
        depth = a.obs.groupby("sample", observed=True).total_umi.median()
        a = a[a.obs["sample"].isin(depth[depth >= MIN_UMI].index).to_numpy()].copy()
        sym = a.var_names.astype(str).to_numpy()
        comp = marker_compartments(a.X, sym, a.obs.total_umi.to_numpy(float))

        def sub(i):
            return lambda: (a.X[i], a.var_names, a.obs.total_umi.to_numpy(float)[i], comp[i], a.obs.index[i], sym)
        out = []
        chem = a.obs["chemistry"].astype(str) if "chemistry" in a.obs else pd.Series("", index=a.obs.index)
        for smp, g in a.obs.groupby("sample", observed=True):
            i = np.flatnonzero(a.obs["sample"].to_numpy() == smp)
            out.append((str(smp), str(g.patient.iloc[0]), str(g.sample_type.iloc[0]), sub(i), chem.iloc[i[0]], len(i)))
        return out, GEO_CANCER.get(study, "")
    s = select(study)
    tech = pd.read_csv(ARC, usecols=["srx_accession", "tech_10x"]).set_index("srx_accession").tech_10x.astype(str)
    cancer = s.cancer_type.mode().iloc[0] if s.cancer_type.notna().any() else ""
    return [(r.srx, str(r.patient), r.sample_type, scbc_loader(r.file), tech.get(r.srx, ""), int(r.n_cells))
            for r in s.itertuples()], cancer


def majority_chemistry(us: list) -> list:
    """Drop the samples of a patient that are not of its majority chemistry (by cells)."""
    cells = {}
    for _, pat, _, _, chem, n in us:
        cells.setdefault(pat, {}).setdefault(chem, 0)
        cells[pat][chem] += n
    keep = {pat: max(c, key=c.get) for pat, c in cells.items()}
    return [u for u in us if u[4] == keep[u[1]]]


def calls(study: str, positions: pd.DataFrame | None = None) -> pd.DataFrame | None:
    """Copy-number calls of the epithelial cells of a study (`cnv.study_calls`: malignant / normal / ambiguous, against the
    normal epithelium of the other patients), recomputed from the counts. None for studies without enough normal
    epithelium. Calls of the cohorts as run for the paper: `read_calls`."""
    us, _ = units(study)
    return cnv.study_calls([(smp, pat, kind, (lambda ld=ld: ld()[:5])) for smp, pat, kind, ld, *_ in us],
                           cnv.positions() if positions is None else positions)


def read_calls(study: str) -> pd.Series:
    """Stored calls `$SCCC_DATA/cnv/<study>.csv.gz`, indexed by (sample, barcode)."""
    return pd.read_csv(CNV / f"{study}.csv.gz").set_index(["sample", "barcode"]).cnv_call
