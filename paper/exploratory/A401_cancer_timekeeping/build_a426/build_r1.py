"""A426 R1: genome-wide donor x compartment pseudobulks (all tumour classes) of the on-disk patients that A403 did not score.

Candidate patients = the A425 inventory (counts and flags only; `output/a426/build/inventory_patients_not_in_A403.csv`) in the
R1 categories of the registration: no cancer-cell pseudobulk >= 50 cells in an A401 source; cancers outside the A403 set; sources
never entered in A401; the 8 CELLxGENE breast primaries labelled "exocrine gland"; the Wu / Pal breast patients (`build_breast`).
Excluded as in A401: copies and dedup removals, treated / metastatic / adjacent-normal samples, CD45-sorted (or otherwise sorted)
samples, normalised-only data. COUNTS ONLY: no clock-gene statistic or immune score.

  units    cohorts of `cancer_common.studies.units` (A388 cohorts and scBaseCount studies never entered): samples of the candidate
           patients (A149 sample rules inside `units`; sorted samples removed here), classes by the cohort's labels (author or
           marker rule, `cancer_sc.build._labels` without copy-number calls), counts summed per patient x class x chemistry
  geo / census / bassez   the A401 builds (`$SCCC_DATA/{geo_tumour_a401,cellxgene_tumour_a401,bassez2021_3ca}`), candidate rows
Output `$SCCC_DATA/a426_r1/<source>/{pseudobulk.parquet,samples.csv}` and `output/a426/build/r1_candidates.csv`.

    build_r1.py --list | --source SRC [--source SRC ...] | --all [--jobs 2]
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import a426_common as C

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
EXP = PAPER / "exploratory"
BUILD = HERE.parent / "output" / "a426" / "build"
OUT = C.DATA / "a426_r1"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(EXP / "A388_malignant_tumour_host/scripts"))
for p in ("A227_tme_compartment_pairing/scripts", "A149_paired_tumour_normal/scripts", "A151_fibroblast_pairing/scripts",
          "A226_stroma_pairing_confirm/scripts", "A336_two_readouts/scripts", "A170_sc_breast_er/scripts"):
    sys.path.insert(0, str(EXP / p))

CATEGORIES = ["no_malignant_ge50", "cancer_not_P401", "cancer_unmapped", "never_entered", "A403_table_no_immune_readout"]
EXOCRINE = ["breastcancer4patie_fbdd8c17", "stromalcelldiversi_68b6114f"]
NOT_REGISTERED = {"GBM", "GBM/glioma", "NB", "hepatoblastoma", "lymphoma", "metastatic_label"}  # glioma (no fibroblasts), paediatric / haematological
UNSORTED = {"unknown", "unsorted", "live CD235a-negative (paper); not CD45-sorted"}
SORTED_SOURCES = {"GSE155109": "FACS fraction (EC / stroma enriched)", "GSE183219": "stroma-sorted (CD45-depleted)",
                  "SRP327758": "CD45+ / CD45- sorted", "SRP301444": "CD45+ / CD45- sorted", "SRP329496": "CD45+ / CD45- sorted",
                  "SRP335254": "sorted (other)"}


def candidates() -> pd.DataFrame:
    A = pd.read_csv(BUILD / "inventory_patients_not_in_A403.csv")
    A["patient"] = A.patient.astype(str)
    c = A[A.category.isin(CATEGORIES) & (A.n_stroma_classes >= 1) & (A.n_immune_cells >= 1)]
    x = A[(A.category == "metastatic_site") & A.source.isin(EXOCRINE)].assign(category="exocrine_gland_breast_primary")
    c = pd.concat([c, x], ignore_index=True)
    c["registered_cancer"] = ~c.cancer_label.isin(NOT_REGISTERED)
    c["sorted_source"] = c.source.map(SORTED_SOURCES).fillna("")
    keep = ["source", "patient", "kind", "category", "cancer_label", "n_stroma_classes", "n_immune_cells", "chem", "registered_cancer", "sorted_source"]
    return c[keep]


# ---- cohorts of cancer_common.studies.units -------------------------------------------------------------------------

def build_units(study, patients):
    import a149_lib as L
    import compute_a226 as K226
    K226.install()
    warnings.filterwarnings("ignore")
    inv = pd.read_csv(BUILD / "inventory_samples.csv")
    sm = inv[inv.study == study].set_index("sample")
    if study in ("RCC_5af90777", "Qian2020_CRC", "Qian2020_LUAD") or study.startswith("GSE"):
        pass
    us, _ = L.ST.units(study)
    us = [u for u in L.ST.majority_chemistry(us) if u[1] in set(patients)]
    parts, rows = [], []
    for sample, patient, kind, load, chem, ncell in us:
        if kind != "tumour":
            continue
        sort = sm.sorting.get(sample, "unknown") if len(sm) else "unknown"
        if study in SORTED_SOURCES or (isinstance(sort, str) and sort not in UNSORTED):
            rows.append(dict(sample=sample, patient=patient, chem=chem, n_cells=ncell, used=False, why=f"sorted: {sort}"))
            continue
        X, _g, umi, comp, barcodes, sym = load()
        lab = L.B._labels(comp, kind, sample, barcodes, None)
        cls = np.array([x.split(":", 1)[1] for x in lab], dtype=object)
        pb = C.sample_pseudobulk(X, sym, umi, cls, f"{study}|{patient}", sample, f"{study}|{chem}")
        parts.append(pb)
        rows.append(dict(sample=sample, patient=patient, chem=chem, n_cells=len(umi), used=True, why="", sample_median_umi=float(np.median(umi))))
    return (C.concat_frames(parts) if parts else None), pd.DataFrame(rows)


def build_table(source, patients, kind):
    """geo / census / bassez: candidate rows of the A401 builds."""
    don = [f"{source}|{p}" for p in patients]
    if kind == "geo":
        pb = pd.read_parquet(C.DATA / "geo_tumour_a401" / source / "pseudobulk.parquet", filters=[("donor", "in", don)])
        pb = pb[pb.compartment.str.startswith("tumour:")]
        pb["sample_median_umi"] = np.nan
    elif kind == "bassez":
        pb = pd.read_parquet(C.DATA / "bassez2021_3ca" / "pseudobulk.parquet", filters=[("donor", "in", [f"bassez2021|{p}" for p in patients])])
        pb = pb[(pb.timepoint == "Pre") & pb.compartment.str.startswith("tumour:")].drop(columns=["timepoint", "expansion", "cohort", "BC_type"], errors="ignore")
        pb["units"] = "counts"
        pb["sample_median_umi"] = np.nan
    else:
        pb = pd.read_parquet(C.DATA / "cellxgene_tumour_a401" / f"{source}.parquet", filters=[("donor", "in", don)])
        pb = pb[pb.compartment.str.startswith("tumour:")]
        ok_organ = source in EXOCRINE or True  # the inventory's candidate patients are the primary-site ones
        pb["sample"] = pb.donor
        pb["sample_median_umi"] = pb.median_umi
        pb["median_umi"] = np.nan
        pb["units"] = "counts"
        pb = pb.drop(columns=[c for c in ("n_cells_sample", "disease", "tissue", "tissue_general", "suspension", "tcga", "obs_raw_sum") if c in pb])
    meta = ["donor", "sample", "compartment", "stratum", "n_cells", "depth", "median_umi", "units", "sample_median_umi"]
    return pb[meta + [c for c in pb.columns if c not in meta]].reset_index(drop=True)


def do_source(args):
    source, kind, patients = args
    d = OUT / source
    if (d / "pseudobulk.parquet").exists():
        return source, "exists"
    try:
        if kind in ("a388", "scbc_never_entered"):
            pb, samples = build_units(source, patients)
        else:
            pb, samples = build_table(source, patients, kind), pd.DataFrame()
        if pb is None or not len(pb):
            return source, "no rows"
        d.mkdir(parents=True, exist_ok=True)
        pb.to_parquet(d / "pseudobulk.parquet", index=False)
        samples.to_csv(d / "samples.csv", index=False)
        return source, f"ok {pb.donor.nunique()} donors {len(pb)} rows"
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        return source, f"FAILED {type(e).__name__}: {e}"[:300]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--source", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    c = candidates()
    c.to_csv(BUILD / "r1_candidates.csv", index=False)
    todo = c[c.registered_cancer & (c.sorted_source == "")]
    if a.list:
        print(todo.groupby(["kind", "source"]).size().to_string())
        return
    srcs = a.source or (sorted(todo.source.unique()) if a.all else [])
    jobs = [(s, todo[todo.source == s].kind.iloc[0], list(todo[todo.source == s].patient)) for s in srcs]
    from joblib import Parallel, delayed
    for r in Parallel(n_jobs=a.jobs)(delayed(do_source)(j) for j in jobs):
        C.log(*r)


if __name__ == "__main__":
    main()
