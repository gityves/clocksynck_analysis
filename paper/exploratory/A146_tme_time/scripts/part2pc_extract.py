"""A146 part 2 (pan-cancer) blind preparation, step 1: per-cohort cell counts and program-gene sums of tumour samples.

Blind: sums of raw counts per (patient, compartment, copy-number call) and cell counts only; no time, phase or
program-time statistic. Uses the builder's own sample selection (`cancer_common.studies.units`, majority chemistry) and
compartment rule (`marker_compartments`), tumour samples only. For each patient x marker compartment (incl. `unassigned`)
x copy-number call of epithelial cells: n_cells, genome-wide UMI (`depth`), and the counts of the program genes
(A146 part 1 table). Genes absent from a cohort are NaN.

    python part2pc_extract.py [--jobs 4] [--studies S1 S2 ...]
Output: output/part2pc_cache/<study>.csv.gz
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from joblib import Parallel, delayed

REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO), str(REPO / "paper"), str(REPO / "paper" / "common" / "datasets")]
from clocksync.genes import gene_columns  # noqa: E402
from cancer_common import studies as S  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output" / "part2pc_cache"
PROGRAMS = {
    "P1": "MKI67 CCNB1 CENPA TOP2A AURKB PCNA MCM2 BUB1B BRCA2 MRE11", "P2": "B2M PSMB8 TAPBP TAP1 HLA-A HLA-B HLA-C",
    "P3": "STAT1 IRF1 IFNGR1 IFNGR2 IFNAR2 STAT3", "P4": "TCF7 LEF1 RUNX3 BCL6", "P5": "CIITA HLA-DRA HLA-DRB1 CD74",
    "P6": "DNM2 CDC42 SYK FCGR1A FCGR2A", "P7": "CD3D CD3E CD2 TRAC", "P8": "CD68 CD14 CSF1R LYZ", "P9": "CD274",
    "P10": "PDCD1"}
GENES = sorted({g for v in PROGRAMS.values() for g in v.split()})
ALL = sys.modules[__name__]


def one(study: str) -> dict:
    from clocksync.config import limit_threads
    limit_threads(1)
    t0 = time.time()
    us, cancer = S.units(study)
    us = [u for u in S.majority_chemistry(us) if u[2] == "tumour"]
    calls = S.read_calls(study) if study in S.CNV.glob("*") or (S.CNV / f"{study}.csv.gz").exists() else None
    parts = []
    for sample, patient, kind, load, chem, _ in us:
        X, _g, umi, comp, barcodes, sym = load()
        X = sp.csr_matrix(X)
        cols = gene_columns(sym, GENES)
        if (cols < 0).any() and "MRE11" in GENES and cols[GENES.index("MRE11")] < 0:  # alias of the old symbol
            cols[GENES.index("MRE11")] = gene_columns(sym, ["MRE11A"])[0]
        call = np.full(len(comp), "", dtype=object)
        if calls is not None:
            c = calls.reindex(pd.MultiIndex.from_arrays([[sample] * len(barcodes), list(barcodes)])).to_numpy()
            e = comp == "epithelial"
            call[e] = pd.Series(c[e]).fillna("nocall").astype(str).to_numpy()
        key = pd.Series(comp.astype(str) + "|" + call.astype(str))
        codes, uniq = pd.factorize(key)
        G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
        ok = cols >= 0
        Y = np.full((len(uniq), len(GENES)), np.nan)
        Y[:, ok] = (G @ X[:, cols[ok]]).toarray()
        d = pd.DataFrame(Y, columns=GENES)
        d.insert(0, "depth", G @ np.asarray(umi, float))
        d.insert(0, "n_cells", np.bincount(codes, minlength=len(uniq)))
        k = np.array([u.split("|") for u in uniq])
        d.insert(0, "call", k[:, 1]); d.insert(0, "compartment", k[:, 0])
        d.insert(0, "sample", sample); d.insert(0, "patient", patient)
        parts.append(d)
    r = pd.concat(parts, ignore_index=True)
    agg = {c: "sum" for c in ["n_cells", "depth"]} | {g: (lambda x: x.sum(min_count=1)) for g in GENES}
    r = r.groupby(["patient", "compartment", "call"], as_index=False).agg(agg)
    r.insert(0, "cancer", cancer); r.insert(0, "study", study)
    r.to_csv(OUT / f"{study}.csv.gz", index=False)
    print(f"[{time.time() - t0:7.0f} s] {study}: {r.patient.nunique()} patients, {len(us)} tumour samples", flush=True)
    return dict(study=study, ok=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=4); ap.add_argument("--studies", nargs="*")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    studies = a.studies or [s for s in ALL.S.__dict__.get("STUDIES", [])]
    if not a.studies:
        sys.path.insert(0, str(REPO / "paper" / "common" / "datasets" / "cancer_sc"))
        import importlib.util
        spec = importlib.util.spec_from_file_location("cancer_sc_build", REPO / "paper/common/datasets/cancer_sc/build.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        studies = [s for s in m.STUDIES if not (OUT / f"{s}.csv.gz").exists()]
    print("studies:", len(studies), flush=True)
    res = Parallel(n_jobs=a.jobs)(delayed(one)(s) for s in studies)
