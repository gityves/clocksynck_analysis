"""A140 Q1 revision: genome-wide pseudobulks of the tumour T_NK and tumour myeloid compartments of the 21 gate-passing A17
cancer cohorts (raw UMI sums over all genes), the input of the factor calibration (`calibrate_factor.py`).

Registered in paper/docs/registered_experiments.md (A140, "Q1 revision"). Same cells, labels and filters as the main
`cancer_sc` tables (`cnv=False`): `cancer_common.studies.units` (eligible samples), `studies.majority_chemistry`,
`cancer_sc.build._labels` with calls=None (marker-gene compartments), sum over the samples of a patient, >= 50 cells. A
(patient, compartment) is kept only if it is in the main table (which also applies the cohort-level rule on
compartments); the differences between the two are reported in `build_check.csv`. Processed unit by unit (one sample, or
one GEO cohort at a time) to keep memory low.

Run: RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A140_dip_cyclops/build_tumour_immune_genome.py --jobs 4
Output (git-ignored): output/q1/tumour_immune_genome/<study>.npz (sparse patients x genes), <study>_rows.csv (patient,
compartment, n_cells, depth = genome-wide UMI), <study>_genes.txt, build_check.csv.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path[:0] = [str(PAPER), str(PAPER / "common" / "datasets"), str(PAPER / "common" / "datasets" / "cancer_sc")]
OUT = HERE / "output" / "q1" / "tumour_immune_genome"
COMPS = ("tumour:T_NK", "tumour:myeloid")
MIN_CELLS = 50


def targets() -> list[str]:
    d = pd.read_csv(PAPER / "results" / "figure_5" / "output" / "cancer" / "a17_cancer.csv")
    g = d.drop_duplicates("target")
    return sorted(g[g.gate_passed].target.str.removeprefix("cancer_"))


def build_study(study: str) -> dict:
    from common import data
    from cancer_common import studies as S
    from cancer_sc.build import _labels
    from clocksync.config import limit_threads
    limit_threads(1)
    t0 = time.time()
    main, _ = data.load("cancer_sc", study=study, cnv=False)
    main = main[main.compartment.isin(COMPS)]
    us, _cancer = S.units(study, None)
    us = S.majority_chemistry(us)
    acc: dict[tuple, dict] = {}  # (patient, compartment) -> {symbol: count}, plus n_cells, depth
    for sample, patient, kind, load, _chem, _n in us:
        X, _g, umi, comp, barcodes, sym = load()
        lab = _labels(comp, kind, sample, barcodes, None)
        sel = np.flatnonzero(np.isin(lab, COMPS))
        if len(sel) == 0:
            continue
        X = sp.csr_matrix(X[sel] if sp.issparse(X) else np.asarray(X)[sel]).astype(np.float64)
        labs = lab[sel]
        # one column per unique symbol (duplicated symbols are summed)
        usym, inv = np.unique(sym.astype(str), return_inverse=True)
        M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(usym)))
        for c in COMPS:
            m = labs == c
            if not m.any():
                continue
            v = np.asarray(X[m].sum(0)).ravel() @ M
            a = acc.setdefault((patient, c), dict(v=sp.csr_matrix((1, 0)), g=np.array([], dtype=str), n=0, depth=0.0))
            tot = pd.Series(np.concatenate([a["v"].toarray().ravel() if a["v"].shape[1] else np.zeros(0), v]),
                            index=np.concatenate([a["g"], usym])).groupby(level=0).sum()
            a.update(v=sp.csr_matrix(tot.to_numpy()[None]), g=tot.index.to_numpy(str), n=a["n"] + int(m.sum()),
                     depth=a["depth"] + float(umi[sel][m].sum()))
        del X
    keep = {k: a for k, a in acc.items() if a["n"] >= MIN_CELLS}
    genes = sorted(set().union(*[set(a["g"]) for a in keep.values()])) if keep else []
    gi = {g: i for i, g in enumerate(genes)}
    rows, mats = [], []
    for (patient, c), a in sorted(keep.items()):
        r = np.zeros(len(genes))
        r[[gi[g] for g in a["g"]]] = a["v"].toarray().ravel()
        mats.append(sp.csr_matrix(r))
        rows.append(dict(donor=patient, compartment=c, n_cells=a["n"], depth=a["depth"]))
    rt = pd.DataFrame(rows)
    ref = {(d, c) for d, c in zip(main.donor.astype(str), main.compartment)}
    got = {(d, c) for d, c in zip(rt.donor.astype(str), rt.compartment)} if len(rt) else set()
    # restrict to the main table's (patient, compartment) pairs
    use = np.array([(d, c) in ref for d, c in zip(rt.donor.astype(str), rt.compartment)], bool) if len(rt) else np.zeros(0, bool)
    OUT.mkdir(parents=True, exist_ok=True)
    sp.save_npz(OUT / f"{study}.npz", sp.vstack(mats).tocsr()[np.flatnonzero(use)] if use.any() else sp.csr_matrix((0, len(genes))))
    rt[use].to_csv(OUT / f"{study}_rows.csv", index=False)
    (OUT / f"{study}_genes.txt").write_text("\n".join(genes))
    chk = dict(study=study, main_rows=len(ref), built_rows=len(got), kept=int(use.sum()), only_main=len(ref - got),
               only_built=len(got - ref), n_genes=len(genes), seconds=round(time.time() - t0))
    print(chk, flush=True)
    return chk


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--studies", default="")
    a = ap.parse_args()
    st = a.studies.split(",") if a.studies else targets()
    print(f"{len(st)} cohorts: {' '.join(st)}", flush=True)
    todo = [s for s in st if not (OUT / f"{s}_rows.csv").exists()]
    res = Parallel(n_jobs=a.jobs, backend="loky", return_as="generator_unordered")(delayed(build_study)(s) for s in todo)
    chk = pd.DataFrame(list(res))
    old = pd.read_csv(OUT / "build_check.csv") if (OUT / "build_check.csv").exists() else pd.DataFrame()
    pd.concat([old, chk]).drop_duplicates("study", keep="last").sort_values("study").to_csv(OUT / "build_check.csv", index=False)
    print("finished", flush=True)


if __name__ == "__main__":
    main()
