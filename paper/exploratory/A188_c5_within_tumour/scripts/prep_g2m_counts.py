"""A188 prep, step 1: genome-wide counts of the Hallmark G2M genes (clock genes removed, A97 set) in the malignant
pseudobulk of every patient of every H138.1 cohort, rebuilt from the per-cell counts with the rules of
`cancer_sc.build` (majority chemistry, cnv calls, >= 50 cells). The malignant compartment of a cohort is that of
H138.1 (`q3_lib.study_table`: called malignant cells if >= 5 patients have them, else all tumour epithelium).

    python -I prep_g2m_counts.py [study ...]
Outputs (../output/g2m/): g2m_<study>.csv (patient, n_cells, depth, one column per G2M gene) and check_<study>.csv
(rebuilt vs stored n_cells / depth of the malignant pseudobulk).
"""
import sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A138_cancer_estimability" / "scripts"))
import numpy as np, pandas as pd, scipy.sparse as sp  # noqa: E402
import q3_lib as q  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402

OUT = HERE.parent / "output" / "g2m"
S = q.sh.sc_builder()


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def one(study):
    stored, mal = q.study_table(study)
    cnv = mal == q.MAL
    src = "tumour:malignant" if cnv else "tumour:epithelial"
    stored = stored[stored.compartment == mal].set_index("donor")
    g2m = q.sh.hallmark_g2m()
    from common.datasets.cancer_common import studies as CS
    us, _ = CS.units(study)
    us = CS.majority_chemistry(us)
    calls = CS.read_calls(study) if cnv else None
    rows = {}
    for sample, patient, kind, load, chem, _ in us:
        if kind != "tumour":
            continue
        X, _g, umi, comp, barcodes, sym = load()
        lab = S._labels(comp, kind, sample, barcodes, calls)
        m = lab == src
        if not m.any():
            continue
        cols = gene_columns(sym, g2m)
        X = sp.csr_matrix(X)[np.flatnonzero(m)]
        r = rows.setdefault(patient, dict(n_cells=0, depth=0.0, cnt=np.zeros(len(g2m)), present=cols >= 0))
        r["n_cells"] += int(m.sum())
        r["depth"] += float(np.asarray(umi, float)[m].sum())
        r["cnt"] += np.where(cols >= 0, np.asarray(X[:, np.where(cols >= 0, cols, 0)].sum(0)).ravel(), 0.0)
    t = pd.DataFrame([dict(patient=p, n_cells=r["n_cells"], depth=r["depth"], **dict(zip(g2m, r["cnt"]))) for p, r in rows.items()])
    t = t[t.n_cells >= S.MIN_CELLS].set_index("patient")
    chk = pd.DataFrame({"n_cells_rebuilt": t.n_cells, "depth_rebuilt": t.depth}).join(
        stored[["n_cells", "depth"]].rename(columns={"n_cells": "n_cells_stored", "depth": "depth_stored"}), how="outer")
    chk.insert(0, "study", study)
    chk.to_csv(OUT / f"check_{study}.csv")
    t.reset_index().to_csv(OUT / f"g2m_{study}.csv", index=False)
    ok = chk.dropna()
    say(study, f"patients stored {len(stored)}, rebuilt {len(t)}, matched {len(ok)}; "
               f"n_cells equal {np.mean(ok.n_cells_rebuilt == ok.n_cells_stored):.3f}, depth rel.diff max {np.max(np.abs(ok.depth_rebuilt / ok.depth_stored - 1)):.2e}")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    gate = q.a17_gate()
    todo = sys.argv[1:] or [d["study"] for d in (q.q1_design(s) for s in sorted(s for s, g in gate.items() if g)) if d]
    say("cohorts", len(todo), todo)
    for s in todo:
        if (OUT / f"g2m_{s}.csv").exists():
            continue
        one(s)
    say("done")
