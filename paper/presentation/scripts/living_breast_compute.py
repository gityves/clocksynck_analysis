"""Slide data (results, R1.3 A): breast surgery cohort (GSE233242, 42 patients, tumour + adjacent normal, resection 08-18 h).

Re-runs the A200 readout (registered, scored 2026-11-03; `exploratory/A200_breast_surgery_rerun/scripts/score_a200.py`, same
functions and seeds) to keep what the registered script does not write: per patient, the recorded resection hour and the
internal time with the 13:00 rule (estimates' circular mean set to 13:00, no recorded time used), and the Spearman of each of
the 200 expression-matched pseudo-clock sets (the null of H200.1). Recorded hours of a living cohort are not controlled access.

    paper/common/run.sh paper/presentation/scripts/living_breast_compute.py
Outputs: paper/presentation/output/living/{breast_patients,breast_null,breast_summary}.csv
"""
import importlib.util
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("a200", PAPER / "exploratory/A200_breast_surgery_rerun/scripts/score_a200.py")
A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = PAPER / "presentation" / "output" / "living"


def main():
    A.limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    x, key, donors = A.tables()
    rho, e13, est, truth = A.readout(A.pseudobulk(x, key, A.CLOCK), donors)
    off = A.IT.cdiff(A.PRIOR, (np.angle(np.exp(1j * A.IT.W * est).mean()) / A.IT.W) % 24)
    pd.DataFrame(dict(resection_h=truth, internal_h=(est + off) % 24,
                      error_h=A.IT.cdiff(est + off, truth))).to_csv(OUT / "breast_patients.csv", index=False)
    null = []
    for s in A.pseudo_sets(x, key, 200, np.random.default_rng([A.SEED, 1])):
        try:
            null.append(A.readout(A.pseudobulk(x, key, s), donors)[0])
        except Exception:
            null.append(np.nan)
    pd.DataFrame(dict(spearman=null)).to_csv(OUT / "breast_null.csv", index=False)
    pd.DataFrame([dict(n=len(truth), spearman=rho, median_abs_err_h=e13, null_median=np.nanmedian(null),
                       null_q95=np.nanquantile(null, .95), null_failed=int(np.isnan(null).sum()))]).to_csv(
        OUT / "breast_summary.csv", index=False)
    print(f"rho {rho:.3f}, err13 {e13:.2f} h, null median {np.nanmedian(null):.3f}", flush=True)


if __name__ == "__main__":
    main()
