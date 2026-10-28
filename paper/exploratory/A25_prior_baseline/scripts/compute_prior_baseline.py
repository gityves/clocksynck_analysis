"""Roadmap 4.4 (descriptive): prior-only baseline for V25.5, breast per-patient internal time (method median |error| 1.1 h).

Same patients as the V25.5 primary run (paper/results/figure_1/output/windows/v25_5_patients.csv: est_day = tangent hour with the
12 h flip resolved by the hospital prior, true_h = resection hour). Baselines as A123c (compute_truth_iii.py):
  always 13:00 (the hospital prior's window centre), best constant (median recorded hour; uses truth), window centre of the
  patients' own recorded hours (circular mean; uses truth). Method: median |error| after (a) the offset minimising it (0.05 h grid)
  and (b) the offset setting the circular mean of the estimates to 13:00 (no recorded time). Aggregates only.
Output: output/prior_baseline.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
P = PAPER / "results/figure_1/output/windows/v25_5_patients.csv"
OUT = Path(__file__).resolve().parents[1] / "output"
W = 2 * np.pi / 24
cdiff = lambda a, b: (np.asarray(a) - np.asarray(b) + 12) % 24 - 12
cmean = lambda h: float((np.angle(np.exp(1j * W * np.asarray(h)).mean()) / W) % 24)

d = pd.read_csv(P)
p, tr = d.est_day.to_numpy(float), d.true_h.to_numpy(float) % 24
grid = np.arange(0, 24, 0.05)
o_min = grid[np.argmin([np.median(np.abs(cdiff(p + g, tr))) for g in grid])]
o_13 = float(cdiff(13.0, cmean(p)))
med = lambda e: float(np.median(np.abs(e)))
rows = [
    dict(estimator="method, offset minimising the median error (as V25.5)", offset_or_constant_h=o_min, median_abs_error_h=med(cdiff(p + o_min, tr))),
    dict(estimator="method, circular mean of estimates set to 13:00", offset_or_constant_h=o_13, median_abs_error_h=med(cdiff(p + o_13, tr))),
    dict(estimator="always 13:00 (prior window centre)", offset_or_constant_h=13.0, median_abs_error_h=med(cdiff(13.0, tr))),
    dict(estimator="always the circular mean of recorded hours (window centre, uses truth)", offset_or_constant_h=cmean(tr), median_abs_error_h=med(cdiff(cmean(tr), tr))),
    dict(estimator="best constant (minimises the median error, grid 0.05 h)", offset_or_constant_h=0.0, median_abs_error_h=0.0),
]
best = grid[np.argmin([np.median(np.abs(cdiff(g, tr))) for g in grid])]
rows[-1].update(offset_or_constant_h=float(best), median_abs_error_h=med(cdiff(best, tr)))
r = pd.DataFrame(rows).assign(n=len(d), recorded_hour_sd_h=float(np.sqrt(np.mean(cdiff(tr, cmean(tr)) ** 2))))
r.to_csv(OUT / "prior_baseline.csv", index=False)
print(r.round(3).to_string())
