"""Post hoc development check (2026-11-19, gap 3): does restricting intervals to hosts with high estimated reliability
(df_R >= r0) restore coverage of the df-corrected c? Reads A163 validation raw rows (no new simulation).

    python dev_reliability_rule.py
Output: ../output/validate/dev_reliability_rule.csv (cells x r0: coverage, n, median bias).
"""
from pathlib import Path

import numpy as np
import pandas as pd

V = Path(__file__).resolve().parents[1] / "output" / "validate"
d = pd.concat([pd.read_csv(f) for f in sorted((V / "raw").glob("*.csv"))])
d = d[(d.role == "target") & ~d.cell.str.startswith("eq")].copy()
d["hit"] = (d.df_lo <= d.true_c) & (d.true_c <= d.df_hi)
rows = []
for r0 in (0.0, 0.8, 0.85, 0.9):
    for (des, cell), x in d[d.df_R >= r0].groupby(["design", "cell"]):
        rows.append(dict(r0=r0, design=des, cell=cell, coverage=x.hit.mean(), n=len(x), median_bias=np.median(x.df_est - x.true_c)))
t = pd.DataFrame(rows)
t.to_csv(V / "dev_reliability_rule.csv", index=False)
print(t.groupby("r0").apply(lambda g: pd.Series(dict(cells=len(g), cov_ge_090=int((g.coverage >= 0.9).sum()),
                                                   min_cov=g.coverage.min(), max_abs_bias=g.median_bias.abs().max())),
                            include_groups=False).round(3).to_string())
