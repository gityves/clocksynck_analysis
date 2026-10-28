"""Holm over the two responsive-gene tests of the A109 / A138 note (2026-10-28): H109.1-R (ventilated donors) and
H138.4-R. Reads the outputs of `compute_h109_responsive.py` and `score_h138_responsive.py`; secondary family, it does not
alter the registered verdicts.

    python holm_responsive.py
Output: output/responsive/holm_responsive.csv (this folder) .
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
A109 = HERE.parent / "A109_coupling_nulls" / "output" / "responsive" / "h109_1_responsive.csv"
A138 = HERE / "output" / "responsive" / "h138_4_responsive.csv"


def main():
    a = pd.read_csv(A109)
    a = a[(a.group == "vent") & a.test.str.startswith("H109.1-R")].iloc[0]
    b = pd.read_csv(A138)
    b = b[b.test.str.startswith("H138.4-R")].iloc[0]
    t = pd.DataFrame([dict(hypothesis="H109.1-R (ventilated)", statistic=a.C_stat, responsive_median=a.null_median,
                           responsive_q95=a.null_q95, p=a.p_one_sided, n_sets=a.n_sets),
                      dict(hypothesis="H138.4-R", statistic=b.statistic, responsive_median=b.null_median,
                           responsive_q95=b.null_q95, p=b.p_one_sided, n_sets=b.n_null)])
    t = t.sort_values("p").reset_index(drop=True)
    m, run = len(t), 0.0
    adj = []
    for r, p in enumerate(t.p):
        run = max(run, min(1.0, p * (m - r)))
        adj.append(run)
    t["rank"] = range(1, m + 1)
    t["holm_adj"] = adj
    t["above_responsive_q95"] = t.statistic > t.responsive_q95
    t["reject_0_05"] = t.holm_adj < 0.05
    t["verdict"] = t.reject_0_05.map({True: "above responsive sets: not a GC / immediate-early response",
                                      False: "clock-gene-shaped shared state, not separable from a shared acute response"})
    t.to_csv(HERE / "output" / "responsive" / "holm_responsive.csv", index=False)
    print(t.round(5).to_string(index=False))


if __name__ == "__main__":
    main()
