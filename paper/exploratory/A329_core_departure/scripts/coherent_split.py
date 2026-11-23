"""A329 erratum (post hoc diagnosis, 2026-11-23): split each tissue's personal departure into a coherent clock-state
part, invariant to the projector (L X = I), and a gene-level part, amplified by the projector's noise gain.

Model per tissue: D_full = C + N, D_core = C + g N, with g the tissue's noise gain from the matched sets
(median own departure of the matched sets at the core positions over that at all 14 positions). Solved for C and N.
The unweighted noise gain of the projectors, trace(L L') core / full, is reported beside it.

    python coherent_split.py
Input: ../output/tissues.csv (A329). Output: ../output/coherent_split.csv (per tissue), coherent_summary.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync._core.plane import projector

OUT = Path(__file__).resolve().parents[1] / "output"


def main():
    t = pd.read_csv(OUT / "tissues.csv", index_col="tissue")
    g = t.own_matched_median_core / t.own_matched_median_full
    C = (g * t.own_full - t.own_core) / (g - 1)
    s = pd.DataFrame(dict(noise_gain=g, coherent=C, noise_clock=t.own_full - C,
                          noise_matched=t.own_matched_median_full, coherent_over_matched=C / t.own_matched_median_full))
    s.to_csv(OUT / "coherent_split.csv")
    w = cs.gtex()
    Lf, Lc = projector(w), projector(w.drop(["PER1", "PER2", "NR1D1"]))
    summ = s.median().to_frame("median").T.assign(
        tissues=len(s), coherent_positive=int((s.coherent > 0).sum()), coherent_min=s.coherent.min(),
        coherent_max=s.coherent.max(), projector_gain_unweighted=np.trace(Lc @ Lc.T) / np.trace(Lf @ Lf.T),
        own_full_median=t.own_full.median(), own_core_median=t.own_core.median(),
        matched_full_median=t.own_matched_median_full.median(), matched_core_median=t.own_matched_median_core.median())
    summ.to_csv(OUT / "coherent_summary.csv", index=False)
    print(summ.round(3).T.to_string())


if __name__ == "__main__":
    main()
