"""A203 scoring: V203.1 (cells saturate) and V203.2 (donors limit) from the split-half errors of `compute_a203.py`.

    python score_a203.py [--subsample]
Writes (../output): errors.csv (arm x level: median |offset difference| (h) and |c difference| over readout compartments
and splits, with the between-half circular correlation of offsets and Pearson of c, medians over splits), m95.csv,
families.csv (cells arm per compartment), summary.csv (V203.1-2).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "output"
HOURS = 24 / (2 * np.pi)
SAT, DON = 1.10, 1.25


def circ_r(a, b):
    """Circular correlation (Fisher-Lee / Jammalamadaka) of two angle vectors in hours."""
    a, b = np.asarray(a) / HOURS, np.asarray(b) / HOURS
    sa, sb = np.sin(a - np.angle(np.mean(np.exp(1j * a)))), np.sin(b - np.angle(np.mean(np.exp(1j * b))))
    d = np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2))
    return float(np.sum(sa * sb) / d) if d > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--subsample", action="store_true")
    args = ap.parse_args()
    raw = pd.read_csv(OUT / ("raw_subsample.csv" if args.subsample else "raw.csv"))
    tag = "_subsample" if args.subsample else ""
    per_split = []
    for (arm, lv, s), g in raw.groupby(["arm", "level", "split"]):
        g = g.dropna(subset=["off_1", "off_2", "c_1", "c_2"])
        per_split.append(dict(arm=arm, level=lv, split=s, r_off=circ_r(g.off_1, g.off_2),
                              r_c=float(np.corrcoef(g.c_1, g.c_2)[0, 1]) if len(g) > 2 else np.nan))
    ps = pd.DataFrame(per_split)
    E = raw.groupby(["arm", "level"]).agg(err_offset_h=("d_off", "median"), err_c=("d_c", "median"),
                                          n=("d_off", "size")).reset_index()
    E = E.merge(ps.groupby(["arm", "level"])[["r_off", "r_c"]].median().reset_index(), on=["arm", "level"])
    E.to_csv(OUT / f"errors{tag}.csv", index=False)

    cells = E[E.arm == "cells"].set_index("level")
    donors = E[E.arm == "donors"].set_index("level")
    rows = []
    for metric in ("err_offset_h", "err_c"):
        r1 = cells.loc[100, metric] / cells.loc[-1, metric]
        rows.append(dict(criterion="V203.1", metric=metric, ratio=r1, threshold=f"<= {SAT}", passes=bool(r1 <= SAT),
                         note="err(m = 100) / err(all cells)"))
        n_full = donors.index.max()
        r2 = donors.loc[100, metric] / donors.loc[n_full, metric]
        rows.append(dict(criterion="V203.2", metric=metric, ratio=r2, threshold=f">= {DON}", passes=bool(r2 >= DON),
                         note=f"err(N = 100) / err(N = {n_full})"))
    S = pd.DataFrame(rows)
    S.to_csv(OUT / f"summary{tag}.csv", index=False)

    m95 = []
    lv = [x for x in cells.index if x > 0]
    for metric in ("err_offset_h", "err_c"):
        full = cells.loc[-1, metric]
        x = np.array(sorted(lv), float)
        y = cells.loc[sorted(lv), metric].to_numpy(float)
        hit = np.flatnonzero(y <= 1.05 * full)
        if not len(hit):
            m = np.nan
        elif hit[0] == 0:
            m = x[0]
        else:
            j = hit[0]  # log-linear interpolation between the last level above and the first within 5 %
            t = (y[j - 1] - 1.05 * full) / (y[j - 1] - y[j])
            m = float(np.exp(np.log(x[j - 1]) + t * (np.log(x[j]) - np.log(x[j - 1]))))
        m95.append(dict(metric=metric, m95=m, err_all=full))
    pd.DataFrame(m95).to_csv(OUT / f"m95{tag}.csv", index=False)

    fam = raw[raw.arm.isin(["cells", "thin"])].groupby(["arm", "level", "compartment"]).agg(
        err_offset_h=("d_off", "median"), err_c=("d_c", "median")).reset_index()
    fam.to_csv(OUT / f"families{tag}.csv", index=False)
    pd.set_option("display.width", 200)
    print(E.round(4).to_string(index=False))
    print(S.round(3).to_string(index=False))
    print(pd.DataFrame(m95).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
