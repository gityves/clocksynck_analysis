"""Slide data (results, drives, 4.2): human muscle after 3.5 days of 12 h-shifted behaviour (Wefers 2018, GSE106800).

Per subject, delta = z(19:00) - z(07:00) in the aligned (dc) and misaligned (dm) condition, as A257 / A335 part 1a
(same loading, projector and genes; 14 clock genes and the core loop). Summary: the complex ratio
r = sum dm conj(dc) / sum |dc|^2, whose modulus is the excursion ratio (amplitude of the misaligned 07:00-to-19:00
change relative to the aligned one) and whose argument is the A257 shift s (> 0: the misaligned clock state is further
along its cycle, i.e. peaks earlier). The registered R (A257 H257.3) projects dm on the aligned chord and so contains the
rotation; |r| does not (correction note of 2026-11-19: about R / cos s). Interval: 90 % subject bootstrap (2000 draws,
seed 4303). Descriptive, post hoc display quantity (user, 2026-11-25). Open data (GEO); per-subject values written.

    python scripts/nightshift_compute.py
Output: presentation/output/nightshift/{summary,per_subject}.csv
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("a335", PAPER / "exploratory/A335_review_batch/scripts/a335_1a_nightshift.py")
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)
A = B.A  # A257 module

OUT = PAPER / "presentation" / "output" / "nightshift"
N_BOOT, SEED = 2000, 4303


def deltas(smp, y, wf):
    d = A.deltas(smp, A.zstates(y, wf))
    both = sorted(s for (s, c) in d if c == "ctl" and (s, "mis") in d)
    return both, np.array([d[(s, "ctl")] for s in both]), np.array([d[(s, "mis")] for s in both])


def ratio(dm, dc):
    return np.sum(dm * np.conj(dc)) / np.sum(abs(dc) ** 2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    smp, y = A.load()
    base = A.WF.subset([x for x in A.WF.genes if x in y.columns])
    core = A.WF.drop(B.DROP)
    rng = np.random.default_rng(SEED)
    rows, subs = [], []
    for name, wf in (("14 genes", base), ("core loop", core)):
        both, dc, dm = deltas(smp, y, wf)
        r = ratio(dm, dc)
        s = float(np.angle(r) / A.W)
        bs = np.array([ratio(dm[i], dc[i]) for i in (rng.integers(0, len(both), len(both)) for _ in range(N_BOOT))])
        sb = (np.angle(bs) / A.W - s + 12) % 24 - 12 + s
        rows.append(dict(gene_set=name, subjects=len(both), excursion_ratio=float(abs(r)),
                         excursion_lo90=float(np.quantile(abs(bs), 0.05)),
                         excursion_hi90=float(np.quantile(abs(bs), 0.95)), shift_h=s,
                         shift_lo90=float(np.quantile(sb, 0.05)), shift_hi90=float(np.quantile(sb, 0.95))))
        subs.append(pd.DataFrame(dict(gene_set=name, subject=both, shift_h=np.angle(dm * np.conj(dc)) / A.W,
                                      excursion_ratio=abs(dm) / abs(dc))))
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "summary.csv", index=False)
    pd.concat(subs).to_csv(OUT / "per_subject.csv", index=False)
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
