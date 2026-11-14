"""Scout for roadmap item T3 (rotation between two same-day biopsies, program level). Simulation only.

The arms with nothing done between biopsies (design inventory 2026-11-15, `README.md`): GSE40589 saline arm (adipose,
7 persons, 4.0 h apart, clock time of the first biopsy not recorded, 08:00 assumed) and GSE59088 control arm (muscle, 6
persons, 2.5 h apart, first biopsy 13:00). A program of G genes with reference amplitudes A_g and acrophases phi_g (timed
muscle, GSE108539; `../T1_nightshift/output/reference_genes.csv`, p < 0.01) rotating with the clock changes by
delta_g = A_g [cos w(t2 - phi_g) - cos w(t1 - phi_g)]. Statistic: per person b = sum_g delta_g d_g / sum_g delta_g^2
(1 = rotates as the reference, 0 = flat); one-sample t over persons. Noise: within-person difference SD per gene of
0.32 log2 (the median residual SD of the GSE106800 muscle differences, T1 scout), genes independent (optimistic: real
genes share noise). Power to reject b = 0 (one-sided 0.05) for a program that rotates at the reference amplitude, by
program size, in each arm and with both arms pooled (13 persons, assumes the muscle reference also holds in fat).

    python scout_t3.py
Output: output/power.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
W = 2 * np.pi / 24
SD = 0.32
ARMS = {"GSE40589 saline (adipose)": (7, 8.0, 12.0), "GSE59088 control (muscle)": (6, 13.0, 15.5)}


def main():
    ref = pd.read_csv(HERE.parent / "T1_nightshift/output/reference_genes.csv")
    ref = ref[ref.p < 0.01].reset_index(drop=True)
    rng = np.random.default_rng(1115)
    rows = []
    for G in (10, 20, 50, 100, 200):
        for arm in list(ARMS) + ["pooled"]:
            hits = []
            for _ in range(2000):
                idx = rng.choice(len(ref), G, replace=False)
                a, phi = ref.amp.to_numpy()[idx], ref.acro_h.to_numpy()[idx]
                b = []
                for name, (n, t1, t2) in ARMS.items():
                    if arm not in (name, "pooled"):
                        continue
                    dl = a * (np.cos(W * (t2 - phi)) - np.cos(W * (t1 - phi)))
                    d = dl[None, :] + rng.normal(0, SD, (n, G))
                    b += list(d @ dl / (dl @ dl))
                b = np.asarray(b)
                t = b.mean() / (b.std(ddof=1) / np.sqrt(len(b)))
                hits.append(stats.t.sf(t, len(b) - 1) < 0.05)
            rows.append(dict(G=G, arm=arm, persons=sum(ARMS[k][0] for k in ARMS) if arm == "pooled" else ARMS[arm][0],
                             power_rotation_vs_flat=float(np.mean(hits))))
    out = pd.DataFrame(rows)
    (HERE / "output").mkdir(exist_ok=True)
    out.to_csv(HERE / "output/power.csv", index=False)
    print(out.pivot(index="G", columns="arm", values="power_rotation_vs_flat").round(3).to_string())


if __name__ == "__main__":
    main()
