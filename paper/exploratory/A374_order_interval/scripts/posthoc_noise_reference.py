"""A374 post hoc (2026-11-24, descriptive; on the user's request for the slide 3.last figure): how similar to the fast-death
arrangement would a group of n fast deaths be, from sampling noise alone?

n fast deaths drawn without replacement (n = each group's donor count) against the remaining fast deaths as reference;
Fisher-Lee rho of the 66 pair phases (pairs with >= 20 donors in both); 200 draws per n (seed [374, 2]). 14 genes and core
loop. Aggregates only.

    python posthoc_noise_reference.py
Output: ../output/noise_reference.csv
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a374_interval as A  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

NS = {"slow deaths (Hardy 4)": 103, "Hardy 3 (descriptive)": 50, "ventilated, living brain": 79}


def main():
    limit_threads(4)
    T, fast, _ = A.C.tissues()
    rng = np.random.default_rng([374, 2])
    rows = []
    for gs, kind in (("14 genes", "clock"), ("core", "core")):
        Z, M = A.S.raw(T, kind, 0, T.R)
        for g, n in NS.items():
            r = []
            for _ in range(200):
                p = rng.permutation(fast)
                sub, ref = p[:n], p[n:]
                mask = A.B.pair_mask(M, ref, sub)
                r.append(float(ST.fisher_lee(A.phases(Z, M, ref)[mask], A.phases(Z, M, sub)[mask])))
            rows.append(dict(gene_set=gs, group=g, n=n, rho_median=float(np.median(r)), rho_q05=float(np.quantile(r, .05)),
                             rho_q95=float(np.quantile(r, .95))))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(A.OUT / "noise_reference.csv", index=False)


if __name__ == "__main__":
    main()
