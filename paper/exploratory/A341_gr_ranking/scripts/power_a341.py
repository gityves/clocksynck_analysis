"""A341 power stage (registered in `paper/docs/registered_experiments.md`): detectable Spearman rho across k organs.

One-sided permutation Spearman at alpha 0.05 (exact over all k! orders); power under a bivariate normal latent with
correlation r between outcome and predictor (20,000 draws per r). The detectable rho is the smallest r with power >= 0.8.
If it exceeds 0.8, H341 is reported descriptively only (registered rule).

    python power_a341.py
Output: ../output/power.csv
"""
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT = Path(__file__).resolve().parents[1] / "output"


def critical(k: int) -> float:
    base = np.arange(k)
    null = np.array([spearmanr(base, p)[0] for p in permutations(base)])
    return float(np.quantile(null, 0.95, method="higher")), null


def main():
    rng = np.random.default_rng(341)
    rows = []
    for k in (6, 7, 8):  # human organs (7 registered units: adrenal cortex, medulla, lung, liver, heart, kidney, brain)
        crit, null = critical(k)
        for r in np.arange(0.5, 1.0, 0.05):
            x = rng.multivariate_normal([0, 0], [[1, r], [r, 1]], size=(20000, k))
            rho = np.array([spearmanr(a[:, 0], a[:, 1])[0] for a in x])
            p = np.array([(np.sum(null >= v) / len(null)) for v in rho])
            rows.append(dict(k=k, r=round(r, 2), critical_rho=crit, power=float(np.mean(p < 0.05))))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "power.csv", index=False)
    print(t.pivot(index="r", columns="k", values="power").round(3).to_string())
    for k, g in t.groupby("k"):
        ok = g[g.power >= 0.8]
        print(k, "organs: critical rho", round(g.critical_rho.iloc[0], 3), "| detectable r",
              ok.r.min() if len(ok) else "> 0.95")


if __name__ == "__main__":
    main()
