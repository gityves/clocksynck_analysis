"""A271 Part B, post hoc (not registered; chosen after H271.2 was scored): is "stronger clocks are more coupled" (A104) a
count artefact, and which variance term carries it?

Per GTEx tissue of the A104 table (`results/figure_2/output/backbone/a104_gtex_tissues.csv`, fast deaths, >= 40 donors):
median over donors of the summed clock-gene counts, of the smallest clock-gene count of a row and of the library depth
(rows as A104, `gtex.load(min_donors=40)`). Per block, Spearman rho of own amplitude A = sqrt(z2 V + 2 psi) and of coupling
c with each count summary and with the shared (z2 V) and own (2 psi) terms, and the partial rho(A, c | clock counts) on
ranks. Per-tissue aggregates only (GTEx controlled access).

    python posthoc_a104_counts.py
Outputs (`../output/`): a104_counts_tissues.csv (one row per tissue), a104_counts_rho.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
A104 = HERE.parents[3] / "results" / "figure_2" / "output" / "backbone" / "a104_gtex_tissues.csv"


def partial_rho(x, y, z):
    """Spearman partial correlation of x and y given z (Pearson on rank residuals)."""
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    res = [r - np.polyval(np.polyfit(rz, r, 1), rz) for r in (rx, ry)]
    return float(np.corrcoef(*res)[0, 1])


def main():
    pb, _ = gtex.load(min_donors=40)
    genes = list(cs.gtex().genes)
    x = pb[genes].to_numpy(float)
    rows = pd.DataFrame(dict(compartment=pb.compartment, clock_counts=x.sum(1), min_gene_count=x.min(1), depth=pb.depth))
    counts = rows.groupby("compartment").median()
    t = pd.read_csv(A104).set_index("tissue").join(counts, how="inner")
    t["shared_z2V"], t["own_2psi"] = t.z2 * t.V, 2 * t.psi
    OUT.mkdir(exist_ok=True)
    t.to_csv(OUT / "a104_counts_tissues.csv")
    out = []
    for block, g in t.groupby("block"):
        g = g.dropna(subset=["own_amplitude", "coupling"])
        for a in ("own_amplitude", "coupling"):
            for b in ("shared_z2V", "own_2psi", "clock_counts", "min_gene_count", "depth"):
                out.append(dict(block=block, x=a, y=b, n=len(g), rho=float(spearmanr(g[a], g[b]).statistic)))
        out.append(dict(block=block, x="own_amplitude", y="coupling | clock_counts", n=len(g),
                        rho=partial_rho(g.own_amplitude, g.coupling, g.clock_counts)))
    R = pd.DataFrame(out)
    R.to_csv(OUT / "a104_counts_rho.csv", index=False)
    print(R.to_string(index=False))
    print(t.groupby("block")[["shared_z2V", "own_2psi", "clock_counts"]].agg(["min", "max"]).T.to_string())


if __name__ == "__main__":
    main()
