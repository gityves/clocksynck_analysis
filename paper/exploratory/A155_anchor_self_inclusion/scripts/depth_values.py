"""Depth values of the A154 / A155 simulations, read from the real tables (no controlled per-donor rows are written).

Single cell: median clock-gene UMI per (donor, compartment) pseudobulk, pooled over the A17 single-cell pseudobulk tables
present in `$SCCC_DATA/pseudobulk` that the A17 anchor table lists. Bulk: median clock-gene count per GTEx sample of the
fast deaths with a time of death (aggregate only). The simulator's `depth` (genome-wide count) is then set so that the
expected clock-gene UMI of a pseudobulk equals that median (`depth_for`).

    python depth_values.py      -> ../output/depth_values.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common.paths import PSEUDOBULK  # noqa: E402

OUT = HERE.parent / "output" / "depth_values.csv"
MESOR_SD = 0.5  # simulate.Compartment default: E exp(offset) = exp(sd^2 / 2)


def depth_for(clock_umi: float) -> float:
    """Genome-wide depth for `simulate` giving `clock_umi` expected clock-gene counts per pseudobulk."""
    w = cs.gtex()
    return float(clock_umi / (np.exp(np.asarray(w.mesor)).sum() * np.exp(MESOR_SD ** 2 / 2)))


def load() -> dict:
    return pd.read_csv(OUT).set_index("depth").clock_umi.to_dict()


if __name__ == "__main__":
    G = list(cs.gtex().genes)
    anchors = pd.read_csv(PAPER / "common" / "anchors_a17.csv")
    rows = []
    for n in sorted(set(anchors.target)):
        f = PSEUDOBULK / f"{n}.parquet"
        if f.exists():
            d = pd.read_parquet(f)
            rows.append(pd.DataFrame({"table": n, "clock_umi": d[G].sum(axis=1).to_numpy(), "depth": d.depth.to_numpy()}))
    sc = pd.concat(rows)
    from common import gtex
    pb, _ = gtex.load()
    bulk = pb[G].sum(axis=1)
    res = pd.DataFrame({"depth": ["single_cell", "bulk"], "clock_umi": [sc.clock_umi.median(), bulk.median()],
                        "genome_depth_median": [sc.depth.median(), pb.depth.median()],
                        "n_rows": [len(sc), len(pb)], "n_tables": [sc.table.nunique(), np.nan]})
    res["sim_depth"] = [depth_for(x) for x in res.clock_umi]
    OUT.parent.mkdir(exist_ok=True)
    res.to_csv(OUT, index=False)
    print(res.to_string())
