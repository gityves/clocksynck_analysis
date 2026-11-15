"""Scout of A265 (no clock statistic, no COFE output on real data): the design per cancer and the COFE runtime.

Design: tumours of the S2 Table histology, pairs, genes after COFE's filters. Runtime: one COFE run as published
(`common.benchmark.cofe.run_cofe`) on the largest cancer with every gene permuted on its own across the tumours (no
structure left), so that no real COFE label is produced before the registration.

    python scout_a265.py [--jobs 8]
Outputs (`../output/`): scout_design.csv, scout_runtime.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import lib_a265 as L  # noqa: E402

OUT = HERE.parents[1] / "output"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rows, big = [], None
    for c in L.CANCERS:
        d = L.cancer(c)
        rows.append(dict(cancer=c, tumours=len(d["cofe_t"]), pairs=len(d["pairs"]), genes_cofe=d["cofe_t"].shape[1],
                         genes_protein_coding=d["counts_t"].shape[1]))
        print(rows[-1], flush=True)
        if big is None or len(d["cofe_t"]) > len(big):
            big = d["cofe_t"]
    pd.DataFrame(rows).to_csv(OUT / "scout_design.csv", index=False)
    from common.benchmark.cofe import run_cofe
    rng = np.random.default_rng(265)
    x = pd.DataFrame(np.column_stack([rng.permutation(v) for v in big.to_numpy().T]), index=big.index, columns=big.columns)
    t0 = time.time()
    r = run_cofe(x, seed=0, ncores=a.jobs)
    pd.DataFrame([dict(tumours=len(x), genes=x.shape[1], jobs=a.jobs, seconds=time.time() - t0,
                       interior_min=r["interior_min"])]).to_csv(OUT / "scout_runtime.csv", index=False)
    print(f"COFE on permuted data {x.shape}: {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
