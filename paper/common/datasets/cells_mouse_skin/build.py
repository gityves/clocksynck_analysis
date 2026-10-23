"""Per-cell clock-gene counts of the mouse dermis time course (GSE223109, 18 samples, ZT 2-22 every 4 h): the positive control
of the per-cell amplitude vs synchrony split (A18 part 2; V10).

Source:   `skin.h5ad` of the CoPhaser training set (`$SCCC_SKIN_H5AD`, default
          /mnt/max8tb/data/CoPhaser_data/training_set/circadian_cycle/skin.h5ad; raw counts in `layers["counts"]`, `sample_name`,
          `ZTmod`, `annotation`), read in row chunks (backed): the whole matrix is never in memory.
Output:   `$SCCC_DATA/cells/mouse_skin_gse223109/cells.parquet`: per cell the 14 clock-gene counts (GTEx order), `umi`
          (genome-wide counts), `sample`, `zt`, `ct` (annotation); cells with `umi` >= 2,000 only (as V10).
Replaces: `load()` of legacy/scripts/v10_cell_dispersion.py.
Check:    `python build.py check`: equals `v10_cell_dispersion.load()` on the first 20,000 cells passing the filter (run 2026-10-22).

    python build.py        # builds the table if missing and prints its size
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from clocksync import gtex  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402
from common.paths import DATA  # noqa: E402

SKIN = Path(os.environ.get("SCCC_SKIN_H5AD", "/mnt/max8tb/data/CoPhaser_data/training_set/circadian_cycle/skin.h5ad"))
OUT = DATA / "cells" / "mouse_skin_gse223109" / "cells.parquet"
GENES = list(gtex().genes)
MIN_UMI = 2000


def build(chunk: int = 50_000, max_cells: int | None = None) -> pd.DataFrame:
    import anndata as ad
    a = ad.read_h5ad(SKIN, backed="r")
    cols = gene_columns(a.var_names, GENES)
    if (cols < 0).any():
        raise KeyError("clock genes missing from the skin table")
    obs = a.obs
    parts = []
    n = a.n_obs if max_cells is None else min(a.n_obs, max_cells)
    for i in range(0, n, chunk):
        sl = slice(i, min(i + chunk, n))
        X = sp.csr_matrix(a.layers["counts"][sl])
        umi = np.asarray(X.sum(1)).ravel()
        ok = umi >= MIN_UMI
        C = X[ok][:, cols].toarray().astype(np.int32)
        o = obs.iloc[sl][ok]
        parts.append(pd.DataFrame(C, columns=GENES).assign(
            umi=umi[ok], sample=o.sample_name.astype(str).to_numpy(), zt=o.ZTmod.astype(float).to_numpy(),
            ct=o.annotation.astype(str).to_numpy()))
    return pd.concat(parts, ignore_index=True)


sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/
from common.paths import ARCHIVE  # noqa: E402
def load() -> pd.DataFrame:
    if not OUT.exists():
        OUT.parent.mkdir(parents=True, exist_ok=True)
        build().to_parquet(OUT, index=False)
    return pd.read_parquet(OUT)


def check(n: int = 20_000) -> None:
    sys.path.insert(0, str(ARCHIVE / "legacy" / "scripts"))
    new = build(max_cells=n)
    import anndata as ad
    a = ad.read_h5ad(SKIN, backed="r")
    X = sp.csr_matrix(a.layers["counts"][:n])
    umi = np.asarray(X.sum(1)).ravel()
    ok = umi >= MIN_UMI
    cols = gene_columns(a.var_names, GENES)
    C = X[ok][:, cols].toarray()
    print("cells", len(new), int(ok.sum()), "counts equal:", bool((new[GENES].to_numpy() == C).all()),
          "umi equal:", bool(np.allclose(new.umi, umi[ok])))


if __name__ == "__main__":
    check() if sys.argv[1:] == ["check"] else print(len(load()))
