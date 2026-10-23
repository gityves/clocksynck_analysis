"""Shared reader of the skin single-cell pseudobulk tables (psoriasis, healthy-skin differentiation bins, body sites).

The tables are the processed output of the legacy raw-data builders (`paired_damping/scripts/build_psoriasis_sc.py`,
`build_kc_states.py`, `build_healthy_skin_bins.py`, `build_bins_both.py`, `build_skin_sites.py`): QC, a marker-and-Leiden
class rule (random_state 0), keratinocyte differentiation quartile bins and raw count sums per sample x class. Those steps
need the GEO / Census / ArrayExpress downloads, scanpy clustering and hours of CPU and are NOT re-implemented here (local:
the standard-table step only). Everything after them is: counts of the 14 clock genes, `n_cells`, `depth` (all-gene counts)
per unit x status x class. The all-gene sums (`*_all.h5ad`, used by A60 / A61 for class signatures) are read by `allgenes()`.

Standard table: one row per unit x status x class: `donor` (unit key), `compartment` (class), `status`
(lesional / nonlesional / healthy), `stratum` (study), `n_cells`, `depth`, the 14 clock genes (`ARNTL` named as in the
panel; the files may carry `BMAL1`). Donor table: one row per unit with its study, patient, visit (and site, volar).
"""

from __future__ import annotations

import pandas as pd

from common.paths import DATA

CLOCK = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"]


def read(path, key=("study", "patient", "visit"), donor_cols=("study", "patient", "visit")):
    """(pb, donors) of a processed pseudobulk csv. `key`: columns joined by `|` into the unit key `donor`."""
    d = pd.read_csv(path)
    if "ARNTL" not in d:
        d = d.rename(columns={"BMAL1": "ARNTL"})
    d["donor"] = d[list(key)].astype(str).agg("|".join, axis=1)
    d["stratum"] = d.study.astype(str)
    cols = ["donor", "compartment", "status", "stratum", "sample", "n_cells", "depth"] + CLOCK
    pb = d[cols].reset_index(drop=True)
    donors = d.drop_duplicates("donor").set_index("donor")[[c for c in donor_cols if c in d] +
                                                           [c for c in ("site", "volar") if c in d]]
    return pb, donors


def allgenes(name: str):
    """AnnData of all-gene count sums (rows as the standard table: obs study, patient, visit, status, sample,
    compartment, n_cells, depth) for `psoriasis_sc/<name>.h5ad`."""
    import anndata as ad
    return ad.read_h5ad(DATA / "psoriasis_sc" / f"{name}.h5ad")
