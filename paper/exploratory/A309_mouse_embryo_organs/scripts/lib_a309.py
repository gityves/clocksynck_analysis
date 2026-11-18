"""A309: the human fetal organ pattern in the late mouse embryo (Qiu 2024, E15.0-E18.75). Shared by power and scoring.

Registered in paper/docs/registered_experiments.md (A309). Rows: pseudobulks with >= 50 cells of five cell types in the
16 late-window embryos; block = batch (embryonic day; one embryo at 00, 06, 12, 18 h in each); covariates = batch
indicators; permutations within batch.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import data  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 309
WF = cs.gtex()
GENES = list(WF.genes)
SI = ["PER1", "PER2", "NR1D1"]
CORE = [g for g in GENES if g not in SI]
GC = ["FKBP5", "TSC22D3", "ZBTB16", "KLF9"]
QIU = DATA / "natexp" / "qiu2024"
GENOME = QIU / "genome_late" / "pseudobulk_celltype.h5ad"
READOUTS = {"adrenal": "Adrenocortical cells", "liver": "Hepatocytes", "lung": "Lung progenitor cells",
            "kidney": "Nephron progenitors", "heart": "Ventricular cardiomyocytes"}
MIN_CELLS = 50
LATE = (15.0, 19.0)


def cov(g):
    """Batch indicators (the embryonic day)."""
    return dev.design(g, None, 0, ("block",))


def late_rows() -> pd.DataFrame:
    """Rows of the five readouts in the late window (14-gene table): donor, compartment, block, day, t_clock, depth, genes."""
    pb, don = data.load("mouse_qiu")
    pb = pb.join(don[["day", "t_clock", "batch"]], on="donor")
    keep = pb.compartment.isin(READOUTS.values()) & (pb.n_cells >= MIN_CELLS) & pb.day.between(*LATE, inclusive="left")
    x = pb[keep].copy()
    x["block"] = x.batch
    x["readout"] = x.compartment.map({v: k for k, v in READOUTS.items()})
    return x.reset_index(drop=True)


def genome(rows: pd.DataFrame) -> pd.DataFrame:
    """Genome-wide counts of `rows` (index aligned), columns = upper-cased symbols (clock genes by the panel name; of
    duplicated symbols the most expressed kept); clock counts checked against the 14-gene table."""
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(GENOME)
    key = pd.MultiIndex.from_arrays([a.obs.donor.astype(str), a.obs.compartment.astype(str)])
    want = pd.MultiIndex.from_arrays([rows.donor.astype(str), rows.compartment.astype(str)])
    i = pd.Series(np.arange(a.n_obs), index=key).reindex(want)
    assert i.notna().all()
    X = a.X[i.to_numpy(int)]
    X = np.asarray(X.toarray() if sp.issparse(X) else X, float)
    panel = pd.read_csv(QIU / "prep" / "panel.csv")
    clock = {m: canonical(f) for f, m in zip(panel.feature_name, panel.mouse_symbol) if canonical(f) in GENES}
    sym = np.array([clock.get(x, str(x).upper()) for x in a.var["symbol"].astype(str)])
    order = np.argsort(-X.mean(0), kind="stable")
    keep = np.sort(order[~pd.Series(sym[order]).duplicated().to_numpy()])
    out = pd.DataFrame(X[:, keep], columns=sym[keep], index=rows.index)
    assert np.array_equal(out[GENES].to_numpy(float), rows[GENES].to_numpy(float)), "clock counts differ"
    return out
