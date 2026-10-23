"""Local helpers of method. Everything marked "local: not in clocksync" is a short statistic the package has no public call for."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import gtex  # noqa: E402

W = 2 * np.pi / 24
HOURS = 1 / W
WF = cs.gtex()
CLOCK = list(WF.genes)


def timed_loading(beta: np.ndarray, hours, groups=None) -> complex:
    """local: not in clocksync (`timed.cosinor` is for real series). Complex loading z of the clock-plane state on clock time,
    beta(t) = mu(group) + z e^{iWt} (the conformal part of the regression on (cos, sin) of the time), arg z = acrophase in
    the clock plane."""
    t = np.asarray(hours, float) * W
    D = [np.ones((len(t), 1))] if groups is None else [pd.get_dummies(pd.Series(groups).astype(str)).to_numpy(float)]
    D = np.column_stack(D + [np.cos(t), np.sin(t)])
    B, *_ = np.linalg.lstsq(D, np.asarray(beta, float), rcond=None)
    (cx, cy), (sx, sy) = B[-2], B[-1]
    return complex((cx + sy) / 2, (cy - sx) / 2)


def cdiff(a, b):
    """Circular difference a - b in hours, in [-12, 12)."""
    return cs.stats.wrap(np.asarray(a, float) - np.asarray(b, float))


def load(tissues, genes="all") -> dict[str, tuple[pd.DataFrame, pd.Series]]:
    """{tissue: (counts: donors x gene symbols, time of death h)} of the fast deaths (`common.gtex.load`; duplicate symbols summed)."""
    pb, donors = gtex.load(tissues=tissues, genes=genes, depth="sum")
    meta = {"donor", "compartment", "depth"}
    out = {}
    for t, g in pb.groupby("compartment"):
        c = g.set_index("donor")[[x for x in g.columns if x not in meta]]
        if c.columns.duplicated().any():
            c = c.T.groupby(level=0).sum().T
        out[t] = (c, donors.tod.reindex(c.index).rename("tod"))
    return {t: out[t] for t in tissues if t in out}


def donor_times(data) -> pd.Series:
    return pd.concat([tod for _, tod in data.values()]).groupby(level=0).first()


