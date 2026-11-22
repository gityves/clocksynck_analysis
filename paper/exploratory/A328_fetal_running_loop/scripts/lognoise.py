"""Counting noise of log(count + 0.5) for A339 / A340.

The delta-method variance 1 / (count + 0.5) overstates the variance of log(X + 0.5) when X ~ Poisson(lam) is small
(the pseudocount compresses low counts), which over-corrects noise at fetal intron depths. `var_log` gives the exact
variance as a function of the expected count lam (interpolated on a log grid), and `expected` the expected counts
from cohort median fractions x depth.
"""
import numpy as np
from scipy.stats import poisson

_GRID = np.logspace(-3, 4.5, 400)


def _exact(lam: float) -> float:
    k = np.arange(int(lam + 12 * np.sqrt(lam) + 30))
    p = poisson.pmf(k, lam)
    y = np.log(k + 0.5)
    m = np.sum(p * y)
    return float(np.sum(p * (y - m) ** 2))


_VAR = np.array([_exact(x) for x in _GRID])


def var_log(lam) -> np.ndarray:
    """Var[log(X + 0.5)], X ~ Poisson(lam); 1 / lam beyond the grid."""
    lam = np.asarray(lam, float)
    v = np.interp(np.log(np.clip(lam, _GRID[0], None)), np.log(_GRID), _VAR)
    return np.where(lam > _GRID[-1], 1.0 / np.maximum(lam, 1e-12), v)


def expected(x, depth, groups) -> np.ndarray:
    """Expected counts (n x G): cohort median fraction of each gene x the row's depth (depth n or n x G)."""
    x = np.asarray(x, float)
    depth = np.broadcast_to(np.asarray(depth, float).reshape(len(x), -1), x.shape)
    out = np.empty_like(x)
    for g in np.unique(groups):
        m = groups == g
        out[m] = np.median(x[m] / depth[m], axis=0)[None] * depth[m]
    return out
