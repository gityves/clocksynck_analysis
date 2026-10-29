"""Bounded Poisson-corrected |complex correlation| of paired clock-plane states (descriptive / post hoc, 2026-10-29).

Wraps the public `clocksync.stats.ccorr` (`return_flag`, `nan_on_flag`): NaN plus a flag when a corrected variance is
<= CCORR_K standard errors. The old statistic (`q3_lib.corr_stat`) is the same formula without the flag, with only a
non-positive corrected variance giving NaN. A flagged unit is NaN in observed value and null, so pooled means drop it
(weights renormalised, `q3_lib.pool`).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER / "exploratory" / "A138_cancer_estimability" / "scripts"))
import q3_lib as q  # noqa: E402,F401
from clocksync import stats as S  # noqa: E402


def _beta(z):
    z = np.asarray(z)
    return np.column_stack([z.real, z.imag])


def stat(za, vA, zb, vB) -> dict:
    """Old |ccorr| (`old`), bounded |ccorr| (`new`, NaN when flagged) and the flag, from centred complex states."""
    ba, bb = _beta(za), _beta(zb)
    val, flag = S.ccorr(ba, bb, vA, vB, return_flag=True)
    new = S.ccorr(ba, bb, vA, vB, nan_on_flag=True)
    return dict(old=float(abs(val)), new=float(abs(new)), flag=bool(flag))


def new_stat(u) -> float:
    return stat(u["za"], u["vA"], u["zb"], u["vB"])["new"]


def perm_null(u, n_perm, rng) -> np.ndarray:
    """Permutation null of the bounded statistic: the registered null (same rng stream) or all NaN when flagged."""
    nl = q.perm_null(u["za"], u["vA"], u["zb"], u["vB"], n_perm, rng)
    return np.full(n_perm, np.nan) if stat(u["za"], u["vA"], u["zb"], u["vB"])["flag"] else nl


def estimable(u, null) -> bool:
    """Registered rule (corrected variances > 0, null mean <= 1) plus not flagged."""
    return bool(np.isfinite(null).any() and q.estimable(u["za"], u["vA"], u["zb"], u["vB"], null))


def perm_p(null, obs) -> float:
    x = null[np.isfinite(null)]
    return float((1 + np.sum(x >= obs)) / (1 + len(x)))
