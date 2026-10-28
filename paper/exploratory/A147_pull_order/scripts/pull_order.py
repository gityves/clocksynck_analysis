"""Draft clocksync readout `pull_order` (A147): pull versus order of compartment-pair phases against a reference group.

Move to `clocksync` once A147 passes. Inputs are Hermitian cross-covariances E[z_a conj z_b] (compartment x compartment)
of the reference group (e.g. fast deaths) and of the test group, as returned by the clocksync pseudobulk moments.

    deviation(theta)       nearest pole (0 or 12 h) and signed deviation from it (h) of pair phases (rad)
    pull_order(C_ref, C_grp, pairs=None, n_perm=1000, rng=None, min_dev_h=1.0)
        order  Spearman of the group's against the reference's signed deviations over pairs;
               one-sided p from `n_perm` joint permutations of the compartment labels in the group
        ratio  median of group / reference deviation over same-pole pairs with |reference deviation| >= min_dev_h
               (1 = no pull, 0 = full pull)
"""

from __future__ import annotations

import numpy as np
from scipy.stats import rankdata, spearmanr

HOURS = 24.0 / (2 * np.pi)


def deviation(theta: np.ndarray, period_h: float = 24.0):
    """Nearest pole (0 or period/2) and signed deviation (h) of pair phases theta (rad)."""
    h = np.asarray(theta) * period_h / (2 * np.pi)
    half = period_h / 2
    pole = np.abs(h) > half / 2
    return pole, np.where(pole, h - np.sign(h) * half, h)


def _pair_index(T: int, pairs):
    if pairs is None:
        return np.triu_indices(T, 1)
    pairs = np.asarray(pairs)
    return pairs[:, 0], pairs[:, 1]


def pull_order(C_ref, C_grp, pairs=None, n_perm: int = 1000, rng=None, min_dev_h: float = 1.0) -> dict:
    """Order (Spearman, permutation p) and shrink ratio of pair-phase deviations from the poles, group vs reference."""
    C_ref, C_grp = np.asarray(C_ref), np.asarray(C_grp)
    T = C_ref.shape[0]
    ia, ib = _pair_index(T, pairs)
    pf, df = deviation(np.angle(C_ref[ia, ib]))
    pg, dg = deviation(np.angle(C_grp[ia, ib]))
    big = (pf == pg) & (np.abs(df) >= min_dev_h)
    out = dict(order=float(spearmanr(df, dg)[0]),
               ratio=float(np.median(dg[big] / df[big])) if big.any() else np.nan,
               same_pole=float((pf == pg).mean()), n_pairs=int(len(df)), n_ratio=int(big.sum()))
    if n_perm:
        rng = np.random.default_rng() if rng is None else rng
        perms = np.argsort(rng.random((n_perm, T)), axis=1)
        dev = deviation(np.angle(C_grp[perms[:, ia], perms[:, ib]]))[1]  # (n_perm, pairs)
        rf = rankdata(df) - rankdata(df).mean()
        rg = np.apply_along_axis(rankdata, 1, dev)
        rg = rg - rg.mean(1, keepdims=True)
        null = (rg @ rf) / np.sqrt((rg ** 2).sum(1) * (rf ** 2).sum())
        out["p_order"] = float((1 + (null >= out["order"]).sum()) / (1 + n_perm))
    return out
