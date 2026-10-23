"""Clock correlation distance (CCD; Shilts et al. 2018, PeerJ 6:e4327), as in the deltaccd R package.

A group of samples is scored by the Euclidean distance between its Spearman correlation matrix of 12 clock genes
and a pan-tissue reference (8 mouse circadian datasets; `ccd_reference.csv`), over the upper triangle. Large CCD
is read as a disorganised clock; delta-CCD is the CCD of a test group minus that of a control group. The p-value
compares the CCD to that of random gene sets of the same size (deltaccd's gene permutation).
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


@cache
def reference() -> pd.DataFrame:
    return pd.read_csv(Path(__file__).with_name("ccd_reference.csv"), index_col=0, comment="#")


def _distance(ref: np.ndarray, X: np.ndarray) -> float:
    iu = np.triu_indices(ref.shape[0], 1)
    r = spearmanr(X).statistic
    return float(np.sqrt(np.nansum((ref[iu] - r[iu]) ** 2)))


def ccd(logexpr: pd.DataFrame, n_perm: int = 0, seed: int = 0) -> tuple[float, float]:
    """`logexpr`: samples x genes of one group, ideally all measured genes (for the permutation null).
    Returns (CCD, p-value; NaN without permutations)."""
    ref = reference()
    genes = [g for g in ref.index if g in logexpr.columns]
    obs = _distance(ref.loc[genes, genes].to_numpy(), logexpr[genes].to_numpy(float))
    if n_perm < 1:
        return obs, np.nan
    rng = np.random.default_rng(seed)
    X = logexpr.loc[:, logexpr.std() > 0].to_numpy(float)
    R = ref.loc[genes, genes].to_numpy()
    null = np.array([_distance(R, X[:, np.sort(rng.choice(X.shape[1], len(genes), replace=False))])
                     for _ in range(n_perm)])
    return obs, float((1 + (null <= obs).sum()) / (n_perm + 1))
