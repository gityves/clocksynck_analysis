"""CHIRAL (Talamanca, Gobet and Naef 2023, Science 379:478) through its Python port, pychiral (A. Salati,
github.com/AndreaSalati/pychiral). CHIRAL infers one phase per sample from a few clock genes: a mean-field
(spin-glass) initialisation on the sample-sample correlation matrix, then EM on a harmonic regression with a
two-state (rhythmic or flat) model per gene. Inputs are log expression, standardised per gene by CHIRAL.
"""

from __future__ import annotations

import contextlib
import io

import numpy as np
import pandas as pd

# pychiral's default clock genes (mouse symbols in the port; Bmal1 is Arntl), as human symbols
CHIRAL_GENES = ("ARNTL", "NPAS2", "CRY1", "CRY2", "PER1", "PER2", "NR1D1", "NR1D2", "TEF", "DBP", "CIART", "PER3")


def run_chiral(logexpr: pd.DataFrame, genes=CHIRAL_GENES, seed: int = 0, attempts: int = 10) -> tuple[pd.Series, int]:
    """`logexpr`: samples x genes (log scale). Returns phase (h, arbitrary origin and direction) per sample and the
    number of failed starts: the EM's Lagrange step can hit a singular system, and the run is then restarted from a
    new random mean-field start (NaN phases if all `attempts` fail)."""
    from pychiral import CHIRAL  # the mean-field start draws from numpy's global generator

    X = logexpr.loc[:, [g for g in genes if g in logexpr.columns]]
    X = X.loc[:, X.std() > 0]
    for k in range(attempts):
        np.random.seed(seed * attempts + k)
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                phi = np.asarray(CHIRAL(X.to_numpy(float), clockgenes=list(X.columns))["phi"])
            return pd.Series(phi * 24 / (2 * np.pi) % 24, index=logexpr.index, name="phase"), k
        except np.linalg.LinAlgError:
            continue
    return pd.Series(np.nan, index=logexpr.index, name="phase"), attempts
