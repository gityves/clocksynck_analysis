"""COFE's TCGA analysis (Ananthasubramaniam et al. 2025, Fig. 2D) and checks of its premise.

`consensus_clock` is the 'eigen-core clock': the clock genes' rhythms in each cancer, as complex numbers
A e^{i phase}, are combined across cancers by the rank-1 complex SVD after choosing each cancer's direction of time
(an ordering's direction is arbitrary) to maximise the variance explained (S4 Text C). Phases are then referenced
to ARNTL.
"""

from __future__ import annotations

from itertools import product

import numpy as np
import pandas as pd

W = 2 * np.pi / 24

# The 15 expressed core clock genes of COFE Fig. 2D.
COFE_CLOCK = ("ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "TEF",
              "RORA", "NFIL3", "BHLHE40")

# Closest GTEx waveform per TCGA project (sccc.reference.gtex keys); no GTEx bladder clock, so the global one.
GTEX_MATCH = {"BRCA": "Breast", "KIRC": "Kidney - Cortex", "KIRP": "Kidney - Cortex", "PRAD": "Prostate",
              "OV": "Ovary", "UCEC": "Uterus", "LUAD": "Lung", "COAD": "Colon - Transverse", "THCA": "Thyroid",
              "LIHC": "Liver", "BLCA": None}


def phasors(cos: pd.DataFrame) -> pd.Series:
    """Cosinor table (`evaluate.cosinor`) -> A e^{i w phase}."""
    return cos.amplitude * np.exp(1j * W * cos.phase)


def consensus_clock(Z: pd.DataFrame, anchor: str = "ARNTL") -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Z: genes x groups (cancers) of complex rhythms. Returns per gene the consensus amplitude and phase (h,
    `anchor` at 0), per group its relative strength (max 1), direction and rotation (h), and the variance
    explained by the rank-1 approximation."""
    Z = Z.dropna()
    M0 = Z.to_numpy(complex)
    best = None
    for signs in product((1, -1), repeat=M0.shape[1] - 1):  # first group fixes the direction
        d = np.array((1,) + signs)
        M = np.where(d > 0, M0, M0.conj())
        s = np.linalg.svd(M, compute_uv=False)
        frac = s[0] ** 2 / (s**2).sum()
        if best is None or frac > best[0]:
            best = (frac, d)
    frac, d = best
    M = np.where(d > 0, M0, M0.conj())
    U, s, Vh = np.linalg.svd(M, full_matrices=False)
    u, v = U[:, 0] * s[0], Vh[0].conj()
    rot = v / np.abs(v)
    u = u * np.abs(v).max()  # amplitude in the units of the strongest group
    ref = np.angle(u[list(Z.index).index(anchor)]) if anchor in Z.index else 0.0
    genes = pd.DataFrame({"amplitude": np.abs(u), "phase": ((np.angle(u) - ref) / W) % 24}, index=Z.index)
    groups = pd.DataFrame({"strength": np.abs(v) / np.abs(v).max(), "direction": d,
                           "rotation_h": (np.angle(rot) / W) % 24}, index=Z.columns)
    return genes, groups, float(frac)
