"""Scores for cyclic orderings against recorded times, and the clock-gene arrangement they imply.

An untimed ordering is defined up to a rotation and a reflection. `align` gives it both from the truth (an oracle
orientation, favourable to the method), so every error below is the method's best case.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

W = 2 * np.pi / 24


def wrap(h, period: float = 24.0):
    """Circular difference in hours, in [-period/2, period/2)."""
    return (np.asarray(h) + period / 2) % period - period / 2


def circ_corr(a_h, b_h) -> float:
    """Jammalamadaka-SenGupta circular correlation of two sets of times (h); its sign is the direction."""
    a, b = np.asarray(a_h) * W, np.asarray(b_h) * W
    sa = np.sin(a - np.angle(np.exp(1j * a).mean()))
    sb = np.sin(b - np.angle(np.exp(1j * b).mean()))
    return float((sa * sb).sum() / np.sqrt((sa**2).sum() * (sb**2).sum()))


def align(pred_h, true_h, step: float = 0.05) -> tuple[np.ndarray, float, int]:
    """Rotation and reflection of `pred_h` minimising the median absolute error to `true_h`.
    Returns aligned times, that median error (h) and the direction (+1 or -1)."""
    pred_h, true_h = np.asarray(pred_h, float), np.asarray(true_h, float)
    best = (np.inf, 0.0, 1)
    for d in (1, -1):
        for off in np.arange(0, 24, step):
            e = np.median(np.abs(wrap(d * pred_h + off - true_h)))
            if e < best[0]:
                best = (e, off, d)
    e, off, d = best
    return (d * pred_h + off) % 24, float(e), d


def skill(mae: float, true_h, step: float = 0.05) -> float:
    """1 - mae / (median absolute error of the best constant prediction): 0 for no information beyond the sampling
    window, 1 for perfect times. The constant's error is 6 h over a uniform day but about 2.5 h within a 10 h window,
    so raw errors are not comparable across designs."""
    t = np.asarray(true_h, float)
    const = min(np.median(np.abs(wrap(t - c))) for c in np.arange(0, 24, step))
    return float(1 - mae / const)


def spread(h) -> float:
    """Circular resultant length of times (1: all equal, 0: uniform over the day)."""
    return float(np.abs(np.exp(1j * W * np.asarray(h)).mean()))


def plane_anisotropy(beta, V) -> float:
    """Ratio of the small to the large eigenvalue of the clock-plane covariance across samples, net of the mean
    measurement covariance V. Phases uniform over the day give 1 whatever the amplitudes (when independent of
    phase); a 10 h window gives 0.13 at a constant amplitude (`arc_anisotropy`), more when amplitudes vary but well
    below 1. Free of the mesor, which only moves the centre."""
    S = np.cov(np.asarray(beta, float).T) - np.asarray(V, float).mean(0)
    lo, hi = np.linalg.eigvalsh(S)
    return float(max(lo, 0.0) / hi)


def arc_anisotropy(width_h: float) -> float:
    """`plane_anisotropy` of noise-free samples spread uniformly over `width_h` hours at a constant amplitude."""
    a = np.pi * width_h / 24
    s2 = np.sin(2 * a) / (2 * a)
    var_cos = (1 + s2) / 2 - (np.sin(a) / a) ** 2
    var_sin = (1 - s2) / 2
    return float(min(var_cos, var_sin) / max(var_cos, var_sin))


def cosinor(logexpr: pd.DataFrame, times_h, genes=None) -> pd.DataFrame:
    """Per gene: y = m + A cos(w (t - phase)). Amplitude A in the units of `logexpr` (log2 -> peak-to-trough fold
    2^(2A)), phase in h, and R^2."""
    X = logexpr if genes is None else logexpr.loc[:, [g for g in genes if g in logexpr.columns]]
    t = np.asarray(times_h, float) * W
    D = np.column_stack([np.ones_like(t), np.cos(t), np.sin(t)])
    Y = X.to_numpy(float)
    coef, *_ = np.linalg.lstsq(D, Y, rcond=None)
    res = Y - D @ coef
    r2 = 1 - (res**2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)
    return pd.DataFrame({"amplitude": np.hypot(coef[1], coef[2]),
                         "phase": (np.arctan2(coef[2], coef[1]) / W) % 24, "r2": r2}, index=X.columns)


def arrangement_error(est: pd.DataFrame, truth: pd.DataFrame, reflect: bool = False) -> pd.DataFrame:
    """Gene acrophases of `est` after the amplitude-weighted rotation (and optional reflection) that best matches
    `truth`, i.e. errors of the relative arrangement only. Columns: phase error (h), log2 amplitude ratio."""
    g = est.index.intersection(truth.index)
    e, t = est.loc[g], truth.loc[g]
    w = np.sqrt(e.amplitude * t.amplitude)
    best = None
    for d in ((1, -1) if reflect else (1,)):
        rot = np.angle((w * np.exp(1j * W * (t.phase - d * e.phase))).sum()) / W
        err = wrap(d * e.phase + rot - t.phase)
        if best is None or np.abs(err).mean() < np.abs(best).mean():
            best = err
    return pd.DataFrame({"phase_err": best, "log2_amp_ratio": np.log2(e.amplitude / t.amplitude)}, index=g)


# --- Tissue estimands without the host-referenced model --------------------------------------------------------

def gene_rhythms(logexpr: pd.DataFrame, times_h, genes) -> pd.Series:
    """Complex cosinor coefficient per gene, c = A exp(-i w phase): a tissue whose clock reads later by d hours
    has every coefficient multiplied by exp(i w d)."""
    c = cosinor(logexpr, times_h, genes)
    return c.amplitude * np.exp(-1j * W * c.phase)


def csvd_tissues(C: pd.DataFrame, ref: str) -> pd.DataFrame:
    """Rank-1 complex SVD of the genes x tissues matrix of cosinor coefficients (as Talamanca et al. 2023 for
    tissue rhythms): C ~ s u v^H, so tissue t scales every gene by conj(v_t). Amplitude and phase offset (h,
    positive: the tissue's clock reads later) relative to `ref`, and the share of |C|^2 the first component
    explains. Free of any waveform or clock-plane model: a truth independent of the sccc estimators."""
    C = C.dropna()
    U, s, Vh = np.linalg.svd(C.to_numpy(complex), full_matrices=False)
    v = Vh[0]  # conj(v_t), the tissue factor
    z = v / v[list(C.columns).index(ref)]
    return pd.DataFrame({"amplitude": np.abs(z), "offset_h": np.angle(z) / W, "explained": s[0] ** 2 / (s**2).sum()},
                        index=pd.Index(C.columns, name="compartment"))


def median_tissues(C: pd.DataFrame, ref: str) -> pd.DataFrame:
    """Per-gene version: median log2 amplitude ratio to `ref` and amplitude-weighted circular mean phase difference."""
    C = C.dropna()
    r = C.div(C[ref], axis=0)
    w = np.sqrt(np.abs(C).mul(np.abs(C[ref]), axis=0))
    return pd.DataFrame({"amplitude": 2 ** np.log2(np.abs(r)).median(),
                         "offset_h": np.angle((w * r / np.abs(r)).sum()) / W},
                        index=pd.Index(C.columns, name="compartment"))


def centred(est: pd.DataFrame, tissues=None) -> pd.DataFrame:
    """Reference-free tissue estimands: log2 amplitude minus its mean over `tissues`, and offset (h) minus their
    circular mean. The reference tissue's zero has no intrinsic meaning; differences between tissues do."""
    e = est if tissues is None else est.loc[list(tissues)]
    la = np.log2(e.amplitude)
    off = np.angle(np.exp(1j * W * e.offset_h).mean()) / W
    return pd.DataFrame({"log2_amp": la - la.mean(), "offset_h": wrap(e.offset_h - off)}, index=e.index)


def tissue_errors(est: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """Reference-free estimates and truth (`centred`, over the tissues both have) and their errors: log2 amplitude
    error and offset error (h), per tissue."""
    k = est.dropna(subset=["amplitude", "offset_h"]).index.intersection(truth.index)
    e, t = centred(est, k), centred(truth, k)
    return pd.DataFrame({"log2_amp": e.log2_amp, "log2_amp_true": t.log2_amp, "offset_h": e.offset_h,
                         "offset_h_true": t.offset_h, "log2_amp_err": e.log2_amp - t.log2_amp,
                         "offset_err_h": wrap(e.offset_h - t.offset_h)}, index=k)
