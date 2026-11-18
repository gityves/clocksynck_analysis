"""A310: mean-field population model of a tissue clock (Ott-Antonsen reduction of driven, coupled phase oscillators).

Cells of a tissue are phase oscillators (fixed single-cell amplitude a_c) with Lorentzian intrinsic frequencies (centre
omega0, half-width gamma), mean-field coupling K and a systemic drive F e^{i Omega t}. In the frame of the drive, the
order parameter w = R e^{i phi} obeys (Ott and Antonsen 2008)

    dw/dt = (i nu - gamma) w + (H - conj(H) w^2) / 2,   H = K w + F,   nu = omega0 - Omega.

A person's tissue departs from the locked state w* by additive white noise of intensity q; its stationary covariance
around w* is the solution C of J C + C J^T + q I = 0 (J the Jacobian at w*), split into the radial (along w*) and
tangential (across w*) components. The linear alternative (damped, driven, no synchronisation) is
dw/dt = (i nu - lam) w + F, whose covariance is isotropic.

Model-only functions: no data are read here.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov


def rhs(w: complex, gamma: float, K: float, F: complex, nu: float = 0.0) -> complex:
    H = K * w + F
    return (1j * nu - gamma) * w + 0.5 * (H - np.conj(H) * w * w)


def locked(gamma: float, K: float, F: complex, nu: float = 0.0, w0: complex = 0.5, t_max: float = 400.0) -> complex:
    """Stable locked state w* (integrated from w0; NaN if not settled)."""
    def f(_, y):
        d = rhs(y[0] + 1j * y[1], gamma, K, F, nu)
        return [d.real, d.imag]
    s = solve_ivp(f, [0, t_max], [w0.real, w0.imag], rtol=1e-10, atol=1e-12)
    w = s.y[0, -1] + 1j * s.y[1, -1]
    return w if abs(rhs(w, gamma, K, F, nu)) < 1e-6 else complex(np.nan, np.nan)


def jacobian(w: complex, gamma: float, K: float, F: complex, nu: float = 0.0, h: float = 1e-6) -> np.ndarray:
    J = np.zeros((2, 2))
    for k, d in enumerate((h, 1j * h)):
        df = (rhs(w + d, gamma, K, F, nu) - rhs(w - d, gamma, K, F, nu)) / (2 * h)
        J[:, k] = df.real, df.imag
    return J


def departure(gamma: float, K: float, F: complex, nu: float = 0.0, q: float = 1.0) -> dict:
    """R = |w*|, phase of w*, and the radial / tangential variance of the stationary departure."""
    w = locked(gamma, K, F, nu)
    C = solve_continuous_lyapunov(jacobian(w, gamma, K, F, nu), -q * np.eye(2))
    e = np.array([w.real, w.imag]) / abs(w)
    t = np.array([-e[1], e[0]])
    vr, vt = float(e @ C @ e), float(t @ C @ t)
    return dict(R=abs(w), phase=float(np.angle(w)), var_rad=vr, var_tan=vt, share=(vt - vr) / (vt + vr))
