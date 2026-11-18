"""A311: the forced Hopf normal form and the three coordinates it predicts (model only, no data).

In the frame of the drive, w = z e^{-i theta(t)} obeys dw/dt = (mu + i nu) w - |w|^2 w + F (+ noise). With F = 1
(units: mu, nu in F^{2/3}, gain in F^{1/3}) the locked state has x = |w|^2 solving x((mu - x)^2 + nu^2) = 1 and
w* = -1 / (mu + i nu - x). Departure: stationary covariance C of additive isotropic noise (J C + C J^T + I = 0, J the
Jacobian at w*), split along w* (radial) and across it (tangential); share = (C_tt - C_rr) / (C_tt + C_rr).

    python hopf.py      -> ../output/hopf_grid.csv (mu, nu, locked, gain, lag_h, var_rad, var_tan, share, elasticity)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_continuous_lyapunov

OUT = Path(__file__).resolve().parents[1] / "output"


def locked_states(mu: float, nu: float) -> list:
    """Stable locked states (w*, J) of the normal form at F = 1."""
    xs = np.roots([1.0, -2 * mu, mu * mu + nu * nu, -1.0])
    out = []
    for x in xs[(abs(xs.imag) < 1e-9) & (xs.real > 0)].real:
        w = -1.0 / (mu + 1j * nu - x)
        a, b = mu + 1j * nu - 2 * x, -w * w  # d(dw/dt) = a dw + b conj(dw)
        J = np.array([[a.real + b.real, -a.imag + b.imag], [a.imag + b.imag, a.real - b.real]])
        if np.max(np.linalg.eigvals(J).real) < 0:
            out.append((w, J))
    return out


def coordinates(mu: float, nu: float = 0.0) -> dict:
    s = locked_states(mu, nu)
    if not s:
        return dict(mu=mu, nu=nu, locked=False)
    w, J = max(s, key=lambda t: abs(t[0]))
    C = solve_continuous_lyapunov(J, -np.eye(2))
    e = np.array([w.real, w.imag]) / abs(w)
    t = np.array([-e[1], e[0]])
    vr, vt = float(e @ C @ e), float(t @ C @ t)
    return dict(mu=mu, nu=nu, locked=True, gain=abs(w), lag_h=float(np.angle(w)) * 12 / np.pi, var_rad=vr,
                var_tan=vt, share=(vt - vr) / (vt + vr))


def elasticity(mu: float, nu: float = 0.0, h: float = 0.01) -> float:
    """d log gain / d log F at F = 1 (1: gain proportional to drive; 0: saturated)."""
    def g(F):
        xs = np.roots([1.0, -2 * mu, mu * mu + nu * nu, -F * F])
        return np.sqrt(xs[(abs(xs.imag) < 1e-9) & (xs.real > 0)].real.max())
    return float((np.log(g(1 + h)) - np.log(g(1 - h))) / (np.log(1 + h) - np.log(1 - h)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    g = pd.DataFrame([coordinates(m, n) for n in (0.0, 0.25, 0.5) for m in np.arange(-4.0, 8.01, 0.5)])
    g["elasticity"] = [elasticity(m, n) if l else np.nan for m, n, l in zip(g.mu, g.nu, g.locked)]
    g.to_csv(OUT / "hopf_grid.csv", index=False)
    z = g[(g.nu == 0) & g.mu.between(-4, 4)]
    slope = np.polyfit(np.log(z.gain ** 2), np.log(z.var_rad + z.var_tan), 1)[0]
    print(g.round(3).to_string())
    print(f"nu 0, mu -4..4: slope of log departure on log gain^2 {slope:.3f}")


if __name__ == "__main__":
    main()
