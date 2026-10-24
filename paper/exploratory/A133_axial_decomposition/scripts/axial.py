"""Circular / axial decomposition of a host covariance (A133).

Model. Centred clock-plane states z_a (complex, one per tissue a) of a donor are
    z_a = A_a e^{i(theta + phi_a)} + s u_a + noise,
theta the host phase, s a real donor-level displacement (variance sigma^2), u_a a complex direction.
For uniformly covered theta and circular noise, with C = E[z conj z^T] and P = E[z z^T]:
    C = C_circ + sigma^2 u u^H,   P = sigma^2 u u^T  (+ 0 from the rotating part, + noise on the diagonal of C only).
So P is rank 1 and complex symmetric, with Takagi factorisation P = sigma^2 u u^T; the axial part of C is the same
rank-1 matrix w w^H with w = sigma u, and C_circ = C - w w^H holds the running clock: arg C_circ[a, b] = phi_a - phi_b.

Takagi factor (rank 1). For symmetric P = U S U^T (U unitary), the leading singular triplet of the SVD P = W S V^H has
V = conj(W) up to a phase e^{-i alpha}; w = sqrt(s_1) W[:, 0] e^{i alpha / 2} satisfies P ~ w w^T (w defined up to sign,
which cancels in w w^H).

Missing pairs. Pairs with too few shared donors are unobserved in P. We use iterative imputation: unobserved entries
are filled with the current rank-1 fit (start: zero), the Takagi factor is refit on the completed matrix, repeated to
convergence. The diagonal of P (always observed) is noise-free in expectation for circular noise and is used in the fit.
Fitted entries equal the data at the observed entries only up to the rank-1 constraint; exact for a noiseless rank-1 P.
"""

from __future__ import annotations

import numpy as np


def takagi_rank1(P: np.ndarray) -> np.ndarray:
    """w with P ~ w w^T for a complex symmetric P (leading Takagi pair, sigma included in w)."""
    W, S, Vh = np.linalg.svd(P)
    # P = s W0 Vh0 (Vh0 = conj(V0)^T); symmetry gives Vh0 = e^{i alpha} W0^T, so alpha = angle(sum Vh0 conj(W0))
    alpha = np.angle(np.sum(Vh[0] * W[:, 0].conj()))
    return np.sqrt(S[0]) * W[:, 0] * np.exp(1j * alpha / 2)


def fit_axial(P: np.ndarray, observed: np.ndarray | None = None, n_iter: int = 500, tol: float = 1e-10) -> np.ndarray:
    """Rank-1 axial factor w (P ~ w w^T) from the pseudo-covariance P over the entries flagged in `observed`
    (symmetric boolean matrix; default: all). Iterative imputation, see the module docstring."""
    P = np.asarray(P, complex)
    P = (P + P.T) / 2
    if observed is None or observed.all():
        return takagi_rank1(P)
    obs = observed | observed.T
    np.fill_diagonal(obs, True)
    X = np.where(obs, P, 0)
    w = takagi_rank1(X)
    for _ in range(n_iter):
        X = np.where(obs, P, np.outer(w, w))
        w_new = takagi_rank1(X)
        if np.linalg.norm(np.outer(w_new, w_new) - np.outer(w, w)) <= tol * max(np.linalg.norm(np.outer(w, w)), 1e-300):
            w = w_new
            break
        w = w_new
    return w


def circular_part(C: np.ndarray, P: np.ndarray, observed: np.ndarray | None = None):
    """(C_circ, w): C_circ = C - w w^H with w the rank-1 Takagi factor of P (observed: pairs with enough donors)."""
    w = fit_axial(P, observed)
    return C - np.outer(w, w.conj()), w
