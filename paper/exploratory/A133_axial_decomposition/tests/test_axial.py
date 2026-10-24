import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from axial import circular_part, fit_axial  # noqa: E402


def _truth(T=8, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.uniform(1, 3, T)
    phi = rng.uniform(-0.5, 0.5, T)
    u = rng.uniform(0.5, 1.5, T) * np.exp(1j * rng.uniform(0, np.pi, T))
    return a, phi, u


def test_noiseless_exact_recovery():
    a, phi, u = _truth()
    v = 2.0
    clock = 3.0 * np.outer(a * np.exp(1j * phi), (a * np.exp(1j * phi)).conj())  # rank-1 Hermitian circular part
    C = clock + v * np.outer(u, u.conj())
    P = v * np.outer(u, u)
    Cc, w = circular_part(C, P)
    assert np.allclose(np.outer(w, w), P, atol=1e-10)
    assert np.allclose(Cc, clock, atol=1e-10)
    iu = np.triu_indices(len(a), 1)
    assert np.allclose(np.angle(Cc[iu]), (phi[:, None] - phi[None])[iu])


def test_missing_pairs_recovered():
    a, phi, u = _truth(10, 1)
    P = 2.0 * np.outer(u, u)
    rng = np.random.default_rng(2)
    obs = rng.uniform(size=P.shape) < 0.7
    obs = np.triu(obs, 1)
    obs = obs | obs.T
    np.fill_diagonal(obs, True)
    w = fit_axial(np.where(obs, P, 0), obs)
    assert np.allclose(np.outer(w, w), P, atol=1e-6)


def test_no_displacement_leaves_c():
    a, phi, _ = _truth()
    z = a * np.exp(1j * phi)
    C = np.outer(z, z.conj())
    Cc, w = circular_part(C, np.zeros_like(C))
    assert np.allclose(Cc, C)
