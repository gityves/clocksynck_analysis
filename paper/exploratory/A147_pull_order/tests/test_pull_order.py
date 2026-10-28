import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from pull_order import deviation, pull_order  # noqa: E402


def _C(phases_h, T=8):
    th = np.asarray(phases_h) * 2 * np.pi / 24
    return np.exp(1j * (th[:, None] - th[None, :]))


def test_deviation_poles():
    pole, d = deviation(np.array([1.0, 11.0, -11.0, -2.0]) * 2 * np.pi / 24)
    assert list(pole) == [False, True, True, False]
    assert np.allclose(d, [1, -1, 1, -2])


def test_identity_and_shrink():
    rng = np.random.default_rng(0)
    ph = rng.uniform(0, 24, 8)
    ref = _C(ph)
    r = pull_order(ref, ref, n_perm=200, rng=rng)
    assert np.isclose(r["order"], 1) and np.isclose(r["ratio"], 1) and r["p_order"] < 0.02
    # halving every pair phase deviation: pole kept, order kept, ratio 0.5 (pairs within +-3 h of a pole)
    th = np.angle(ref)
    iu = np.triu_indices(8, 1)
    pole, d = deviation(th[iu])
    half = np.zeros((8, 8), complex)
    # build directly: new phase = pole phase + dev/2
    newh = np.where(pole, np.sign(th[iu]) * 12, 0) + d / 2
    half[iu] = np.exp(1j * newh * 2 * np.pi / 24)
    half = half + half.conj().T
    r = pull_order(ref, half, n_perm=0)
    assert np.isclose(r["ratio"], 0.5) and np.isclose(r["order"], 1)


def test_shuffled_relations_lose_order():
    rng = np.random.default_rng(1)
    ref = _C(rng.uniform(-3, 3, 12))
    grp = _C(rng.uniform(-3, 3, 12))
    r = pull_order(ref, grp, n_perm=500, rng=rng)
    assert abs(r["order"]) < 0.5 and r["p_order"] > 0.01
