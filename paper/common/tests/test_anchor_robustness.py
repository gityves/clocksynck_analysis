"""`common.anchor_robustness` on a synthetic cohort with known structure (about 1 min)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "clocksync"))
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import anchor_robustness as ar  # noqa: E402

limit_threads(2)
WF = cs.gtex()
N = 60


@pytest.fixture(scope="module")
def ds():
    z = lambda a, h: a * np.exp(1j * h / cs.waveform.HOURS)
    comps = {k: cs.Compartment(z(1.0, h), psi=0.002, mesor_sd=0.0) for k, h in
             (("fibroblast", 0.0), ("endothelial", 1.0), ("pericyte", -1.0), ("B cell", 0.5), ("T cell", 0.0),
                      ("NK", 0.3), ("mono", -0.4), ("epi", 0.8))}
    df, t = cs.simulate(WF, comps, n=N, seed=4, window_h=(4, 22), jitter_h=1.0, tau=0.0)
    df = df.rename(columns={"resection": "donor"})
    donors = pd.DataFrame({"lab": np.where(np.arange(N) % 2, "case", "control")}, index=[f"r{i:03d}" for i in range(N)])
    return cs.Dataset(df, anchors=["fibroblast", "endothelial", "pericyte", "B cell"], donors=donors)


def test_block_agreement_one_host(ds):
    r = ar.block_agreement(ds, n_boot=0, n_split=10)
    assert r.n_anchors == 4 and r.n_rest == 4 and r.r_star > 0.8


def test_block_agreement_small_blocks(ds):
    r = ar.block_agreement(ds, anchors=["fibroblast", "endothelial", "pericyte"], n_boot=0, n_split=5)  # rest of 5
    assert np.isfinite(r.r) and r.r > 0.8
    small = ar.block_agreement(ds, anchors=ds.compartments[:7])  # one compartment left: undefined
    assert np.isnan(small.r_star) and "block" in small.note
    two = ar.block_agreement(ds, anchors=ds.compartments[:6], n_boot=0, n_split=5)  # rest of 2: r only
    assert np.isnan(two.r_star) and np.isfinite(two.r) and "r* undefined" in two.note


def test_loo_dataset_all_anchors(ds):
    assert ar.loo_dataset(ds).anchors == ds.compartments and ds.anchors != ds.compartments


def test_sign_change():
    s = ar.sign_change([1.0, -1.0, 0.5, np.nan], [-1.0, -2.0, 0.3, 1.0])
    assert s.tolist()[:3] == [True, False, False] and pd.isna(s.iloc[3])


def test_host_contrast_columns_and_agreement(ds):
    t = ar.host_contrast(ds, "lab", ["T cell", "fibroblast"], n_perm=10, n_boot=0)
    assert set(t.target) == {"T cell", "fibroblast"} and {"d_anchored", "d_loo", "sign_change", "p_star_loo"} <= set(t)
    # no case-control effect in the simulation and one shared host: both hosts give the same d within noise
    assert (t.d_anchored - t.d_loo).abs().max() < 0.2
    # the anchored host of an anchor target excludes it, as does the LOO host: identical host for fibroblast only if
    # the non-anchors add nothing; at least the numbers must be finite
    assert np.isfinite(t[["d_anchored", "d_loo"]].to_numpy()).all()
