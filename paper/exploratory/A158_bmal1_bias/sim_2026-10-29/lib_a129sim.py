"""Mechanism simulator for the A129 / D145.2 gap (clocksync -0.293 vs timed oracle -0.186), roadmap 4.6, DESCRIPTIVE.

Extends `A158_bmal1_bias/scripts/lib158` (Manella-2021 AL design: real ZTs, genes after QC, depth 3.2e6, 24 mice per genotype,
Liver read on WAT + Lung) with two ingredients the A158 simulator lacked:
  * a per-mouse phase deviation `eps_m ~ N(0, jitter_sd)` shared by every tissue of the mouse, weighted per tissue (`w`):
    the host state follows it, the true ZT does not;
  * an independent per-mouse phase noise of a tissue (`phase_sd`) and a per-tissue amplitude multiplier.
Estimators are the public `clocksync.Dataset.fit` (host = WAT + Lung, leave-one-out for anchors) and the A145 timed oracle.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "scripts"))
from lib158 import ANCH, HOURS, Compartment, fit_cs, oracle, pop_v, psi_of, real_inputs  # noqa: E402,F401


def sim_table(wf, comps, t_h, amplitude, depth, rng, jitter_sd=0.0, w=None, phase_sd=None, amp_mult=None, tau=0.05, shift_sd=0.2):
    """As `lib158.sim_table`, plus shared mouse jitter (weight `w[name]`, default 1), independent phase noise
    `phase_sd[name]` (rad) and amplitude multiplier `amp_mult[name]`. All coupled tissues follow ZT + w * eps_m."""
    w, phase_sd, amp_mult = w or {}, phase_sd or {}, amp_mult or {}
    G, n = len(wf.genes), len(t_h)
    base = np.full(G, -10.0) if wf.mesor is None else np.asarray(wf.mesor, float)
    eps = rng.normal(0.0, jitter_sd, n) if jitter_sd > 0 else np.zeros(n)
    rows = []
    for name, c in comps.items():
        offset = rng.normal(0.0, c.mesor_sd, G)
        for si in range(n):
            ph = t_h[si] / HOURS + w.get(name, 1.0) * eps[si] + rng.normal(0.0, phase_sd.get(name, 0.0))
            zc = amplitude * amp_mult.get(name, 1.0) * complex(c.z) * np.exp(1j * ph)
            beta = np.array([zc.real, zc.imag]) + rng.normal(0.0, np.sqrt(c.psi), 2)
            lr = base + offset + wf.X @ beta + rng.normal(0.0, tau, G) + rng.normal(0.0, shift_sd)
            rows.append([f"m{si:03d}", name, depth, *rng.poisson(depth * np.exp(lr))])
    return pd.DataFrame(rows, columns=["donor", "compartment", "depth", *wf.genes])
