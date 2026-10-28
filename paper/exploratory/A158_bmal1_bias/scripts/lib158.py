"""A158 library: simulator at the Manella-2021 AL liver-vs-WAT+lung design (local copy of clocksync.simulate with the sampling
times given instead of drawn), clocksync fit and A145 timed oracle (two-step on the true ZT)."""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A145_two_step/scripts"))
import clocksync as cs  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402
import lib as A145  # noqa: E402  (A145 library: two_step, states, regress)

HOURS = 24.0 / (2 * np.pi)
ANCH = ["WAT", "Lung"]


def real_inputs():
    spec = importlib.util.spec_from_file_location("a111c", PAPER / "exploratory/A111_feeding_controls/scripts/_common.py")
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
    from common import data

    pb, donors = data.load("manella2021")
    qc = C.gene_qc(pb, "manella2021")
    wf = C.waveform_without(sorted(qc[~qc.kept].gene.unique()))
    return pb, donors, wf


def sim_table(wf, comps: dict, t_h: np.ndarray, amplitude: float, depth: float, rng, tau=0.05, shift_sd=0.2):
    """Long table (donor, compartment, depth, genes) at given true times (h); same generative model as `cs.simulate`."""
    G, n = len(wf.genes), len(t_h)
    base = np.full(G, -10.0) if wf.mesor is None else np.asarray(wf.mesor, float)
    rows = []
    for name, c in comps.items():
        offset = rng.normal(0.0, c.mesor_sd, G)
        for si in range(n):
            phase = t_h[si] / HOURS if c.coupled else rng.uniform(0, 2 * np.pi)
            zc = amplitude * complex(c.z) * np.exp(1j * phase)
            beta = np.array([zc.real, zc.imag]) + rng.normal(0.0, np.sqrt(c.psi), 2)
            lr = base + offset + wf.X @ beta + rng.normal(0.0, tau, G) + rng.normal(0.0, shift_sd)
            rows.append([f"m{si:03d}", name, depth, *rng.poisson(depth * np.exp(lr))])
    return pd.DataFrame(rows, columns=["donor", "compartment", "depth", *wf.genes])


def pop_v(t_h) -> float:
    """Trace of Cov(cos, sin) of the design's times (V of a unit-amplitude host)."""
    p = np.asarray(t_h) / HOURS
    return float(np.trace(np.cov(np.column_stack([np.cos(p), np.sin(p)]).T)))


def psi_of(c, z2, amp, V1):
    """psi giving population coupling c = a^2 z2 V1 / (a^2 z2 V1 + 2 psi) (clocksync definition) for host amplitude `amp`."""
    return amp * amp * z2 * V1 * (1 - c) / (2 * c)


def fit_cs(tab, wf, anchors, genes=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(tab, anchors=anchors, waveform=wf, name="a158")
        f = ds.fit(n_boot=0, gate=False, targets=["Liver"], force=True)
    t = f.table.loc["Liver"]
    return dict(c=float(t.coupling), z2=float(t.z2), boundary=bool(t.boundary), psi=float(t.psi), V=float(t.V))


def oracle(tab, wf, t_h):
    phi = pd.Series(np.asarray(t_h) / HOURS, index=sorted(tab.donor.unique()))
    r = A145.two_step(tab, phi, targets=["Liver"], anchors=ANCH, genes=list(wf.genes))
    return float(r.coupling["Liver"])
