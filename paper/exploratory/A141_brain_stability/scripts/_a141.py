"""Shared helpers of A141 (brain contrast stability): contrast sets, fits, host statistics, state-level simulator.

The brain contrast sets are those of WP0 item 8 (`paper/exploratory/WP0_integrity/scripts/anchor_robustness.py`,
`forest_one`): the A17 registry rows with a brain tissue, A17 learned anchors, targets = the rows of
`anchor_robustness.csv` of that contrast (figure 4, `brain:` targets). Cohort = dataset; the control set of a cohort is
that of the contrast with the most controls (A135), see `COHORTS`.
"""

from __future__ import annotations

import copy
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
PAPER = ROOT / "paper"
sys.path.insert(0, str(ROOT / "clocksync" / "src"))
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from common import data  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
OBS = PAPER / "exploratory/WP0_integrity/output/anchor_robustness/anchor_robustness.csv"
TISSUES = PAPER / "results/figure_2/output/organs/tissues.csv"
GROUPS = ("control", "case")
# cohort (dataset) -> contrasts served, and the contrast whose controls (the largest set) define the cohort's controls
COHORTS = {
    "psychad_mssm": (["AD_mssm", "DLBD_mssm", "PD_mssm", "SCZ_mssm", "VaD_mssm"], "SCZ_mssm"),
    "psychad_hbcc": (["ASCVD_hbcc", "BD_hbcc", "SCZ_hbcc"], "ASCVD_hbcc"),
    "psychad_radc": (["AD_radc", "DM_radc"], "DM_radc"),
    "seaad_dlpfc": (["AD_seaad_dlpfc"], "AD_seaad_dlpfc"),
    "parkinson": (["pd"], "pd"),
}
CONTRASTS = [c for v in COHORTS.values() for c in v[0]]


def _mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def observed() -> pd.DataFrame:
    """Rows of the WP0 item 8 table for the brain contrasts of Figure 4 (observed anchored and all-LOO contrasts)."""
    d = pd.read_csv(OBS)
    d = d[d.target.str.startswith("brain:") & (d.figure.astype(str) == "4") & d.contrast.isin(CONTRASTS)]
    return d.reset_index(drop=True)


class Contrast:
    """One brain contrast at its real design: `ds` (A17 anchors), labels, targets (those of the WP0 table)."""

    def __init__(self, name: str, subsample: int | None = None, seed: int = 0):
        res = _mod("res_shared", PAPER / "results/figure_4/scripts/resilience/_shared.py")
        row = pd.Series(res.targets("contrast").set_index("name", drop=False).loc[name], name=name)
        pb, donors, lab, blocks = res.load_contrast(row, subsample, seed)
        self.name, self.lab, self.blocks = name, lab, blocks
        anc = [a for a in res.anchor_set(name, "learned") if a in set(pb.compartment)]
        self.ds = cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab), name=name)
        counts = pb.assign(g=pb.donor.map(lab)).groupby(["compartment", "g"]).donor.nunique().unstack(fill_value=0)
        obs = observed()
        obs = obs[obs.contrast == name]
        self.targets = [k for k in obs.target if k in set(self.ds.compartments)]
        if subsample:  # keep the registered target rule on the subsample
            fam = res.family_map()
            self.targets = [k for k in self.ds.compartments if fam.get(k, "other") != "other"
                            and counts["control"].get(k, 0) >= 5 and counts["case"].get(k, 0) >= 5]
        self.obs = obs.set_index("target")
        self.anchors = anc

    @property
    def n_groups(self):
        return int((self.lab == "control").sum()), int((self.lab == "case").sum())

    def loo(self, ds=None) -> cs.Dataset:
        return with_data(self.ds if ds is None else ds, None, self.ds.compartments)


def with_data(ds: cs.Dataset, d: A.Data | None, anchors=None) -> cs.Dataset:
    """Shallow copy of `ds` holding clock-plane data `d` (simulated states) and/or other anchors."""
    s = copy.copy(ds)
    if d is not None:
        s.data = d
    if anchors is not None:
        s.anchors = list(anchors)
    return s


def d_values(ds: cs.Dataset, targets, blocks) -> pd.Series:
    """Case-minus-control coupling difference per target (no permutations, no bootstrap)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = ds.contrast("lab", targets=targets, levels=GROUPS, blocks=blocks, n_perm=0, n_boot=0).table
    return t.d.reindex(targets).astype(float)


def n_flips(d_a: pd.Series, d_b: pd.Series) -> int:
    """Targets whose d changes sign (non-zero opposite signs, both finite), as `common.anchor_robustness.sign_change`."""
    s = np.sign(d_a.to_numpy()) * np.sign(d_b.to_numpy())
    return int(np.nansum(s < 0))


# ---- host statistics ----------------------------------------------------------------------------------------------

def host_stats(d: A.Data, anchors, donors=None) -> dict:
    """Host fitted on `anchors` (donors of `donors` only): n, Rayleigh R of its phases, spread V = tr(cov f + mean Su)
    (host means centred within stratum, as `on_host`), Rayleigh p (exp(-n R^2)) and the eccentricity of the host
    covariance (smallest / largest eigenvalue: 1 = round, a full cycle). V depends on the frame of the host (the
    factor model fixes the scale only through its reference compartment); R and eccentricity do not."""
    if donors is not None:
        d = d.take(np.flatnonzero(np.isin(d.resection, donors)))
    h = A.host(d, [a for a in anchors if a in set(d.compartment)])
    f = h.f
    st = pd.Series(d.stratum).groupby(pd.Series(d.resection)).first().reindex(h.resections).to_numpy()
    R = float(np.abs(np.mean(np.exp(1j * np.angle(f[:, 0] + 1j * f[:, 1])))))
    S = np.cov(A._centre_by(f, st).T) + h.f_cov.mean(0)
    ev = np.linalg.eigvalsh(S)
    V = float(np.trace(S))
    return dict(n=len(f), R=R, V=V, ecc=float(ev[0] / ev[1]), p_rayleigh=float(np.exp(-len(f) * R ** 2)), converged=bool(h.converged))


# ---- fitted parameters and the simulator ----------------------------------------------------------------------------

def fitted_params(ds: cs.Dataset) -> pd.DataFrame:
    """Per compartment on the contrast's A17 host: absolute coupled amplitude a = sqrt(z2 V) (clock-plane units),
    phase offset (h, common frame: reference-free centring of `Dataset.fit`), own noise psi."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = ds.fit(n_boot=0, gate=False).table
    p = pd.DataFrame(dict(amp=np.sqrt(t.z2 * t.V), offset_h=t.offset_h, psi=t.psi, z2=t.z2, V=t.V, log_amp_rel=t.log_amp_rel,
                          coupling=t.coupling, n=t.n, anchor=t.anchor))
    return p


class Simulator:
    """Clock-plane states at the real design (rows, V, strata, intercept groups, donors) with known truth.

    beta_{s,k} = a_k exp(i theta_k) exp(i phi_s) + N(0, V_row + psi_k I): host phase phi_s uniform on the cycle (24 h
    sampling), amplitude and phase offset of each compartment as fitted, own noise psi_k as fitted, no case effect. Same
    model as `clocksync.simulate.simulate` at the level of the clock-plane state, with the measurement covariance V of
    the real rows instead of a Poisson draw at one depth."""

    def __init__(self, ds: cs.Dataset, par: pd.DataFrame):
        self.ds, self.par = ds, par
        d = ds.data
        self.d = d
        donors = np.unique(d.resection)
        self.pos = pd.Series(np.arange(len(donors)), index=donors)[d.resection].to_numpy()
        self.ndon = len(donors)
        comp = pd.Series(d.compartment)
        p = par.reindex(comp)
        fill = dict(amp=0.0, psi=float(par.psi.median()), offset_h=0.0)
        self.a = p.amp.fillna(fill["amp"]).to_numpy()
        self.th = p.offset_h.fillna(0.0).to_numpy() * 2 * np.pi / HOURS
        self.psi = np.maximum(p.psi.fillna(fill["psi"]).to_numpy(), 1e-8)
        self.L = np.linalg.cholesky(d.V + self.psi[:, None, None] * np.eye(2) + 1e-12 * np.eye(2))

    def draw(self, rng: np.random.Generator) -> A.Data:
        phi = rng.uniform(0, 2 * np.pi, self.ndon)[self.pos]
        z = self.a * np.exp(1j * (self.th + phi))
        b = np.column_stack([z.real, z.imag]) + np.einsum("nab,nb->na", self.L, rng.standard_normal((len(z), 2)))
        return A.Data(b, self.d.V, self.d.resection, self.d.compartment, self.d.stratum, self.d.group)
