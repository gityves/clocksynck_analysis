"""A145 library: simulated design, clocksync fit, two-step regression on clocksync's clock-plane states, orderer runners.

Registered in paper/docs/registered_experiments.md (A145). Technical choices not fixed by the registration are listed
in ../README.md ("Choices not in the registration") and were fixed before any method was run on simulated data.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent))
sys.path.insert(0, str(PAPER))

import clocksync as cs  # noqa: E402
from clocksync import paired as PD  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402

OUT = HERE.parent / "output"
HOURS = 24.0 / (2 * np.pi)
WF = cs.gtex()
GENES = list(WF.genes)
N_DONORS, SEED = 100, 145
ANCHORS = ["A1", "A2", "A3", "A4"]
TARGETS = [f"T{i}" for i in range(5)]
C_LEVELS = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
C_ANCHOR = 0.9
OFFSETS_H = np.array([0.0, 3.0, 6.0, 9.0, 12.0])
AMPS = np.array([0.5, 1.0])
ARMS = {  # window, jitter, shared_sd key (None: 0), amplitude
    "F": dict(window=(0.0, 24.0), jitter=0.0, shared=None, amplitude=1.0),
    "W": dict(window=(8.0, 18.0), jitter=1.0, shared=None, amplitude=1.0),
    "S_small": dict(window=(0.0, 24.0), jitter=0.0, shared="small", amplitude=1.0),
    "S_large": dict(window=(0.0, 24.0), jitter=0.0, shared="large", amplitude=1.0),
    "N": dict(window=(0.0, 24.0), jitter=0.0, shared="large", amplitude=0.0),
}
METHODS = ("chiral", "cyclops", "cofe")


def shared_sd(key) -> float:
    """shared_sd of an arm: 0, or the value fixed by the simulator check (output/simulator_check.csv)."""
    if key is None:
        return 0.0
    chk = pd.read_csv(OUT / "simulator_check.csv")
    return float(chk.loc[chk.quantity == f"shared_sd_{key}", "value"].iloc[0])


def host_variance(window, jitter) -> float:
    """V = trace Cov(cos, sin) of the sampled host time (population value, 400,000 draws)."""
    r = np.random.default_rng(0)
    t = r.uniform(*window, 400_000) + (r.normal(0, jitter, 400_000) if jitter > 0 else 0.0)
    return float(np.trace(np.cov(np.column_stack([np.cos(t / HOURS), np.sin(t / HOURS)]).T)))


def psi_of(c: float, z: float, V: float) -> float:
    """psi that gives population coupling c = z^2 V / (z^2 V + 2 psi) (the clocksync definition)."""
    return z * z * V * (1 - c) / (2 * c) if c > 0 else 0.0


def design(rep: int) -> pd.DataFrame:
    """Targets of replicate `rep`: true c in C_LEVELS, offsets in OFFSETS_H and amplitudes in AMPS assigned by independent
    permutations (the same in every arm), so that c, offset and amplitude are not confounded across replicates."""
    r = np.random.default_rng([SEED, rep])
    amp = r.permutation(np.array([0.5, 0.5, 1.0, 1.0, 0.5 if r.random() < 0.5 else 1.0]))
    return pd.DataFrame(dict(c=C_LEVELS, offset_h=r.permutation(OFFSETS_H), amp=amp), index=TARGETS)


def cohort(arm: str, rep: int, n: int = N_DONORS):
    """(long table with `donor`, design table, true times t (h))."""
    a, d = ARMS[arm], design(rep)
    V = host_variance(a["window"], a["jitter"])
    comps = {k: Compartment(z=1.0, psi=psi_of(C_ANCHOR, 1.0, V)) for k in ANCHORS}
    for k, r in d.iterrows():
        z = r.amp * np.exp(1j * r.offset_h / HOURS)
        comps[k] = Compartment(z=z, psi=psi_of(r.c, r.amp, V), coupled=r.c > 0)
    tab, t = cs.simulate(WF, comps, n=n, window_h=a["window"], jitter_h=a["jitter"], amplitude=a["amplitude"],
                         shared_sd=shared_sd(a["shared"]), seed=SEED * 1000 + rep)
    return tab.rename(columns={"resection": "donor"}), d, t


# ---- clocksync ---------------------------------------------------------------------------------------------------

def clocksync_fit(tab: pd.DataFrame, seed: int = 0) -> dict:
    """Dataset.fit, 4 anchors (leave-one-out host), point estimates (no bootstrap), gate computed then bypassed (`force`)
    so that estimates exist when it fails: gate summary and per-target coupling, offset and relative amplitude."""
    ds = cs.Dataset(tab, anchors=ANCHORS, name="a145")
    f = ds.fit(n_boot=0, seed=seed, force=True)
    g = f.gate
    t = f.table
    return dict(gate=dict(passes=bool(g.passes), p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share),
                table=t[["coupling", "offset_h", "log_amp_rel", "z2", "V", "psi"]])


# ---- two-step -------------------------------------------------------------------------------------------------------

def states(tab: pd.DataFrame, comp: str, wf=WF, genes=None):
    """clocksync's clock-plane state of one compartment (own projector from its median counts, library depth, Poisson
    variance): beta (n x 2), trace of its Poisson covariance (n,), donors."""
    x = tab[tab.compartment == comp].reset_index(drop=True)
    g = list(genes or wf.genes)
    w = wf.subset(g)
    L = PD.projector(x[g].astype(float), w, counts=True)
    b, v = PD.states(x[g].astype(float), L, w, depth=x.depth.to_numpy(float))
    return b, v, x.donor.to_numpy()


def regress(beta: np.ndarray, v: np.ndarray, phi: np.ndarray) -> dict:
    """Regress the centred state on (cos, sin) of phi (rad): coupling proxy (explained variance of the fit over the
    total minus mean Poisson variance, clipped to [0, 1]); coefficients a+ (on exp(+i phi)) and a- (on exp(-i phi))."""
    n = len(beta)
    X = np.column_stack([np.ones(n), np.cos(phi), np.sin(phi)])
    B, *_ = np.linalg.lstsq(X, beta, rcond=None)  # 3 x 2
    fit = X[:, 1:] @ B[1:]
    ev = np.sum((fit - fit.mean(0)) ** 2) / (n - 1)
    tot = np.sum((beta - beta.mean(0)) ** 2) / (n - 1) - np.mean(v)
    u, w = B[1, 0] + 1j * B[1, 1], B[2, 0] + 1j * B[2, 1]  # complex coefficients of cos and sin
    return dict(coupling=float(np.clip(ev / tot, 0, 1)) if tot > 0 else np.nan, a_plus=(u - 1j * w) / 2, a_minus=(u + 1j * w) / 2)


def two_step(tab: pd.DataFrame, phi: pd.Series, targets=TARGETS, anchors=ANCHORS, genes=None) -> pd.DataFrame:
    """Per target: coupling proxy, offset (h, relative to the anchors' circular mean) and log relative amplitude (log of
    the modulus over the anchors' mean modulus), from the regression on the donor phase `phi` (rad, arbitrary origin).
    Direction: phi is reflected if the anchors' states then turn with it (sum |a+|^2 - |a-|^2 larger), the clock
    plane's own rotation sense (clocksync: states turn counter-clockwise as the day advances)."""
    st = {k: states(tab, k, genes=genes) for k in list(anchors) + list(targets)}

    def fit_all(sign):
        out = {}
        for k, (b, v, d) in st.items():
            p = phi.reindex(d).to_numpy(float)
            ok = np.isfinite(p)
            out[k] = regress(b[ok], v[ok], sign * p[ok])
        return out
    best = None
    for sign in (1.0, -1.0):
        r = fit_all(sign)
        score = sum(abs(r[k]["a_plus"]) ** 2 - abs(r[k]["a_minus"]) ** 2 for k in anchors)
        if best is None or score > best[0]:
            best = (score, r)
    r = best[1]
    anc_arg = np.angle(np.sum([r[k]["a_plus"] / abs(r[k]["a_plus"]) for k in anchors]))
    anc_mod = np.mean([abs(r[k]["a_plus"]) for k in anchors])
    return pd.DataFrame({k: dict(coupling=r[k]["coupling"],
                                 offset_h=float(S.wrap((np.angle(r[k]["a_plus"]) - anc_arg) * HOURS)),
                                 log_amp_rel=float(np.log(abs(r[k]["a_plus"]) / anc_mod))) for k in targets}).T


def clocksync_relative(t: pd.DataFrame, targets=TARGETS, anchors=ANCHORS) -> pd.DataFrame:
    """clocksync's offsets and log amplitudes relative to the anchors (circular mean offset, mean log amplitude)."""
    ref = S.circmean(t.loc[anchors, "offset_h"].to_numpy())
    out = pd.DataFrame(index=targets)
    out["coupling"] = t.loc[targets, "coupling"]
    out["offset_h"] = S.wrap(t.loc[targets, "offset_h"].to_numpy() - ref)
    out["log_amp_rel"] = t.loc[targets, "log_amp_rel"] - t.loc[anchors, "log_amp_rel"].mean()
    return out


# ---- orderers ---------------------------------------------------------------------------------------------------------

def orderer_input(tab: pd.DataFrame, anchors=ANCHORS, genes=None) -> pd.DataFrame:
    """Donors x (gene, anchor) matrix of log2 normalised counts, anchors concatenated (multi-tissue CHIRAL input)."""
    from common.benchmark.data import log2_normalised

    g = list(genes or GENES)
    cols = []
    for k in anchors:
        x = tab[tab.compartment == k].set_index("donor")[g].astype(float)
        L = log2_normalised(x)
        L.columns = [f"{c}_{k}" for c in L.columns]
        cols.append(L)
    return pd.concat(cols, axis=1).dropna()


def run_orderer(method: str, L: pd.DataFrame, tmp: Path, threads: int = 1, seed: int = 0) -> tuple[pd.Series | None, str, float]:
    """(phase in radians per donor or None on failure, error text, seconds). Clock-gene seed: every feature is a clock gene."""
    t0 = time.time()
    try:
        if method == "chiral":
            from common.benchmark.chiral import run_chiral

            ph, failed = run_chiral(L, genes=tuple(L.columns), seed=seed)
            ph = ph / 24 * 2 * np.pi
            err = f"failed_starts={failed}" if ph.isna().all() else ""
        elif method == "cyclops":
            from common.benchmark.cyclops import run_cyclops

            r = run_cyclops((2 ** L - 1).T, list(L.columns), tmp, seed=seed, threads=threads)
            ph, err = r["phase_raw"].astype(float).set_axis(L.index), ""
        elif method == "cofe":
            from common.benchmark.cofe import run_cofe

            r = run_cofe(L.loc[:, L.std() > 0], seed=seed, ncores=threads)
            ph, err = (r["phase"] / 24 * 2 * np.pi).astype(float).set_axis(L.index), ""
        else:
            raise ValueError(method)
        if ph.isna().all():
            return None, err or "all phases NaN", time.time() - t0
        return ph, err, time.time() - t0
    except Exception as e:  # a failed fit is data
        return None, f"{type(e).__name__}: {e}"[:300], time.time() - t0
