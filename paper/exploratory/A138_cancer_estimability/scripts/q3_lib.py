"""A138 Q1/Q3 statistics: states, permutation nulls and simulation designs of H138.1, H138.4 and H138.5.

Registered in paper/docs/registered_experiments.md (A138, note 3 of 2026-10-27). Importable by the power step
(`power.py`) and, unchanged, by the later scoring script. This module computes no observed statistic by itself: the
functions return states and nulls; scoring applies `corr_stat` to the unpermuted pairing.

H138.4  per cancer: tumour and normal states on one projector (`clocksync.paired.cohort_states`), centred within sample
        type; statistic = Poisson-corrected |complex correlation| (`clocksync.paired.cross`, `corr`); pooled = mean over
        gated cancers weighted by n_pairs; null = normals' rows permuted across patients within cancer.
H138.5  per single-cell cohort: malignant (copy-number-called, else tumour epithelium) and adjacent-normal epithelium
        states, each replaced by its residual after complex least-squares regression on the host state of the patient
        (`Dataset.host` on the stromal anchors), Poisson-corrected as `cross`; pooled weights n; null = normal residuals
        permuted within cohort. Companion: tumour vs adjacent fibroblasts, host = the other stromal compartments.
H138.1  simulated cohorts at the design of each Q1 cohort (`clocksync.simulate.simulate`), anchored coupling of the
        malignant compartment with donor bootstrap, DerSimonian-Laird over cohorts.

Choices not in the registration are marked `CHOICE` below and listed in the report of the power step.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "paper"))
sys.path.insert(0, str(REPO / "paper" / "results" / "figure_5" / "scripts" / "cancer"))
import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync import simulate as sim  # noqa: E402
from clocksync.paired import cohort_states, projector, reference_gate, states  # noqa: E402
from clocksync.paired import cross as paired_cross  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common import data  # noqa: E402
from common.anchors import TABLE, anchors  # noqa: E402

import _shared as sh  # noqa: E402  (load_tcga, pair_frames, WF, STROMA, MAL, NORM)

WF = sh.WF
GENES = list(WF.genes)
STROMA = sh.STROMA
MAL, NORM, EPI = sh.MAL, sh.NORM, "tumour:epithelial"
A17_TABLE = REPO / "paper" / "results" / "figure_5" / "output" / "cancer" / "a17_cancer.csv"
N_PERM = 10_000
MIN_PAIRED = 6


# ---- statistic ---------------------------------------------------------------------------------------------------

def centred_complex(b: np.ndarray) -> np.ndarray:
    z = S.complex_state(b)
    return z - z.mean()


def tr_corrected(z: np.ndarray, v: np.ndarray) -> float:
    """Poisson-corrected variance (trace) of centred complex states: sum |z|^2 / (n - 1) - mean(v) (as `paired._tr`)."""
    return float(np.sum(np.abs(z) ** 2) / (len(z) - 1) - np.mean(v))


def corr_stat(za: np.ndarray, vA: np.ndarray, zb: np.ndarray, vB: np.ndarray) -> float:
    """|complex correlation| of centred complex states, Poisson-corrected (`clocksync.paired.cross`, `corr`); nan if
    a corrected variance is not positive. The row order of zb is the pairing."""
    ea, eb = tr_corrected(za, vA), tr_corrected(zb, vB)
    if ea <= 0 or eb <= 0:
        return np.nan
    return float(abs(np.sum(za * np.conj(zb)) / (len(za) - 1)) / np.sqrt(ea * eb))


def perm_null(za, vA, zb, vB, n_perm: int, rng) -> np.ndarray:
    """Statistic with the rows of B permuted across units (pairing broken, both marginals kept): `n_perm` uniform
    random permutations. The corrected variances do not depend on the permutation."""
    n = len(za)
    ea, eb = tr_corrected(za, vA), tr_corrected(zb, vB)
    if ea <= 0 or eb <= 0:
        return np.full(n_perm, np.nan)
    perms = np.argsort(rng.random((n_perm, n)), axis=1)
    return np.abs((za[None, :] * np.conj(zb[perms])).sum(1) / (n - 1)) / np.sqrt(ea * eb)


def estimable(za, vA, zb, vB, null) -> bool:
    """A138 note 4 (ii): a cohort is estimable if the Poisson-corrected variance is > 0 in both compartments and the
    mean of its permutation null is <= 1 (the correction has not broken down). Same call at scoring."""
    ok = tr_corrected(za, vA) > 0 and tr_corrected(zb, vB) > 0 and np.isfinite(null).any()
    return bool(ok and np.nanmean(null) <= 1.0)


def pool(stats: np.ndarray, weights) -> np.ndarray:
    """Weighted mean over cohorts (rows: cohorts, columns: permutations); a cohort with an undefined statistic is
    left out of that column (weights renormalised)."""
    s, w = np.asarray(stats, float), np.asarray(weights, float)[:, None] * np.ones_like(stats)
    w = np.where(np.isfinite(s), w, 0.0)
    return np.nansum(s * w, axis=0) / np.where(w.sum(0) > 0, w.sum(0), np.nan)


def summarise(null: np.ndarray) -> dict:
    """Null mean, SD, 95th percentile, smallest attainable p and MDE (q95 + 0.84 sd; one-sided 0.05, power 0.8)."""
    x = null[np.isfinite(null)]
    if len(x) < 2:
        return dict(null_mean=np.nan, null_sd=np.nan, null_q95=np.nan, n_perm=0, min_p=np.nan, mde=np.nan)
    m, sd, q = float(x.mean()), float(x.std(ddof=1)), float(np.quantile(x, 0.95))
    return dict(null_mean=m, null_sd=sd, null_q95=q, n_perm=len(x), min_p=1 / (len(x) + 1), mde=q + 0.84 * sd)


# ---- H138.4: TCGA paired -----------------------------------------------------------------------------------------

def tcga_states(subsample: int | None = None, seed: int = 0, n_gate: int = 200) -> dict:
    """Per cancer: centred complex tumour / normal states (`za`, `zb`), Poisson variances, n_pairs and the gate of the
    normals (`reference_gate`, A52 rule: p < 0.05). One projector for both sample types."""
    samples = sh.load_tcga(subsample, seed)
    out = {}
    for i, cancer in enumerate(sh.CANCERS):
        case, ref = sh.pair_frames(samples, cancer)
        if len(case) < 10:
            continue
        yr = np.log(np.maximum(ref[GENES].astype(float), 1e-300)).sub(np.log(ref.depth.to_numpy(float)), axis=0).to_numpy()
        g = reference_gate(yr, WF, n_gate, np.random.default_rng([52, seed, i]))
        bt, vt, bn, vn = cohort_states(case, ref, True, WF)
        out[cancer] = dict(za=centred_complex(bt), vA=vt, zb=centred_complex(bn), vB=vn, n=len(case), gate_p=g["p"],
                           gated=bool(g["p"] < 0.05), raw=(bt, vt, bn, vn))
    return out


# ---- H138.5: single cell ------------------------------------------------------------------------------------------

def study_table(study: str) -> tuple[pd.DataFrame, str]:
    """Main table of a cohort (cnv=False) and the name of its malignant compartment. CHOICE: where the cohort has
    copy-number calls for >= 5 patients, the tumour epithelium is replaced by the called malignant compartment (cnv=True
    table, malignant rows only); otherwise tumour epithelium."""
    pb, _ = data.load("cancer_sc", study=study, cnv=False)
    if study in sh.sc_builder().CNV_STUDIES:
        pc, _ = data.load("cancer_sc", study=study, cnv=True)
        m = pc[pc.compartment == MAL]
        if m.donor.nunique() >= 5:
            return pd.concat([pb[pb.compartment != EPI], m], ignore_index=True), MAL
    return pb, EPI


def cls(c: pd.Series) -> pd.Series:
    return c.str.split(":").str[1]


def a17_gate() -> dict:
    """Cohort -> A17 gate passed (gate columns only of the A17 table, learned rule)."""
    a = pd.read_csv(A17_TABLE, usecols=["target", "rule", "gate_passed"]).drop_duplicates()
    a = a[a.rule == "learned"]
    return {t.removeprefix("cancer_"): bool(g) for t, g in zip(a.target, a.gate_passed)}


def stromal_anchors(study: str, present, exclude_class=()) -> tuple[list, str]:
    """Stromal host compartments of a cohort. CHOICE: the learned A17 anchors of the cohort in the stromal classes;
    if fewer than 2 (or none registered), every stromal compartment present. `exclude_class`: classes left out
    (fibroblast for the companion). Returns (compartments, source)."""
    present = [c for c in present if cls(pd.Series([c])).iloc[0] in STROMA and cls(pd.Series([c])).iloc[0] not in exclude_class]
    try:
        a = [c for c in anchors(f"cancer_{study}", "learned") if c in present]
    except KeyError:
        a = []
    return (a, "a17_learned") if len(a) >= 2 else (present, "all_stromal")


def host_of(pb: pd.DataFrame, comps: list) -> pd.DataFrame | None:
    """Host state per patient (`Dataset.host()`: x, y, posterior covariance) on `comps`; None if not estimable."""
    t = pb[pb.compartment.isin(comps)]
    if len(comps) < 1 or t.donor.nunique() < 5:
        return None
    try:
        return cs.Dataset(t, anchors=list(comps)).host()
    except Exception:  # noqa: BLE001  (too few donors / compartments, host not converged)
        return None


def residual(m: np.ndarray, vm: np.ndarray, h: np.ndarray, vh: np.ndarray):
    """Host-referenced state: residual of the centred complex state m after its complex least-squares regression on the
    centred host state h, r = m - a h, a = sum conj(h) m / sum |h|^2. CHOICE (Poisson correction): the measurement
    variance (trace) of r is v_m + |a|^2 v_h (host posterior trace); the residual is scaled by sqrt((n-1)/(n-2)) for
    the fitted complex coefficient (one complex parameter on n-1 centred units, as `Dataset`'s psi correction), so
    that sum |r|^2 / (n - 1) - mean(v) is unbiased for the residual clock-plane variance."""
    n = len(m)
    a = np.sum(np.conj(h) * m) / np.sum(np.abs(h) ** 2)
    r = m - a * h
    r = (r - r.mean()) * np.sqrt((n - 1) / (n - 2))
    return r, vm + np.abs(a) ** 2 * vh


def sc_pair(pb: pd.DataFrame, host: pd.DataFrame, a: str, b: str, strata=True) -> dict | None:
    """Host-referenced residual states of compartments a and b (one row per patient) over the patients having both and
    a host. States on one projector (`clocksync.paired.projector`), centred within stratum; None if fewer than
    MIN_PAIRED patients."""
    A_, B_ = (pb[pb.compartment == c].drop_duplicates("donor").set_index("donor") for c in (a, b))
    pts = sorted(set(A_.index) & set(B_.index) & set(host.index))
    if len(pts) < MIN_PAIRED:
        return None
    A_, B_ = A_.loc[pts], B_.loc[pts]
    L = projector(pd.concat([A_, B_]), WF, counts=True)
    stratum = A_.stratum.astype(str).to_numpy() if strata and "stratum" in A_ else np.zeros(len(pts), int)

    def z(x):
        beta, v = states(x, L, WF, depth=x.depth.to_numpy(float))
        zz = S.complex_state(beta)
        return zz - pd.Series(zz).groupby(stratum).transform("mean").to_numpy(), v

    (za, va), (zb, vb) = z(A_), z(B_)
    H = host.loc[pts]
    h = (H.x + 1j * H.y).to_numpy()
    h = h - pd.Series(h).groupby(stratum).transform("mean").to_numpy()
    vh = (H.xx + H.yy).to_numpy()
    ra, vra = residual(za, va, h, vh)
    rb, vrb = residual(zb, vb, h, vh)
    return dict(za=ra, vA=vra, zb=rb, vB=vrb, n=len(pts))


def sc_units(study: str, pb: pd.DataFrame, mal: str) -> dict:
    """H138.5 units of one cohort: primary (malignant vs normal epithelium, host = stromal anchors) and companion
    (tumour vs adjacent fibroblasts, host = the other stromal compartments). Entries are None when not estimable."""
    out = dict(primary=None, fibroblast=None)
    comp = pb.compartment
    pre = list(comp.unique())
    anc, src = stromal_anchors(study, pre)
    if mal in pre and NORM in pre and len(anc) >= 1:
        h = host_of(pb, anc)
        if h is not None:
            u = sc_pair(pb, h, mal, NORM)
            if u:
                out["primary"] = u | dict(host_source=src, n_host=len(anc))
    fa, fsrc = stromal_anchors(study, pre, exclude_class=("fibroblast",))
    if "tumour:fibroblast" in pre and "adjacent_normal:fibroblast" in pre and len(fa) >= 1:
        h = host_of(pb, fa)
        if h is not None:
            u = sc_pair(pb, h, "tumour:fibroblast", "adjacent_normal:fibroblast")
            if u:
                out["fibroblast"] = u | dict(host_source=fsrc, n_host=len(fa))
    return out


# ---- H138.1: simulated cohorts at the observed design -------------------------------------------------------------

WINDOW_H, JITTER_H = (8.0, 18.0), 1.0  # CHOICE: the simulator's default hospital-hours window
C_GRID = (0.0, 0.1, 0.2, 0.3, 0.5)
C_ANCHOR = 0.8  # CHOICE: true coupling of the simulated anchors (A17 learned anchors: strongly coupled)
HOST_C, TARGET = "host", "MAL"


def q1_design(study: str) -> dict | None:
    """Design of a Q1 cohort read from the counts only (no coupling is computed): patients with the malignant
    compartment, anchors present (learned A17 rule), median depth of the malignant and of each anchor compartment,
    share of the malignant patients with each anchor, mean clock-gene fractions (the simulator's mesor)."""
    pb, mal = study_table(study)
    try:
        anc_all = [c for c in anchors(f"cancer_{study}", "learned") if c in set(pb.compartment)]
    except KeyError:
        return None
    # CHOICE: the malignant compartment is the target, not a host anchor (tumour epithelium is a learned anchor of some
    # cohorts: its fit is then the leave-one-out one, i.e. the host of the other anchors); >= 3 learned anchors
    # present as in A17 / A18.3, >= 2 of them host anchors
    anc = [c for c in anc_all if c != EPI]
    m = pb[pb.compartment == mal]
    n = m.donor.nunique()
    if len(anc_all) < 3 or len(anc) < 2 or n < 5:
        return None
    pts = set(m.donor)
    rows = pb[pb.compartment.isin(anc + [mal]) & pb.donor.isin(pts)]
    frac = rows[GENES].sum() / rows.depth.sum()
    return dict(study=study, mal=mal, n=n, anchors=anc,
                depth={c: float(pb[pb.compartment == c].depth.median()) for c in anc + [mal]},
                presence={c: float(pb[(pb.compartment == c) & pb.donor.isin(pts)].donor.nunique() / n) for c in anc},
                mesor=np.log(np.maximum(frac.to_numpy(float), 1e-9)))


def _coupled_variance(rng, n: int = 200_000) -> float:
    """Trace of the covariance of (cos, sin) of the host phase over the sampling window: the clock-plane variance S of
    a unit-amplitude coupled compartment, so that psi = S (1 - c) / (2 c) gives true coupling c (`on_host`)."""
    t = rng.uniform(*WINDOW_H, n) + rng.normal(0.0, JITTER_H, n)
    return float(1.0 - abs(np.mean(np.exp(1j * t / sh.HOURS))) ** 2)


S_PHASE = _coupled_variance(np.random.default_rng(0))


def simulate_design(d: dict, c: float, seed: int) -> pd.DataFrame:
    """One simulated cohort at the design `d` with true malignant coupling `c` (0: decoupled, uniform phase). CHOICE:
    depth is per compartment (median of the real compartment): counts are simulated at the largest depth and binomially
    thinned to each compartment's own depth (a Poisson sum thins to a Poisson sum). Anchors: amplitude 1, coupling
    C_ANCHOR, presence as observed; malignant: amplitude 1, presence 1."""
    wf = Waveform(WF.genes, WF.amp, WF.psi, d["mesor"], WF.sd, "design")
    psi = lambda cc: S_PHASE * (1 - cc) / (2 * cc)
    comps = {a: sim.Compartment(z=1.0, psi=psi(C_ANCHOR), presence=d["presence"][a]) for a in d["anchors"]}
    comps[TARGET] = sim.Compartment(z=1.0, psi=psi(c) if c > 0 else 0.0, coupled=c > 0)
    dmax = max(d["depth"].values())
    t, _ = sim.simulate(wf, comps, n=d["n"], window_h=WINDOW_H, jitter_h=JITTER_H, depth=dmax, seed=seed)
    rng = np.random.default_rng([seed, 1])
    dep = np.array([d["depth"][d["mal"] if k == TARGET else k] for k in t.compartment])
    p = (dep / dmax)[:, None]
    t[GENES] = rng.binomial(t[GENES].to_numpy().astype(np.int64), p)
    t["depth"] = dep
    return t.rename(columns={"resection": "donor"})


def fit_coupling(t: pd.DataFrame, d: dict, n_boot: int, seed: int) -> tuple[float, float]:
    """Anchored coupling of the malignant compartment and its donor-bootstrap SD (host refitted in each replicate)."""
    ds = cs.Dataset(t, anchors=d["anchors"])
    f = ds.fit(n_boot=n_boot, seed=seed, targets=[TARGET], gate=False)
    c = float(f.table.loc[TARGET, "coupling"])
    rep = f.replicates[TARGET].to_numpy(float)
    rep = rep[np.isfinite(rep)]
    return c, float(np.std(rep, ddof=1)) if len(rep) > 2 else np.nan
