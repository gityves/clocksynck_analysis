"""Helpers of A109 (coupling nulls). Data access (GTEx tissues, pseudo-clock gene sets) and statistics that are not
yet in `clocksync` ("local: not in clocksync"; proposed additions are listed in the README).

GTEx time of death and protected phenotypes are controlled access: functions here keep them in memory and return
aggregates only.
"""

from __future__ import annotations

import importlib
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_3" / "scripts" / "gtex_cues"))
import clocksync as cs  # noqa: E402
import _shared as S  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import gtex  # noqa: E402
from common.paths import GTEX_PHENOTYPES  # noqa: E402

P = importlib.import_module("clocksync.paired")  # the module (the package namespace holds the function `paired`)

HOURS = 24.0 / (2 * np.pi)
BODY = [t for t in S.TISSUES if t != "Pituitary" and not t.startswith("Brain")]  # the 12 A92 tissues
GROUPS = {"fast": (1, 2), "H3": (3,), "H4": (4,), "vent": (0,)}
MIN_N = 20  # shared donors per pair (A93)
GC_H = 17.7  # A27 D6c pull phase of glucocorticoids in ventilated donors (diag6c_gc_summary.csv: 17.68 h)
# D6 clock-output genes (compute_glucocorticoid.classes, "clock output (non-frame)"), excluded from pseudo-clock sets
CLOCK_OUTPUT = {"HLF", "NFIL3", "BHLHE40", "BHLHE41", "RORC", "NPAS2", "ARNTL2", "CRY2", "TSPAN4", "USP2", "LONRF3",
                "TEF", "NAMPT", "KLF10", "DDIT4L", "WEE1", "PPARA", "CIPC"}
WF = cs.gtex()
PAIRS = list(itertools.combinations(range(len(BODY)), 2))


# ---- data ------------------------------------------------------------------------------------------------------------

class Tissues:
    """The 12 tissues of GTEx (Hardy 0-4, time of death recorded, cultured cells excluded, first sample of a donor).

    donors: global donor index (>= 2 tissues; `subsample` random donors for checks); hardy, tod (controlled): arrays
    over donors; R[t]: row of donor d in tissue t (-1: absent); rows[t]: DataFrame of clock-gene counts + depth
    (library size); pseudo[t]: counts (rows x genes) of the genes of the pseudo-clock sets of tissue t; sets[t]: the
    pseudo-clock gene sets of tissue t (list of gene lists); hardy_of_row etc. via donors.
    """

    def __init__(self, n_sets: int, seed: int = 0, subsample: int | None = None):
        import anndata as ad

        a = ad.read_h5ad(gtex.LABELLED, backed="r")
        try:
            o = a.obs
            keep = (o.DTHHRDY.isin([0, 1, 2, 3, 4]) & o.tissue.astype(str).isin(BODY) & o.tod.notna()).to_numpy()
            sub = o.iloc[np.flatnonzero(keep)]
            subj = sub.subject.astype(str).to_numpy()
            ok = ~pd.Series(list(zip(subj, sub.tissue.astype(str)))).duplicated().to_numpy()
            sub, subj = sub[ok], subj[ok]
            nt = pd.Series(subj).value_counts()
            donors = sorted(nt.index[nt >= 2])
            if subsample:
                donors = sorted(np.random.default_rng(seed).permutation(donors)[:subsample])
            self.donors = np.asarray(donors)
            gi = {d: i for i, d in enumerate(donors)}
            first = sub.assign(subject=subj).drop_duplicates("subject").set_index("subject")
            d = first.reindex(donors)
            self.hardy = d.DTHHRDY.to_numpy(int)
            self.tod = d.tod.to_numpy(float) % 24
            self.group = np.array([next(g for g, h in GROUPS.items() if x in h) for x in self.hardy])
            self.R, self.rows, self.pseudo, self.sets, self.row_donor = {}, {}, {}, {}, {}
            names = np.asarray(a.var_names)
            clock = list(WF.genes)
            from clocksync.genes import gene_columns
            ccols = gene_columns(a.var_names, clock)
            for k, t in enumerate(BODY):
                m = (sub.tissue.astype(str) == t).to_numpy() & np.isin(subj, donors)
                ids = sub.index[m]
                pos = a.obs_names.get_indexer(ids)
                o_ = np.argsort(pos)
                pos, ids, sd = pos[o_], ids[o_], subj[m][o_]
                X = a.X[pos]
                X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
                lib = a.obs.library_size.to_numpy(float)[pos]
                lm = pd.Series(np.log(X / lib[:, None] * 1e6 + 0.5).mean(0), index=names)
                lm = lm[~lm.index.isin(CLOCK_OUTPUT - set(clock))]
                lm = lm[~lm.index.duplicated()]
                wfs = P.pseudo_waveforms(lm, n_sets, np.random.default_rng(1000 * (seed + 1) + k), WF)
                self.sets[t] = [list(w.genes) for w in wfs]
                u = sorted({g for s in self.sets[t] for g in s})
                col = {g: i for i, g in enumerate(names)}
                self.pseudo[t] = pd.DataFrame(X[:, [col[g] for g in u]].astype(np.float32), columns=u)
                x = pd.DataFrame(X[:, ccols], columns=clock)
                x["depth"] = lib
                self.rows[t] = x
                r = np.full(len(donors), -1)
                r[[gi[s] for s in sd]] = np.arange(len(sd))
                self.R[t] = r
                self.row_donor[t] = sd
        finally:
            a.file.close()
        self.n_sets = n_sets
        self.Z = {}  # kind -> tissue -> (K, n_rows) complex states; V: Poisson variance (trace)

    def build_states(self, acro_rng: np.random.Generator, n_acro: int, real: bool = True):
        """Clock-plane states of every tissue for the real clock genes (K = 1), the pseudo-clock sets (K = n_sets) and
        permuted-acrophase waveforms (K = n_acro); one projector per tissue and set."""
        self.Z = {k: {} for k in ("clock", "pseudo", "acro")}
        self.V = {k: {} for k in ("clock", "pseudo", "acro")}
        for t in BODY:
            x, dep = self.rows[t], self.rows[t].depth.to_numpy(float)
            if real:
                z, v = _state_set(x, [WF], dep)
                self.Z["clock"][t], self.V["clock"][t] = z, v
                ws = [WF.permuted(acro_rng) for _ in range(n_acro)]
                z, v = _state_set(x, ws, dep)
                self.Z["acro"][t], self.V["acro"][t] = z, v
            ps = [cs.Waveform(tuple(g), WF.amp, WF.psi, None, None, f"p{j}")
                  for j, g in enumerate(self.sets[t])]
            z, v = [], []
            for w, g in zip(ps, self.sets[t]):
                xx = self.pseudo[t][list(g)].astype(float)
                L = P.projector(xx, w, counts=True)
                b, vv = P.states(xx, L, w, depth=dep)
                z.append(b[:, 0] + 1j * b[:, 1])
                v.append(vv)
            self.Z["pseudo"][t], self.V["pseudo"][t] = np.array(z), np.array(v)

    def group_donors(self, g: str) -> np.ndarray:
        return np.flatnonzero(self.group == g)


def _state_set(x: pd.DataFrame, waveforms, depth):
    z, v = [], []
    for w in waveforms:
        L = P.projector(x[list(w.genes)], w, counts=True)
        b, vv = P.states(x, L, w, depth=depth)
        z.append(b[:, 0] + 1j * b[:, 1])
        v.append(vv)
    return np.array(z), np.array(v)


def procedure_kits() -> dict:
    """Collection kit (SMCAT of the GTEx sample attributes, controlled access) of every tissue: the majority kit."""
    f = GTEX_PHENOTYPES / "phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt"
    d = pd.read_csv(f, sep="\t", comment="#", usecols=["SMTSD", "SMCAT"], low_memory=False)
    d = d[d.SMTSD.isin(BODY)]
    return d.groupby("SMTSD").SMCAT.agg(lambda s: s.value_counts().index[0]).to_dict()


# ---- coherence C (A93 form) ------------------------------------------------------------------------------------------

def pair_corr(T: Tissues, kind: str, D: np.ndarray, R=None, pair_idx=None) -> np.ndarray:
    """(K, n_pairs) complex Poisson-corrected correlations E[a b*] / sqrt(E|a|^2 E|b|^2) over donors D (global
    indices, duplicates allowed) for all K state sets of `kind`; NaN where fewer than MIN_N donors are shared.
    Centred over the shared donors, as A93 `pair_corr`."""
    R = T.R if R is None else R
    Z, V = T.Z[kind], T.V[kind]
    K = next(iter(Z.values())).shape[0]
    pairs = PAIRS if pair_idx is None else [PAIRS[i] for i in pair_idx]
    out = np.full((K, len(pairs)), np.nan + 0j)
    rr = {t: R[t][D] for t in BODY}
    for j, (a, b) in enumerate(pairs):
        ta, tb = BODY[a], BODY[b]
        m = (rr[ta] >= 0) & (rr[tb] >= 0)
        n = int(m.sum())
        if n < MIN_N:
            continue
        ia, ib = rr[ta][m], rr[tb][m]
        A, B = Z[ta][:, ia], Z[tb][:, ib]
        A, B = A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)
        k = n - 1
        ea = (abs(A) ** 2).sum(1) / k - V[ta][:, ia].mean(1)
        eb = (abs(B) ** 2).sum(1) / k - V[tb][:, ib].mean(1)
        out[:, j] = (A * B.conj()).sum(1) / k / np.sqrt(np.maximum(ea, 1e-12) * np.maximum(eb, 1e-12))
    return out


def coherence(T: Tissues, kind: str, draws: list, R=None, pair_idx=None) -> np.ndarray:
    """(K,) C: mean over pairs of |ccorr|, median over the donor draws (equal-n subsamples)."""
    c = np.array([np.nanmean(abs(pair_corr(T, kind, D, R, pair_idx)), axis=1) for D in draws])
    return np.nanmedian(c, axis=0)


def draws_equal_n(donors: np.ndarray, n: int, n_sub: int, rng) -> list:
    if len(donors) <= n:
        return [donors]
    return [rng.choice(donors, n, replace=False) for _ in range(n_sub)]


def shuffled_R(T: Tissues, donors: np.ndarray, rng) -> dict:
    """Donor labels shuffled within tissue among the donors of one group (pairing destroyed, data kept)."""
    R = {t: r.copy() for t, r in T.R.items()}
    for t in BODY:
        idx = donors[R[t][donors] >= 0]
        R[t][idx] = rng.permutation(R[t][idx])
    return R


# ---- group-centred cross-moments (matrix form) ------------------------------------------------------------------------

def donor_matrix(T: Tissues, D: np.ndarray, kind: str = "clock", k: int = 0):
    """(Z, M, Vm): donor x tissue complex states (zero where absent), presence mask, Poisson variance for donors D,
    states centred within tissue over D."""
    n = len(D)
    Z = np.zeros((n, len(BODY)), complex)
    M = np.zeros((n, len(BODY)))
    Vm = np.zeros((n, len(BODY)))
    for j, t in enumerate(BODY):
        r = T.R[t][D]
        m = r >= 0
        z = T.Z[kind][t][k][r[m]]
        Z[m, j] = z - z.mean()
        M[m, j] = 1
        Vm[m, j] = T.V[kind][t][k][r[m]]
    return Z, M, Vm


def cross_moments(Z: np.ndarray, M: np.ndarray):
    """(E[z_a conj z_b], E[z_a z_b], N) tissue x tissue matrices over shared donors (Z zero where absent)."""
    N = M.T @ M
    k = np.maximum(N - 1, 1)
    return (Z.T @ Z.conj()) / k, (Z.T @ Z) / k, N


def pair_phase(H: np.ndarray) -> np.ndarray:
    """Phase (rad) of E[z_a conj z_b] of the 66 pairs a < b (a minus b)."""
    return np.array([np.angle(H[a, b]) for a, b in PAIRS])


def permute_labels(H: np.ndarray, perm: np.ndarray) -> np.ndarray:
    """Pair phases after the tissue labels are permuted jointly: phase of E[z_p(a) conj z_p(b)]."""
    return np.angle(H[np.ix_(perm, perm)])[np.triu_indices(len(perm), 1)]


def fisher_lee_fast(a: np.ndarray, b: np.ndarray) -> float:
    return ST.fisher_lee(a, b)


def folded_distance(theta: np.ndarray) -> float:
    """Median distance (h) of the pair phases from the nearest of 0 and 12 h (|theta| folded to [0, 12 h])."""
    x = np.abs(theta) * HOURS
    return float(np.median(np.minimum(x, 12 - x)))


def folded_median(theta: np.ndarray) -> float:
    """Median absolute pair phase (h) folded to [0, 12 h]."""
    return float(np.median(np.abs(theta) * HOURS))


# ---- host descriptives and circularity ------------------------------------------------------------------------------

def hermitian(Z, M, Vm):
    """Noise-corrected Hermitian covariance of the tissues (diagonal: corrected variance)."""
    H, Pc, N = cross_moments(Z, M)
    var = (abs(Z) ** 2).sum(0) / np.maximum(N.diagonal() - 1, 1) - (Vm * M).sum(0) / np.maximum(M.sum(0), 1)
    H = H.copy()
    H[np.diag_indices(len(H))] = var
    return H, Pc, N


def host_frame(H: np.ndarray):
    """Leading eigenvector of the Hermitian covariance, phases anchored so the weights sum to a positive real, and the
    leading eigenvalue share."""
    w, U = np.linalg.eigh(H)
    u = U[:, -1]
    u = u * np.exp(-1j * np.angle(u.sum()))
    return u, float(w[-1] / max(w[w > 0].sum(), 1e-12))


def host_states(Z, M, u):
    """Host state of every donor: weighted mean of the tissue states aligned with the leading eigenvector."""
    num = (Z * M * np.conj(u)).sum(1)
    den = (M * np.abs(u)).sum(1)
    return num / np.maximum(den, 1e-12)


def axis_angle_h(w: np.ndarray) -> float:
    """Major axis of a cloud of complex states: angle in h, mod 12."""
    x = np.c_[w.real, w.imag]
    x = x - x.mean(0)
    _, v = np.linalg.eigh(x.T @ x)
    return float((np.arctan2(v[1, -1], v[0, -1]) * HOURS) % 12)


def circ_ratio(Z, M) -> float:
    """Pooled |E[z_a z_b]| / |E[z_a conj z_b]| over pairs (mean magnitudes): 0 for rotating clocks with uniform phase,
    1 for a pure displacement along a fixed direction."""
    H, Pc, _ = cross_moments(Z, M)
    iu = np.triu_indices(len(H), 1)
    return float(np.mean(abs(Pc[iu])) / np.mean(abs(H[iu])))


def simulate_group(Z, M, Vm, c: float, f: float, gamma_h: float, rng, kind: str = "mix"):
    """One simulated donor x tissue state matrix with the observed design (donors, tissues, Poisson variance), the
    observed per-tissue corrected variance a_t^2 and tissue offsets, and the observed coupling c (local: not in
    clocksync).  z = a_t [ sqrt(c f) s e^{i gamma} + sqrt(c (1 - f)) e^{i (phi + delta_t)} + sqrt(1 - c) xi ] + noise
    with s ~ N(0, 1) (shared push), phi ~ U[0, 2 pi) (rotating shared clock), xi complex Gaussian, E|z|^2 = a_t^2.
    f = 0: rotating clocks at coupling c; f = 1: pure displacement of the observed size along gamma."""
    H, _, _ = hermitian(Z, M, Vm)
    a = np.sqrt(np.maximum(H.diagonal().real, 1e-12))
    u, _ = host_frame(H)
    delta = np.angle(u)
    n, T = Z.shape
    g = 2 * np.pi * gamma_h / 24
    s = rng.standard_normal(n)[:, None]
    phi = rng.uniform(0, 2 * np.pi, n)[:, None]
    xi = (rng.standard_normal((n, T)) + 1j * rng.standard_normal((n, T))) / np.sqrt(2)
    z = a * (np.sqrt(c * f) * s * np.exp(1j * g) + np.sqrt(c * (1 - f)) * np.exp(1j * (phi + delta))
             + np.sqrt(1 - c) * xi)
    noise = np.sqrt(Vm / 2) * (rng.standard_normal((n, T)) + 1j * rng.standard_normal((n, T)))
    z = (z + noise) * M
    mu = z.sum(0) / np.maximum(M.sum(0), 1)
    return (z - mu * M) * M


# ---- non-temporal shared direction ------------------------------------------------------------------------------------

def tod_residual_states(T: Tissues, D: np.ndarray):
    """Per tissue, states (b0, b1) of the donors D with (cos, sin) of the time of death regressed out (controlled:
    in memory only). Returns donor x tissue complex residuals, mask."""
    n = len(D)
    Z = np.zeros((n, len(BODY)), complex)
    M = np.zeros((n, len(BODY)))
    w = 2 * np.pi * T.tod[D] / 24
    for j, t in enumerate(BODY):
        r = T.R[t][D]
        m = r >= 0
        z = T.Z["clock"][t][0][r[m]]
        X = np.c_[np.ones(m.sum()), np.cos(w[m]), np.sin(w[m])]
        Y = np.c_[z.real, z.imag]
        e = Y - X @ np.linalg.lstsq(X, Y, rcond=None)[0]
        Z[m, j] = e[:, 0] + 1j * e[:, 1]
        M[m, j] = 1
    return Z, M


def cross_cov_matrix(Z: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Sum over tissue pairs a < b of the real 2 x 2 cross-covariance E[b_a b_b^T] of the clock-plane states,
    symmetrised (donors shared by the pair, states centred within tissue over the donors given)."""
    out = np.zeros((2, 2))
    for a, b in PAIRS:
        m = (M[:, a] > 0) & (M[:, b] > 0)
        if m.sum() < MIN_N:
            continue
        A = np.c_[Z[m, a].real, Z[m, a].imag]
        B = np.c_[Z[m, b].real, Z[m, b].imag]
        A, B = A - A.mean(0), B - B.mean(0)
        C = A.T @ B / (m.sum() - 1)
        out += (C + C.T) / 2
    return out


def share_along(Mx: np.ndarray, angle_h: float) -> float:
    """Share of the summed cross-covariance trace along the axis `angle_h` (h, mod 12) of the clock plane."""
    th = angle_h / HOURS
    d = np.array([np.cos(th), np.sin(th)])
    o = np.array([-d[1], d[0]])
    al, ort = d @ Mx @ d, o @ Mx @ o
    return float(al / (al + ort)) if (al + ort) != 0 else np.nan
