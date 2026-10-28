"""A140 Q1: simulator of the preprint's design, the published CYCLOPS/DIP pipeline and the clocksync gate.

Registered in paper/docs/registered_experiments.md (A140, Q1). Choices not in the registration are listed in README.md
("Q1", "Choices not in the registration") and fixed here before any scoring.

Design: 76 patients x 3 compartments (CD4, CD8, myeloid) x p = 2,000 genes at single-cell pseudobulk depth.
Clock part: `clocksync.simulate` (GTEx waveform of the 14 clock genes, rhythmic followers at GTEx gene phases).
Non-clock part: K patient factors shared by the three compartments, applied to the counts by binomial thinning
(`clocksync.simulate.thin_clock`, the package's count-level manipulation), loading on 10-30 % of the genes, with the
variance of the leading non-clock principal components of real immune pseudobulks (`calibrate`).
"""
from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.simulate import Compartment, thin_clock
from clocksync.waveform import Waveform

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
Q2 = HERE / "output" / "q2"
OUT = HERE / "output" / os.environ.get("A140_OUT", "q1")  # q1 revision: A140_OUT=q1rev
CALIB = HERE / "output" / "q1" / "tumour_immune_genome" / "calibration_lambda.csv"  # used when A140_CALIB=1
COMPS = ["CD4_T", "CD8_T", "myeloid"]
N_PATIENTS, P_GENES = 76, 2000
COUPLING = 0.8
HOURS = 24 / (2 * np.pi)
GTEX_COS = PAPER / "results" / "figure_2" / "output" / "gene_rhythms" / "gtex_cosinor.csv.gz"
ARMS = {"A": dict(clock=True, coupled=True, window=(8.0, 18.0)),
        "B": dict(clock=False, coupled=True, window=(8.0, 18.0)),
        "C": dict(clock=True, coupled=False, window=(8.0, 18.0)),
        "D": dict(clock=True, coupled=True, window=(0.0, 24.0))}
ARM_INDEX = {a: i for i, a in enumerate("ABCD")}
SEED_SIZES = (50, 150, 500, 1500)
import os
N_MODELS = int(os.environ.get("A140_N_MODELS", 10))  # CYCLOPS train_n_models: 10 (template 80; cost, see README)
CORE = os.environ.get("A140_CORE") == "1"  # 2026-10-28 note: gate / primary with the core-clock seed; outputs <step>_core/
SUFFIX = "_core" if CORE else ""
SEED_MTH_GENE = 1000  # CYCLOPS expression cutoff rank: half of the 2,000 genes (template 10,000 of ~20,000)
CLOCK = list(cs.gtex().genes)


# ---- the gene panel and the calibration of the non-clock factor (from GSE266919 pseudobulks) --------------------------

def build_panel() -> dict:
    """Panel (2,000 genes), per-gene mesor (log fraction of the mean real immune pseudobulk), GTEx waveform of the
    clock genes and followers, depth, and the non-clock PC variances of the primary set (44 Pre-treatment samples)."""
    sm = pd.read_csv(Q2 / "samples.csv", index_col=0)
    pre = sm.index[sm.Group == "Pre-treatment"]
    pbs = {}
    for c in COMPS:
        p = pd.read_csv(Q2 / f"pb_{c}.csv.gz", index_col=0)
        p = p[p.n_cells >= 20]
        pbs[c] = p.drop(columns="n_cells")
    genes = sorted(set.intersection(*[set(p.columns) for p in pbs.values()]))
    frac = np.mean([pbs[c][genes].sum(0).to_numpy() / pbs[c][genes].to_numpy().sum() for c in COMPS], axis=0)
    frac = pd.Series(frac, index=genes)
    depth = float(np.median(np.concatenate([p.sum(axis=1).to_numpy() for p in pbs.values()])))

    # followers: genes with GTEx p < 0.01 in >= 3 of 13 tissues; phase and amplitude pooled over those tissues
    g = pd.read_csv(GTEX_COS)
    g = g[g.p < 0.01]
    cnt = g.groupby("gene").tissue.nunique()
    fol = cnt.index[cnt >= 3]
    g = g[g.gene.isin(fol)]
    z = g.assign(c=g.b_cos + 1j * g.b_sin).groupby("gene").c.mean()
    amp_f = g.groupby("gene").amp.mean()
    fol = [x for x in fol if x in frac.index and frac[x] >= 1e-6 and x not in CLOCK]
    clock = [x for x in CLOCK if x in frac.index]
    assert len(clock) == 14
    rest = frac.drop(index=clock + fol).sort_values(ascending=False).index[: P_GENES - len(clock) - len(fol)]
    panel = clock + fol + list(rest)
    w0 = cs.gtex()
    amp = np.zeros(len(panel))
    psi = np.zeros(len(panel))
    for i, x in enumerate(panel):
        if x in clock:
            j = w0.genes.index(x)
            amp[i], psi[i] = w0.amp[j], w0.psi[j]
        elif x in fol:
            amp[i], psi[i] = amp_f[x], np.angle(z[x]) % (2 * np.pi)
    mesor = np.log(frac[panel].to_numpy())

    # calibration: PCs of the log fractions of the non-clock genes, compartments centred, stacked (44 primary samples)
    nc = [x for x in panel if x not in clock]
    rows = []
    for c in COMPS:
        p = pbs[c].loc[pbs[c].index.intersection(pre)]
        y = np.log((p[nc].to_numpy() + 0.5) / p.sum(axis=1).to_numpy()[:, None])
        rows.append(y - y.mean(0))
    Y = np.vstack(rows)
    s = np.linalg.svd(Y, compute_uv=False)
    lam = s ** 2 / (len(Y) - 1)
    if os.environ.get("A140_CALIB") == "1":  # revision: eigenvalues of the tumour immune genome-wide calibration
        lam = np.concatenate([pd.read_csv(CALIB).lam.to_numpy(), lam[3:5]])
    return dict(genes=panel, n_clock=len(clock), n_follow=len(fol), amp=amp, psi=psi, mesor=mesor, depth=depth,
                lam=lam[:5], lam_total=float((Y ** 2).sum() / (len(Y) - 1)), n_rows=len(Y))


def load_panel() -> dict:
    f = OUT / "panel.npz"
    if not f.exists():
        OUT.mkdir(parents=True, exist_ok=True)
        P = build_panel()
        np.savez(f, **{k: (np.array(v) if k != "genes" else np.array(v, dtype=str)) for k, v in P.items()})
        pd.DataFrame({"gene": P["genes"], "amp": P["amp"], "acro_h": P["psi"] * HOURS, "mesor": P["mesor"]}).to_csv(
            OUT / "panel_genes.csv", index=False)
    z = np.load(f, allow_pickle=False)
    P = {k: z[k] for k in z.files}
    P["genes"] = [str(g) for g in P["genes"]]
    for k in ("n_clock", "n_follow", "n_rows"):
        P[k] = int(P[k])
    P["depth"] = float(P["depth"])
    return P


# ---- simulator -----------------------------------------------------------------------------------------------------

def loadings(P: dict, K: int, seed: int) -> np.ndarray:
    """G x K loadings: each axis loads (normal) on a random 10-30 % of the genes, unit norm."""
    rng = np.random.default_rng(seed)
    G = len(P["genes"])
    L = np.zeros((G, K))
    for k in range(K):
        n = int(round(rng.uniform(0.10, 0.30) * G))
        idx = rng.choice(G, n, replace=False)
        v = rng.normal(size=n)
        L[idx, k] = v / np.linalg.norm(v)
    return L


def psi_for(coupling: float, window=(8.0, 18.0), jitter: float = 1.0, z: float = 1.0) -> float:
    """Clock-plane noise psi giving true coupling `coupling` at the sampling spread of `window`
    (psi = z^2 V (1 - c) / (2 c), V = trace of the covariance of (cos, sin) of the sampled time; as the A19/A21 simulations)."""
    t = np.random.default_rng(0).uniform(*window, 200000) + np.random.default_rng(1).normal(0, jitter, 200000)
    V = np.trace(np.cov(np.column_stack([np.cos(t / HOURS), np.sin(t / HOURS)]).T))
    return z ** 2 * V * (1 - coupling) / (2 * coupling)


def simulate_cohort(P: dict, arm: str, K: int, seed: int, L: np.ndarray, n: int = N_PATIENTS, window=None,
                    scale: float = 1.0):
    """(long table of the 3 compartments, true clock time t (h), factor scores S (n x K))."""
    a = ARMS[arm]
    window = a["window"] if window is None else window
    G = len(P["genes"])
    amp = P["amp"] if a["clock"] else np.zeros(G)
    wf = Waveform(tuple(P["genes"]), amp, P["psi"], P["mesor"], None, f"A140_{arm}")
    psi = psi_for(COUPLING, ARMS["A"]["window"])
    comps = {c: Compartment(z=1.0, psi=psi, coupled=a["coupled"]) for c in COMPS}
    tab, t = cs.simulate(wf, comps, n=n, window_h=window, jitter_h=1.0, depth=P["depth"], seed=seed)
    rng = np.random.default_rng(seed + 7_000_003)
    sd = np.sqrt(P["lam"][:K])
    S = rng.normal(size=(n, K)) * sd * scale  # scale != 1 only in the post hoc diagnostic
    pid = tab.resection.str[1:].astype(int).to_numpy()
    counts = tab[P["genes"]].to_numpy(float)
    th = thin_clock(counts, L, S[pid], rng)
    tab = tab.copy()
    tab[P["genes"]] = th
    return tab, t, S


# ---- the published pipeline ----------------------------------------------------------------------------------------

def cpm(tab: pd.DataFrame, genes) -> pd.DataFrame:
    x = tab[list(genes)].to_numpy(float)
    return pd.DataFrame(x / x.sum(1, keepdims=True) * 1e6, columns=list(genes))


def wrap(a):
    return (np.asarray(a) + np.pi) % (2 * np.pi) - np.pi


def circ_corr(a, b) -> float:
    """Jammalamadaka-SenGupta circular correlation (sign-sensitive)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ma, mb = np.angle(np.exp(1j * a).sum()), np.angle(np.exp(1j * b).sum())
    sa, sb = np.sin(a - ma), np.sin(b - mb)
    den = np.sqrt((sa ** 2).sum() * (sb ** 2).sum())
    return float((sa * sb).sum() / den) if den > 0 else np.nan


def rayleigh_R(a) -> float:
    a = np.asarray(a, float)
    return float(np.abs(np.exp(1j * a).mean())) if len(a) else np.nan


def dip(phases: np.ndarray, max_iter: int = 5) -> np.ndarray:
    """Donor internal phase: circular mean of a donor's phases (n x 3), iteratively dropping phases more than pi/2
    from the mean (at most `max_iter` iterations; a drop that would empty the set is not applied)."""
    out = np.empty(len(phases))
    for i, ph in enumerate(phases):
        keep = np.ones(len(ph), bool)
        m = np.angle(np.exp(1j * ph).sum())
        for _ in range(max_iter):
            far = np.abs(wrap(ph - m)) > np.pi / 2
            new = keep & ~far
            if new.sum() == 0 or (new == keep).all():
                break
            keep = new
            m = np.angle(np.exp(1j * ph[keep]).sum())
        out[i] = m % (2 * np.pi)
    return out


def dip_bins(d: np.ndarray, nbin: int = 6):
    """Six DIP bins of equal size: the circle is cut at the largest gap between sorted DIPs, the donors ranked along the
    resulting arc and split in `nbin` groups; bin time = circular mean DIP of the group. Returns bin index and times (rad)."""
    o = np.argsort(d)
    ds = d[o]
    gaps = np.diff(np.concatenate([ds, [ds[0] + 2 * np.pi]]))
    start = (int(np.argmax(gaps)) + 1) % len(d)
    order = np.concatenate([o[start:], o[:start]])
    b = np.empty(len(d), int)
    for k, idx in enumerate(np.array_split(order, nbin)):
        b[idx] = k
    times = np.array([np.angle(np.exp(1j * d[b == k]).sum()) % (2 * np.pi) for k in range(nbin)])
    return b, times


def cosinor_bins(Y: np.ndarray, bin_t: np.ndarray, nbin_min: int = 4):
    """Cosinor on bin means: Y (bins x genes) at bin times (rad). F-test (2, nb - 3) per gene. Returns p, acrophase (rad)."""
    from scipy import stats
    nb = len(bin_t)
    if nb < nbin_min:
        return None, None
    D = np.column_stack([np.ones(nb), np.cos(bin_t), np.sin(bin_t)])
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    rss = ((Y - D @ B) ** 2).sum(0)
    tss = ((Y - Y.mean(0)) ** 2).sum(0)
    df = nb - 3
    F = ((tss - rss) / 2) / np.maximum(rss / df, 1e-300) if df > 0 else np.full(Y.shape[1], np.nan)
    p = stats.f.sf(F, 2, df) if df > 0 else np.full(Y.shape[1], np.nan)
    return p, np.arctan2(B[2], B[1]) % (2 * np.pi)


def pipeline_metrics(tab: pd.DataFrame, ph: pd.DataFrame, P: dict, t_true: np.ndarray) -> dict:
    """Readouts of the published pipeline from the CYCLOPS phases `ph` (index: rows of `tab`, column phase_raw)."""
    res = tab.resection.to_numpy()
    donors = sorted(set(res))
    phase = {c: pd.Series(ph.loc[:, "phase_raw"].to_numpy()[(tab.compartment == c).to_numpy()],
                          index=res[(tab.compartment == c).to_numpy()]).reindex(donors).to_numpy() for c in COMPS}
    M = np.column_stack([phase[c] for c in COMPS])
    pairs = [(0, 1), (0, 2), (1, 2)]
    cc = [circ_corr(M[:, i], M[:, j]) for i, j in pairs]
    t_rad = (t_true % 24) / HOURS
    vs_truth = [abs(circ_corr(M[:, i], t_rad)) for i in range(3)]
    d = dip(M)
    b, bin_t = dip_bins(d, 6)
    occupied = np.arange(len(bin_t))
    Ycpm = np.log2(cpm(tab, P["genes"]).to_numpy() + 1)
    out = dict(cc_cd4_cd8=cc[0], cc_cd4_myel=cc[1], cc_cd8_myel=cc[2], cc_median=float(np.nanmedian(cc)),
               truth_r_cd4=vs_truth[0], truth_r_cd8=vs_truth[1], truth_r_myel=vs_truth[2],
               truth_r_median=float(np.median(vs_truth)), dip_truth_r=abs(circ_corr(d, t_rad)),
               dip_R=rayleigh_R(d), dip_cover=len(np.unique(np.minimum((d / (2 * np.pi) * 24).astype(int), 23))) / 24,
               n_bins=len(occupied))
    n_rhy, accR, orders = [], [], []
    cmask = {c: (tab.compartment == c).to_numpy() for c in COMPS}
    pos = {r: i for i, r in enumerate(donors)}
    bin_of = np.array([b[pos[r]] for r in res])
    for c in COMPS:
        m = cmask[c]
        Yb = np.vstack([Ycpm[m & (bin_of == k)].mean(0) for k in occupied])
        p, acro = cosinor_bins(Yb, bin_t)
        if p is None:
            n_rhy.append(0)
            accR.append(np.nan)
            continue
        called = p < 0.05
        n_rhy.append(int(called.sum()))
        accR.append(rayleigh_R(acro[called]) if called.any() else np.nan)
    out.update(n_rhy_cd4=n_rhy[0], n_rhy_cd8=n_rhy[1], n_rhy_myel=n_rhy[2], n_rhy_min=int(min(n_rhy)),
               acro_R_cd4=accR[0], acro_R_cd8=accR[1], acro_R_myel=accR[2], acro_R_mean=float(np.nanmean(accR)) if not np.all(np.isnan(accR)) else np.nan)
    # clock-gene order: acrophases (cosinor on DIP-bin means pooled over the three compartments) against the GTEx waveform,
    # |mean exp(i(s acro - psi_gtex))| maximised over the direction s (rotation-free, reflection-free)
    full = np.log2(cpm(tab, P["genes"][:14]).to_numpy() + 1)
    Yb = np.vstack([np.mean([full[cmask[c] & (bin_of == k)].mean(0) - full[cmask[c]].mean(0) for c in COMPS], axis=0)
                    for k in occupied])
    _, acro = cosinor_bins(Yb, bin_t)
    if acro is not None:
        w = cs.gtex()
        psi_g = np.array([w.psi[w.genes.index(g)] for g in P["genes"][:14]])
        out["clock_order_R"] = float(max(abs(np.exp(1j * (s * acro - psi_g)).mean()) for s in (1, -1)))
    else:
        out["clock_order_R"] = np.nan
    return out


# ---- seeds from the simulated timed reference ------------------------------------------------------------------------

def seed_ranking(P: dict, K: int, L: np.ndarray, seed: int, scale: float = 1.0) -> list[str]:
    """Genes ranked by their strongest association with the non-clock factor or the clock in a simulated timed reference
    (arm D design: coupled clock, 24 h, same factor loadings; 76 patients): max(multiple R with (cos, sin) of the true
    time, |r| with each factor score), computed on log CPM with compartments centred. Clock genes included in the ranking."""
    tab, t, S = simulate_cohort(P, "D", K, seed + 55_000_000, L, scale=scale)
    X = np.log2(cpm(tab, P["genes"]).to_numpy() + 1)
    comp = tab.compartment.to_numpy()
    for c in COMPS:
        X[comp == c] -= X[comp == c].mean(0)
    pid = tab.resection.str[1:].astype(int).to_numpy()
    tt = t[pid] / HOURS
    D = np.column_stack([np.cos(tt), np.sin(tt)])
    D = D - D.mean(0)
    B, *_ = np.linalg.lstsq(D, X, rcond=None)
    ss = (X ** 2).sum(0) + 1e-12
    r_clock = np.sqrt(np.maximum(((D @ B) ** 2).sum(0) / ss, 0))
    Sc = S[pid] - S[pid].mean(0)
    r_fac = np.zeros(X.shape[1])
    for k in range(K):
        s = Sc[:, k]
        r_fac = np.maximum(r_fac, np.abs(s @ X) / (max(np.linalg.norm(s), 1e-12) * np.sqrt(ss)))
    score = np.maximum(r_clock, r_fac)
    return [P["genes"][i] for i in np.argsort(-score)]


def seed_set(P: dict, ranking: list[str], n: int, core: bool) -> list[str]:
    """Top `n` of the ranking without the 14 core clock genes, or the 14 core clock genes plus the top n - 14 others."""
    others = [g for g in ranking if g not in set(CLOCK)]
    return (list(CLOCK) + others[: n - len(CLOCK)]) if core else others[:n]


def run_fit(tab: pd.DataFrame, P: dict, seeds: list[str], tmp: Path, threads: int, seed: int = 0) -> pd.DataFrame:
    """CYCLOPS 2.0 (port) on the pooled pseudobulks with the cell type as covariate. Returns phase_raw, phase (rad) per row."""
    from common.benchmark.cyclops import run_cyclops
    X = cpm(tab, P["genes"])
    ids = [f"{r}_{c}" for r, c in zip(tab.resection, tab.compartment)]
    X.index = ids
    cov = pd.DataFrame({"celltype": tab.compartment.to_numpy()}, index=ids)
    r = run_cyclops(X.T, seeds, tmp, seed=seed, threads=threads, covariates=cov, seed_mth_Gene=SEED_MTH_GENE, train_n_models=N_MODELS)
    return r.reset_index(drop=True)


def clocksync_gate(tab: pd.DataFrame, seed: int = 0):
    """Free host on the three compartments, gate at the package default (1,000 permutations)."""
    meta = ["resection", "compartment", "depth"]
    ds = cs.Dataset(tab[meta + CLOCK], anchors=None)
    g = ds.gate(seed=seed)
    return dict(gate_pass=bool(g.passes), gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase,
                gate_clock_share=g.clock_share, gate_share_null=g.share_null_mean)


def cell_dir(name: str) -> Path:
    d = Path(tempfile.mkdtemp(prefix=f"a140_{name}_"))
    return d


def cleanup(d: Path):
    shutil.rmtree(d, ignore_errors=True)
