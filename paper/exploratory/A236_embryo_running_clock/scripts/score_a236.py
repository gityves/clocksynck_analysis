"""A236 scoring (registered in `paper/docs/registered_experiments.md`, commit d51cc24): is the time-of-day state of
the late mouse embryo (Qiu 2024, E15.0-E18.75) a running clock?

Design: 16 embryos, 4 batches, one embryo per batch at 0, 6, 12 and 18 h; cell types by the A108 rule. Timed model per
cell type: batch intercepts, a linear within-batch drift (t - 9 h) and cos wt, sin wt.
H236.1  gene content: common cos / sin coefficients (a_g, b_g) of the gene-centred log fractions of the 14 clock
        genes; S1 = share of (a, b) in the clock plane. Null (power note, before scoring): S1 of plasmode data sets
        without a time-locked signal, on the cell types present in all embryos; the registered acrophase-permutation
        null is anti-conservative and is reported beside it. H236.1b: without PER1 and PER2, and the H108.1 statistic
        without them.
H236.2  rotation: polarisation rho of the time-locked clock-plane component (1 adult-direction rotation, 0 axis),
        leave-one-batch-out jackknife t, one-sided, against the simulated critical value of the power note.
Both tests are scored on the cell types present in all 16 embryos, where the plasmode nulls are defined.
H236.3  specificity: H108.1's R against 200 expression-matched and 200 responsive pseudo-clock sets
        (`paired.pseudo_waveforms`) on the genome-wide pseudobulks (`build_genome.py`).

    python score_a236.py power [--n-sim 2000]            plasmode size and power; prints no observed statistic
    python score_a236.py score                           H236.1, H236.1b, H236.2, descriptive and sensitivity
    python score_a236.py companions [--jobs 4] [--n-pseudo 200] [--check]
"""
import argparse
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core.plane import gene_weights, log_fraction, projector  # noqa: E402
from clocksync._core.study import plane  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.timed import W  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

OUT = HERE.parents[1] / "output"
GENOME = DATA / "natexp" / "qiu2024" / "genome_late" / "pseudobulk_celltype.h5ad"
SEED = 236
WF = cs.gtex()
GENES = list(WF.genes)
GC_CLOCK = ["PER1", "PER2"]  # clock genes that glucocorticoids induce directly
ARMS = {"BMAL1 arm": ["ARNTL", "NPAS2", "CLOCK"],
        "E-box / D-box arm": ["PER1", "PER2", "PER3", "NR1D1", "NR1D2", "DBP", "TEF", "CIART"]}
NAMED = ["Hepatocytes", "Cardiomyocytes"]


# ---- data -------------------------------------------------------------------------------------------------------

@dataclass
class Frame:
    """Late-window arrays: Y (embryo x cell type x gene) log fractions, NaN where the cell type is absent."""
    win: pd.DataFrame
    donors: pd.DataFrame
    comps: list
    t: np.ndarray  # clock time (h)
    batch: np.ndarray
    Y: np.ndarray
    L: np.ndarray  # cell type x 2 x gene: clock-plane projector
    wg: np.ndarray  # cell type x gene: `gene_weights`

    @property
    def Z(self):
        """Clock-plane states as complex numbers (embryo x cell type)."""
        b = np.einsum("kcg,skg->skc", self.L, self.Y)
        return b[..., 0] + 1j * b[..., 1]


def load() -> Frame:
    pb, donors = data.load("mouse_qiu")
    win, don = S.window(pb, donors, *S.LATE)
    donors = donors.loc[don].sort_values("day")
    comps, emb = sorted(win.compartment.unique()), list(donors.index)
    Y = np.full((len(emb), len(comps), len(GENES)), np.nan)
    L, wg = np.zeros((len(comps), 2, len(GENES))), np.zeros((len(comps), len(GENES)))
    for k, c in enumerate(comps):
        r = win[win.compartment == c]
        counts = r[GENES].to_numpy(float)
        wg[k] = gene_weights(counts)
        L[k] = projector(WF, wg[k])
        Y[[emb.index(d) for d in r.donor], k] = log_fraction(counts, r.depth.to_numpy(float))
    fr = Frame(win, donors, comps, donors.t_clock.to_numpy(float), donors.batch.astype(str).to_numpy(), Y, L, wg)
    ds = cs.Dataset(win, anchors=comps, donors=donors, name="a236").data  # the package's states, as a check
    mine = pd.DataFrame(fr.Z, index=emb, columns=comps).stack()
    theirs = pd.Series(ds.beta[:, 0] + 1j * ds.beta[:, 1], index=pd.MultiIndex.from_arrays(
        [[str(x).split("#")[0] for x in ds.resection], ds.compartment]))
    assert np.allclose(mine.reindex(theirs.index).to_numpy(), theirs.to_numpy()), "states differ from clocksync"
    return fr


# ---- timed model ------------------------------------------------------------------------------------------------

def design(t, batch, drift=True):
    """Columns: batch intercepts, drift (t - 9 h, optional), cos wt, sin wt."""
    B = (batch[:, None] == np.unique(batch)[None]).astype(float)
    return np.column_stack([B] + ([t - 9.0] if drift else []) + [np.cos(W * t), np.sin(W * t)])


def _patterns(present):
    """Groups of responses (columns) that share a presence pattern over embryos."""
    pats, inv = np.unique(present, axis=1, return_inverse=True)
    return [(pats[:, j], np.flatnonzero(np.ravel(inv) == j)) for j in range(pats.shape[1])]


def _usable(Dp):
    Dp = Dp[:, np.abs(Dp).sum(0) > 0]
    ok = len(Dp) > Dp.shape[1] and np.linalg.matrix_rank(Dp) == Dp.shape[1]
    return Dp, ok


def harmonics(Z, D):
    """OLS of each column of Z (embryo x response; NaN = absent) on D. Returns the cos and sin coefficients, the
    residual mean square, the (cos, sin) block of (D'D)^-1 and the residuals; NaN where the design is not of full
    rank with residual degrees of freedom."""
    m = Z.shape[1]
    c, s = np.full(m, np.nan, Z.dtype), np.full(m, np.nan, Z.dtype)
    ms, C, res = np.full(m, np.nan), np.full((m, 2, 2), np.nan), np.full(Z.shape, np.nan, Z.dtype)
    for pat, cols in _patterns(np.isfinite(Z)):
        Dp, ok = _usable(D[pat])
        if not ok:
            continue
        B = np.linalg.pinv(Dp) @ Z[np.ix_(pat, cols)]
        r = Z[np.ix_(pat, cols)] - Dp @ B
        c[cols], s[cols] = B[-2], B[-1]
        ms[cols] = (np.abs(r) ** 2).sum(0) / (len(Dp) - Dp.shape[1])
        C[cols] = np.linalg.inv(Dp.T @ Dp)[-2:, -2:]
        res[np.ix_(pat, cols)] = r
    return c, s, ms, C, res


# ---- H236.2: rotation -------------------------------------------------------------------------------------------

def polarisation(c, s, w):
    """rho = 2 sum w Im(s conj c) / sum w (|c|^2 + |s|^2), with its numerator and denominator."""
    ok = np.isfinite(c) & np.isfinite(w)
    num = 2 * np.sum(w[ok] * np.imag(s[ok] * np.conj(c[ok])))
    den = np.sum(w[ok] * (np.abs(c[ok]) ** 2 + np.abs(s[ok]) ** 2))
    return num / den, num, den


def rotation(Z, t, batch, drift=True) -> dict:
    """Polarisation of the time-locked component and its leave-one-batch-out jackknife t (weights 1 / residual mean
    square of the full fit, fixed over the jackknife)."""
    c, s, ms, C, _ = harmonics(Z, design(t, batch, drift))
    w = 1.0 / ms
    rho, num, den = polarisation(c, s, w)
    loo = np.array([polarisation(*harmonics(Z[batch != b], design(t[batch != b], batch[batch != b], drift))[:2], w)[0]
                    for b in np.unique(batch)])
    n = len(loo)
    se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    noise = np.nansum(w * ms * (C[:, 0, 0] + C[:, 1, 1]))  # expected noise part of the denominator
    return dict(n_celltypes=int(np.isfinite(c).sum()), rho=rho, se=se, t=rho / se, p=float(stats.t.sf(rho / se, n - 1)),
                rho_noise_corrected=num / (den - noise) if den > noise else np.nan,
                signal_share=1 - noise / den, loo_min=loo.min(), loo_max=loo.max())


# ---- H236.1: gene content ---------------------------------------------------------------------------------------

def gene_harmonics(Y, wg, t, batch, drift=True):
    """Common cos / sin coefficients (a_g, b_g) over cell types of the log fractions centred over genes, with
    cell-type-specific batch intercepts and drift, and the cell type's gene weights."""
    Yc = Y - Y.mean(2, keepdims=True)
    D = design(t, batch, drift)
    G = Y.shape[2]
    A, m = np.zeros((G, 2, 2)), np.zeros((G, 2))
    for pat, cols in _patterns(np.isfinite(Yc[:, :, 0])):
        Dp, ok = _usable(D[pat])
        if not ok:
            continue
        N, H = Dp[:, :-2], Dp[:, -2:]
        Hr = H - N @ (np.linalg.pinv(N) @ H)  # harmonics residualised on the nuisance columns
        A += np.einsum("kg,hj->ghj", wg[cols], Hr.T @ Hr)
        m += np.einsum("kg,nh,nkg->gh", wg[cols], Hr, Yc[np.ix_(pat, cols)])
    ab = np.linalg.solve(A, m[..., None])[..., 0]
    return ab[:, 0], ab[:, 1]


def null_planes(wf, n, rng):
    """Clock-plane projectors of `n` waveforms with acrophases permuted over genes."""
    return np.stack([plane(wf.permuted(rng)) for _ in range(n)])


def share(a, b, P):
    """Share of the centred harmonics (a, b) in the plane(s) P (genes x genes, or a stack)."""
    a, b = a - a.mean(), b - b.mean()
    return (np.einsum("g,...gh,h->...", a, P, a) + np.einsum("g,...gh,h->...", b, P, b)) / (a @ a + b @ b)


def plane_test(a, b, wf, planes) -> dict:
    """S1 against acrophase-permuted planes (the registered null; anti-conservative, see the power note)."""
    s1, null = float(share(a, b, plane(wf))), share(a, b, planes)
    return dict(S1=s1, acro_null_mean=null.mean(), acro_null_q95=np.quantile(null, 0.95),
                p_acro=(1 + np.sum(null >= s1 - 1e-12)) / (1 + len(null)))


def blocks(Y, t, batch, drift=True):
    """Gene-level timed model on cell types present in all embryos: the fitted nuisance part (intercepts, drift) and
    the residual blocks (embryo x cell type x gene), rescaled for the fitted degrees of freedom."""
    D = design(t, batch, drift)
    B = np.einsum("pn,nkg->pkg", np.linalg.pinv(D), Y)
    E = (Y - np.einsum("np,pkg->nkg", D, B)) * np.sqrt(len(D) / (len(D) - D.shape[1]))
    return np.einsum("np,pkg->nkg", D[:, :-2], B[:-2]), E


def draw(nuisance, E, rng, signal=0.0):
    """One plasmode data set: residual blocks reassigned between embryos (random permutation and sign)."""
    return nuisance + signal + E[rng.permutation(len(E))] * rng.choice([-1.0, 1.0], len(E))[:, None, None]


def null_S1(Y, wg, t, batch, wf, drift, n, rng):
    """S1 of `n` plasmode data sets without a time-locked signal (the null of H236.1 after the power note)."""
    nuisance, E = blocks(Y, t, batch, drift)
    P = plane(wf)
    return np.array([share(*gene_harmonics(draw(nuisance, E, rng), wg, t, batch, drift), P) for _ in range(n)])


def content(fr, drop=(), drift=True, n_null=10_000, seed=0, complete=True) -> tuple[dict, pd.DataFrame]:
    """H236.1 on the waveform without `drop`: the test row and the per-gene table. `complete`: cell types present in
    all embryos, the set on which the plasmode null is defined (p_sim); otherwise S1 only."""
    wf = WF.drop(list(drop)) if drop else WF
    g = [GENES.index(x) for x in wf.genes]
    k = np.flatnonzero(np.isfinite(fr.Y[:, :, 0]).all(0)) if complete else np.arange(len(fr.comps))
    Y, wg = fr.Y[:, k][:, :, g], fr.wg[k][:, g]
    a, b = gene_harmonics(Y, wg, fr.t, fr.batch, drift)
    rng = np.random.default_rng([SEED, 1, seed])
    row = dict(n_genes=len(g), n_celltypes=len(k), **plane_test(a, b, wf, null_planes(wf, 10_000, rng)))
    if complete and n_null:
        null = null_S1(Y, wg, fr.t, fr.batch, wf, drift, n_null, rng)
        row.update(sim_null_mean=null.mean(), sim_null_q95=np.quantile(null, 0.95),
                   p_sim=(1 + np.sum(null >= row["S1"] - 1e-12)) / (1 + n_null))
    tab = pd.DataFrame(dict(gene=wf.genes, a=a, b=b, amplitude=np.hypot(a, b),
                            peak_h=np.mod(np.arctan2(b, a) * HOURS, 24), ref_amp=wf.amp,
                            ref_psi_h=np.mod(wf.psi * HOURS, 24)))
    tab["peak_minus_ref_h"] = np.mod(tab.peak_h - tab.ref_psi_h, 24)
    return row, tab


def profiles(fr) -> pd.DataFrame:
    """Descriptive (post hoc): per gene, the mean over batches and complete cell types (gene weights) of the
    gene-centred log fraction at each clock time, centred within batch and cell type."""
    k = np.isfinite(fr.Y[:, :, 0]).all(0)
    Yc = fr.Y[:, k] - fr.Y[:, k].mean(2, keepdims=True)
    for b in np.unique(fr.batch):
        Yc[fr.batch == b] -= Yc[fr.batch == b].mean(0)
    w = fr.wg[k] / fr.wg[k].sum(0)
    return pd.DataFrame({f"{h:g} h": np.einsum("nkg,kg->g", Yc[fr.t == h], w) / np.sum(fr.t == h)
                         for h in np.unique(fr.t)}, index=pd.Index(GENES, name="gene")).reset_index()


def arms(tab) -> dict:
    """Peak time (h) of each arm's summed harmonic and their difference (12 h: antiphase)."""
    h = (tab.a + 1j * tab.b).set_axis(tab.gene)
    peak = {k: np.mod(np.angle(h[v].sum()) * HOURS, 24) for k, v in ARMS.items()}
    ref = pd.Series(WF.amp * np.exp(1j * WF.psi), index=GENES)
    refd = np.mod(np.angle(ref[ARMS["BMAL1 arm"]].sum() / ref[ARMS["E-box / D-box arm"]].sum()) * HOURS, 24)
    return {**{f"peak_h {k}": v for k, v in peak.items()},
            "difference_h": np.mod(peak["BMAL1 arm"] - peak["E-box / D-box arm"], 24), "reference_difference_h": refd}


def timed_R(win, donors, comps, drop=None, n_perm=10_000, seed=SEED) -> dict:
    """H108.1's statistic (host phase against clock time), optionally without some genes."""
    ds = cs.Dataset(win, anchors=comps, donors=donors, name="a236", drop_genes=drop)
    r = ds.timed("t_clock", n_perm=n_perm, seed=seed)
    return dict(R=r["R"], p_R=r["p_R"], offset_h=r["offset_h"])


# ---- power: gene-level plasmode ---------------------------------------------------------------------------------

LAMBDA = (0.25, 0.5, 1.0, 2.0)  # planted signal-to-noise ratios of the 24-h power, beside the observed one
KINDS = {"axis": (0.5, 0.5), "rotation rho 1": (1.0, 0.0), "rotation rho 0.5": (0.75, 0.25)}  # shares of a+, a-


def planted(noise, lam, kind):
    """Cos and sin coefficients (per cell type) of a planted clock-plane signal common in orientation to the cell
    types, with 24-h power lam x the expected noise power `noise` of the cell type's fitted harmonics."""
    fp, fm = KINDS[kind]
    P = lam * noise / 2
    ap, am = np.sqrt(fp * P), np.sqrt(fm * P)
    return (ap + am).astype(complex), 1j * (ap - am)


def power(fr, n_sim, n_null=10_000) -> tuple[pd.DataFrame, dict]:
    """Size and power of H236.1 and H236.2 on the cell types present in all 16 embryos (`blocks`, `draw`). Signals
    are planted in the clock plane at a signal-to-noise ratio lam of the weighted 24-h power (the observed ratio and
    a grid), or on an acrophase-permuted design ("non-clock harmonics"). Prints no observed rho, numerator or S1."""
    rng = np.random.default_rng([SEED, 2])
    full = np.flatnonzero(np.isfinite(fr.Y[:, :, 0]).all(0))
    Y, L, wg, t, batch = fr.Y[:, full], fr.L[full], fr.wg[full], fr.t, fr.batch
    D = design(t, batch)
    nuisance, E = blocks(Y, t, batch)
    state = lambda y: (lambda b: b[..., 0] + 1j * b[..., 1])(np.einsum("kcg,skg->skc", L, y))
    c, s, ms, C, _ = harmonics(state(Y), D)
    noise = ms * (C[:, 0, 0] + C[:, 1, 1])  # expected noise part of |c|^2 + |s|^2
    snr = float(np.sum((np.abs(c) ** 2 + np.abs(s) ** 2) / ms) / np.sum(noise / ms) - 1)
    cos, sin = np.cos(W * t)[:, None], np.sin(W * t)[:, None]

    def signal(kind, lam, X):
        pc, ps = planted(noise, lam, kind)
        z = cos * pc + sin * ps  # embryo x cell type
        return np.einsum("gc,skc->skg", X, np.stack([z.real, z.imag], -1))

    # code check on noiseless signals: rho = 1, 0.5, 0; S1 = 1 for a state common to the cell types. With planted
    # amplitudes that differ between cell types the gene-specific weights pull S1 below 1 (its ceiling, recorded).
    for kind, want in (("rotation rho 1", 1.0), ("rotation rho 0.5", 0.5), ("axis", 0.0)):
        y0 = nuisance + signal(kind, 1.0, WF.X)
        assert abs(polarisation(*harmonics(state(y0), D)[:2], 1 / ms)[0] - want) < 1e-6, kind
    common = np.einsum("gc,sc->sg", WF.X, np.column_stack([np.cos(W * t), np.sin(W * t)]))[:, None, :]
    assert abs(share(*gene_harmonics(nuisance + common, wg, t, batch), plane(WF)) - 1) < 1e-9
    ceiling = float(share(*gene_harmonics(nuisance + signal("axis", 1.0, WF.X), wg, t, batch), plane(WF)))

    planes = null_planes(WF, 2000, rng)
    null = null_S1(Y, wg, t, batch, WF, True, n_null, rng)
    s1_crit = float(np.quantile(null, 0.95))
    tcrit = stats.t.ppf(0.95, len(np.unique(batch)) - 1)
    lams = sorted(set(LAMBDA) | ({round(snr, 3)} if snr > 0.02 else set()))
    scen = [("no signal", None, 0.0)] + [(k, k, lam) for lam in lams for k in KINDS] + \
           [("non-clock harmonics", "axis", lam) for lam in lams]
    draws = {}
    for name, kind, lam in scen:
        sig, out = 0.0 if kind is None else signal(kind, lam, WF.X), []
        for _ in range(n_sim):
            if name == "non-clock harmonics":
                sig = signal(kind, lam, WF.permuted(rng).X)
            y = draw(nuisance, E, rng, sig)
            r = rotation(state(y), t, batch)
            g = plane_test(*gene_harmonics(y, wg, t, batch), WF, planes)
            out.append((r["rho"], r["t"], g["S1"], g["p_acro"]))
        draws[name, lam] = np.array(out)
        print(f"  {name}, snr {lam}: done", flush=True)
    # H236.2: nulls are every scenario without a rotation; the critical value is the largest simulated 95 % quantile
    nulls = [k for k in draws if not k[0].startswith("rotation")]
    worst = max(nulls, key=lambda k: np.quantile(draws[k][:, 1], 0.95))
    t_sim = float(np.quantile(draws[worst][:, 1], 0.95))
    rows = [dict(scenario=name, snr=lam, observed_snr=bool(abs(lam - round(snr, 3)) < 1e-9), n_sim=n_sim,
                 rho_mean=o[:, 0].mean(), rho_sd=o[:, 0].std(), t_q95=np.quantile(o[:, 1], 0.95),
                 reject_H236_2_t3=np.mean(o[:, 1] > tcrit), reject_H236_2_sim=np.mean(o[:, 1] > t_sim),
                 S1_mean=o[:, 2].mean(), S1_q95=np.quantile(o[:, 2], 0.95),
                 reject_H236_1_sim=np.mean(o[:, 2] > s1_crit), reject_H236_1_acro=np.mean(o[:, 3] < 0.05))
            for (name, lam), o in draws.items()]
    np.save(OUT / "null_t.npy", draws[worst][:, 1])
    info = dict(n_celltypes=len(full), observed_snr=snr, t_crit_t3=tcrit, t_crit_sim=t_sim,
                t_crit_sim_scenario=f"{worst[0]}, snr {worst[1]}", S1_crit_sim=s1_crit,
                S1_null_mean=float(null.mean()), S1_ceiling_planted=ceiling, n_null=n_null)
    return pd.DataFrame(rows), info


# ---- H236.3: companions -----------------------------------------------------------------------------------------

G = {}


def genome(fr=None):
    """Genome-wide late-window pseudobulks restricted to the A108 rows: counts (rows x genes), the row table
    (donor, compartment, stratum, n_cells, depth) and gene symbols (clock genes by their waveform name, the others
    upper-cased; of duplicated symbols the most expressed is kept)."""
    if not G:
        import anndata as ad
        import scipy.sparse as sp
        a = ad.read_h5ad(GENOME)
        win = (fr or load()).win
        key = pd.MultiIndex.from_arrays([a.obs.donor.astype(str), a.obs.compartment.astype(str)])
        want = pd.MultiIndex.from_arrays([win.donor.astype(str), win.compartment.astype(str)])
        i = pd.Series(np.arange(a.n_obs), index=key).reindex(want)
        assert i.notna().all(), "late-window rows missing from the genome-wide table"
        X = a.X[i.to_numpy(int)]
        X = np.asarray(X.toarray() if sp.issparse(X) else X, float)
        panel = pd.read_csv(S.QIU / "prep" / "panel.csv")
        clock = {m: canonical(f) for f, m in zip(panel.feature_name, panel.mouse_symbol) if canonical(f) in GENES}
        sym = np.array([clock.get(x, str(x).upper()) for x in a.var["symbol"].astype(str)])
        order = np.argsort(-X.mean(0), kind="stable")
        keep = np.sort(order[~pd.Series(sym[order]).duplicated().to_numpy()])
        rows = win[["donor", "compartment", "stratum", "n_cells", "depth"]].reset_index(drop=True)
        assert np.array_equal(X[:, keep[[list(sym[keep]).index(g) for g in GENES]]], win[GENES].to_numpy(float)), \
            "clock-gene counts differ from the stored table"
        G.update(X=X[:, keep], sym=sym[keep], rows=rows)
    return G["X"], G["rows"], G["sym"]


def log_mean():
    """Log of the mean fraction over the late-window pseudobulks, genes with mean count >= 1 (as A189)."""
    X, rows, sym = genome()
    return pd.Series(np.log((X / rows.depth.to_numpy()[:, None]).mean(0) + 1e-12), index=sym)[X.mean(0) >= 1]


def pseudo_job(args):
    kind, j, genes, n_perm = args
    limit_threads(1)
    X, rows, sym = genome()
    donors = data.load("mouse_qiu")[1]
    col = pd.Series(np.arange(len(sym)), index=sym)
    tab = pd.concat([rows, pd.DataFrame(X[:, col[list(genes)].to_numpy()], columns=list(WF.genes))], axis=1)
    out = dict(kind=kind, set=j, genes=";".join(genes))
    ds = cs.Dataset(tab, anchors=sorted(tab.compartment.unique()), donors=donors, name=f"{kind}{j}")
    for tag, d in (("", ds), ("_resid", S.residualised(ds, donors.day))):
        r = d.timed("t_clock", n_perm=n_perm, seed=SEED + j)
        out.update({f"R{tag}": r["R"], f"p_R{tag}": r["p_R"]})
    return out


def companions(jobs, n_pseudo, n_perm) -> tuple[pd.DataFrame, pd.DataFrame]:
    """R of the clock genes (waveform order: the same code path as the pseudo sets) and of the pseudo sets."""
    lm = log_mean()
    sets = [("clock", 0, tuple(GENES))]
    for kind, tag in (("matched", 3), ("responsive", 4)):
        ws = pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, tag]), WF, kind)
        sets += [(kind, j, tuple(w.genes)) for j, w in enumerate(ws)]
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as pool:
        res, t0 = [], time.time()
        for r in pool.imap_unordered(pseudo_job, [(k, j, g, 10_000 if k == "clock" else n_perm)
                                                  for k, j, g in sets]):
            res.append(r)
            if len(res) % 20 == 0:
                print(f"  {len(res)} / {len(sets)} sets, {time.time() - t0:.0f} s", flush=True)
    t = pd.DataFrame(res).sort_values(["kind", "set"]).reset_index(drop=True)
    obs, rows = t[t.kind == "clock"].iloc[0], []
    for tag in ("", "_resid"):
        row = dict(statistic="R" + tag, R_clock=obs["R" + tag], p_clock=obs["p_R" + tag])
        for kind in ("matched", "responsive"):
            v = t.loc[t.kind == kind, "R" + tag]
            row.update({f"{kind}_n": len(v), f"{kind}_median": v.median(), f"{kind}_q95": v.quantile(0.95),
                        f"{kind}_pct_below_clock": 100 * np.mean(v < obs["R" + tag]),
                        f"{kind}_share_p05": np.mean(t.loc[t.kind == kind, "p_R" + tag] < 0.05)})
        row["clock_specific"] = bool(row["R_clock"] > max(row["matched_q95"], row["responsive_q95"]))
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["n_genes_matching"] = len(lm)
    summary["n_responsive_genes_used"] = len(set(";".join(t.genes[t.kind == "responsive"]).split(";")))
    return t, summary


# ---- stages -----------------------------------------------------------------------------------------------------

def score(fr):
    comps, traj = np.array(fr.comps), S.trajectory_of().reindex(fr.comps).to_numpy()
    # H236.1 and H236.1b
    r1, genes = content(fr)
    r1b, _ = content(fr, drop=GC_CLOCK, seed=1)
    nd, genes_nd = content(fr, drift=False, seed=2)
    tests = [dict(test="H236.1", **r1), dict(test="H236.1b without PER1, PER2", **r1b),
             dict(test="H236.1 sens: no drift term", **nd),
             dict(test="H236.1b sens: no drift term", **content(fr, drop=GC_CLOCK, drift=False, seed=3)[0]),
             dict(test="descriptive: all cell types", **content(fr, complete=False)[0]),
             dict(test="descriptive: all cell types, without PER1, PER2", **content(fr, drop=GC_CLOCK,
                                                                                    complete=False)[0])]
    full = timed_R(fr.win, fr.donors, fr.comps)
    timed = [dict(test="H108.1 (reproduced)", **full),
             dict(test="H108.1 without PER1, PER2", **timed_R(fr.win, fr.donors, fr.comps, GC_CLOCK))]
    logo = []
    for g in GENES:
        logo.append(dict(gene=g, S1_without=content(fr, drop=[g], n_null=0, seed=10)[0]["S1"],
                         R_without=timed_R(fr.win, fr.donors, fr.comps, [g], n_perm=1000)["R"]))
    genes = genes.merge(pd.DataFrame(logo), on="gene")
    # H236.2 (primary on the cell types of the plasmode; p_sim against its worst null, `power`)
    Z, null_t = fr.Z, np.load(OUT / "null_t.npy")
    full = np.isfinite(Z).all(0)
    p_sim = lambda r: dict(r, t_crit_sim=np.quantile(null_t, 0.95), p_sim=(1 + np.sum(null_t >= r["t"])) / (1 + len(null_t)))
    rot = [dict(test="H236.2", **p_sim(rotation(Z[:, full], fr.t, fr.batch))),
           dict(test="H236.2 sens: all cell types", **p_sim(rotation(Z, fr.t, fr.batch))),
           dict(test="H236.2 sens: no drift term", **p_sim(rotation(Z[:, full], fr.t, fr.batch, drift=False)))]
    for name in NAMED:
        rot.append(dict(test=f"descriptive: {name}", **rotation(Z[:, traj == name], fr.t, fr.batch)))
    for name, n in pd.Series(traj).value_counts().items():
        if n >= 5:
            rot.append(dict(test=f"descriptive: {name}", **rotation(Z[:, traj == name], fr.t, fr.batch)))
    t1, t2, t3 = pd.DataFrame(tests), pd.DataFrame(rot), pd.DataFrame(timed)
    t1.to_csv(OUT / "h236_1_content.csv", index=False)
    genes.to_csv(OUT / "h236_1_genes.csv", index=False)
    genes_nd.to_csv(OUT / "h236_1_genes_nodrift.csv", index=False)
    prof = profiles(fr)
    prof.to_csv(OUT / "h236_1_profiles.csv", index=False)
    t3.to_csv(OUT / "h236_1b_timed.csv", index=False)
    pd.DataFrame([arms(genes)]).to_csv(OUT / "h236_1_arms.csv", index=False)
    t2.to_csv(OUT / "h236_2_rotation.csv", index=False)
    for t in (t1, t3, genes, genes_nd, prof, pd.DataFrame([arms(genes)]), t2):
        print(t.round(4).to_string(index=False), "\n")


def primary() -> pd.DataFrame:
    """The three primary tests with Holm-adjusted p values, from the stage outputs. H236.3: the larger of the two
    one-sided ranks of the clock set's R among the matched and the responsive pseudo sets."""
    from clocksync.stats import holm
    c, r = pd.read_csv(OUT / "h236_1_content.csv"), pd.read_csv(OUT / "h236_2_rotation.csv")
    t, k = pd.read_csv(OUT / "h236_3_sets.csv"), pd.read_csv(OUT / "h236_3_companions.csv").iloc[0]
    obs = t.R[t.kind == "clock"].iloc[0]
    rank = {kind: (1 + np.sum(t.R[t.kind == kind] >= obs)) / (1 + np.sum(t.kind == kind))
            for kind in ("matched", "responsive")}
    tab = pd.DataFrame([
        dict(test="H236.1", statistic="S1", value=c.S1[0], critical=c.sim_null_q95[0], p=c.p_sim[0]),
        dict(test="H236.2", statistic="jackknife t of rho", value=r.t[0], critical=r.t_crit_sim[0], p=r.p_sim[0]),
        dict(test="H236.3", statistic="R", value=obs, critical=max(k.matched_q95, k.responsive_q95),
             p=max(rank.values()))])
    tab["p_holm"] = holm(tab.p.to_numpy())
    tab["pass"] = (tab.value > tab.critical) & (tab.p_holm < 0.05)
    tab["p_matched"], tab["p_responsive"] = [np.nan, np.nan, rank["matched"]], [np.nan, np.nan, rank["responsive"]]
    tab["responsive_distinct_genes"] = len(set(";".join(t.genes[t.kind == "responsive"]).split(";")))
    return tab


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["power", "score", "companions", "primary"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--n-sim", type=int, default=2000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--check", action="store_true", help="companions on 4 sets per kind")
    a = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    limit_threads(4 if a.stage != "companions" else 1)
    pd.set_option("display.width", 250)
    fr = load()
    print(f"late window: {len(fr.donors)} embryos, {len(fr.comps)} cell types, "
          f"{int(np.isfinite(fr.Y[:, :, 0]).all(0).sum())} in all embryos", flush=True)
    if a.stage == "power":
        tab, info = power(fr, a.n_sim)
        tab.to_csv(OUT / "power.csv", index=False)
        pd.Series(info).to_csv(OUT / "power_info.csv", header=["value"])
        print(tab.round(4).to_string(index=False))
        print(pd.Series(info).to_string())
    elif a.stage == "score":
        score(fr)
    elif a.stage == "primary":
        tab = primary()
        tab.to_csv(OUT / "h236_primary.csv", index=False)
        print(tab.round(4).to_string(index=False))
    else:
        sfx = "_check" if a.check else ""
        t, summary = companions(a.jobs, 4 if a.check else a.n_pseudo, 200 if a.check else 1000)
        t.to_csv(OUT / f"h236_3_sets{sfx}.csv", index=False)
        summary.to_csv(OUT / f"h236_3_companions{sfx}.csv", index=False)
        print(summary.round(4).T.to_string())
    print(f"runtime {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
