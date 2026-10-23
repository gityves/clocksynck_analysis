"""Helpers of organs: GTEx cohorts, the timed per-tissue fit and the departures from recorded time, factor
analysis of the shared departure, block correlation (brain vs periphery) and the A33 / A42 / A103 tables.

Clock states come from `clocksync` (`Dataset.data`: the clock-plane state beta and its Poisson covariance V of every
(donor, tissue) row, one projector per tissue; `Dataset.fit` for the untimed host). Everything marked "local: not in
clocksync" implements a statistic the package has no public call for (listed in the folder README). GTEx time of
death is controlled access: only aggregates leave these functions' callers.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import gtex  # noqa: E402

W = 2 * np.pi / 24
REF = "Muscle - Skeletal"
COHORTS = {"fast": gtex.FAST, "ventilator": gtex.VENTILATOR}
CLOCK = tuple(cs.gtex().genes)
# the 15 reference tissues of the host in the ventilator / cue analyses (legacy `v2_gtex.TISSUES`)
V2_TISSUES = (
    "Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
    "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
    "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum", "Brain - Cortex", "Whole Blood",
)


# ---- data --------------------------------------------------------------------------------------------------------

def cohort(name: str, min_donors: int = 40, genes="clock", subsample: int | None = None, seed: int = 0):
    """(pseudobulk, donors) of a GTEx cohort: fast deaths (Hardy 1-2) or ventilator cases (Hardy 0), tissues with
    >= `min_donors` donors (01_gtex_organs `cohort_table`), one sample per donor and tissue, depth = counts summed
    over all genes. `subsample`: N random donors (checks), applied before the tissue filter's donor count."""
    pb, donors = gtex.load(hardy=COHORTS[name], genes=genes, depth="sum")
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)]
    n = pb.groupby("compartment").donor.nunique()
    pb = pb[pb.compartment.isin(n.index[n >= (min(min_donors, 10) if subsample else min_donors)])]
    return pb.reset_index(drop=True), donors


def sample_ischaemic_h(pb: pd.DataFrame) -> np.ndarray:
    """Ischaemic time (h) of the sample behind every (donor, tissue) row of `pb` (local: `gtex.load` keeps one value
    per donor; the departures are regressed on the time of the tissue's own sample)."""
    import anndata as ad

    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["subject", "tissue", "ischemic_time"]].copy()
    finally:
        a.file.close()
    o["donor"], o["compartment"] = o.subject.astype(str), o.tissue.astype(str)
    o = o.drop_duplicates(["donor", "compartment"]).set_index(["donor", "compartment"]).ischemic_time.astype(float) / 60
    return o.reindex(pd.MultiIndex.from_frame(pb[["donor", "compartment"]])).to_numpy()


def dataset(pb: pd.DataFrame, waveform: cs.Waveform | None = None, anchors=None) -> cs.Dataset:
    """Dataset of a GTEx cohort (every tissue an anchor unless `anchors`); its rows keep the order of `pb`, so
    `ds.data.beta`, `ds.data.V` (clock-plane state and Poisson covariance, one projector per tissue) align with `pb`.
    `waveform`: a sub-waveform (e.g. without PER1 and PER2)."""
    wf = waveform if waveform is not None else cs.gtex()
    ds = cs.Dataset(pb[["donor", "compartment", "depth", *wf.genes]],
                    anchors=sorted(pb.compartment.unique()) if anchors is None else list(anchors), waveform=wf)
    assert (ds.table.resection.to_numpy() == pb.donor.to_numpy()).all()
    return ds


# ---- timed fit and departures from recorded time (local: not in clocksync) -------------------------------------------

def _psi_step(e2, v, psi, floor=1e-8, n_iter=100, tol=1e-10) -> float:
    for _ in range(n_iter):
        w2 = 1.0 / (psi + v) ** 2
        new = max(float(np.sum(w2 * (e2 - v)) / np.sum(w2)), floor)
        if abs(new - psi) <= tol * max(psi, floor):
            return new
        psi = new
    return psi


def timed_fit(beta, V, t_h, tissue, fit_rows=None, n_iter: int = 200, tol: float = 1e-10):
    """local: not in clocksync. Per tissue beta = mu_k + a_k e^{i w t} + e by Poisson-weighted least squares with the
    own-noise variance psi_k by maximum likelihood (the A36 / 01_gtex_organs timed fit). Rows `fit_rows` (default
    all) enter the fit; returns (labels, mu (K x 2), a (complex K), psi, complex residual of every row)."""
    beta, V, t_h = np.asarray(beta, float), np.asarray(V, float), np.asarray(t_h, float)
    fit_rows = np.ones(len(beta), bool) if fit_rows is None else np.asarray(fit_rows)
    lab = np.unique(tissue[fit_rows])
    f = np.column_stack([np.cos(t_h * W), np.sin(t_h * W)])
    mu, a, psi = np.zeros((len(lab), 2)), np.zeros(len(lab), complex), np.zeros(len(lab))
    I2 = np.eye(2)
    J = np.array([[0.0, -1.0], [1.0, 0.0]])
    for kk, k in enumerate(lab):
        i = (tissue == k) & fit_rows
        b, Vi, fi = beta[i], V[i], f[i]
        G = np.zeros((i.sum(), 2, 4))
        G[:, :, :2] = I2
        G[:, :, 2], G[:, :, 3] = fi, fi @ J.T
        v_eig, U = np.linalg.eigh(Vi)
        p = max(float(np.var(b, 0).mean()), 1e-8)
        for _ in range(n_iter):
            Wt = np.linalg.inv(p * I2 + Vi)
            GtW = np.swapaxes(G, 1, 2) @ Wt
            th = np.linalg.solve((GtW @ G).sum(0), np.einsum("nab,nb->a", GtW, b))
            e = b - np.einsum("nab,b->na", G, th)
            new = _psi_step(np.einsum("nai,na,nb,nbi->ni", U, e, e, U).ravel(), v_eig.ravel(), p)
            done = abs(new - p) <= tol * max(p, 1e-8)
            p = new
            if done:
                break
        mu[kk], psi[kk], a[kk] = th[:2], p, th[2] + 1j * th[3]
    k = np.searchsorted(lab, tissue)
    e = (beta[:, 0] + 1j * beta[:, 1]) - (mu[k, 0] + 1j * mu[k, 1]) - a[k] * np.exp(1j * t_h * W)
    return lab, mu, a, psi, e


def departures(e, a_rows, t_h, exact: bool = False):
    """Phase departure (h; > 0: the tissue's clock is ahead of recorded time) and log-amplitude departure of every row.
    w = (beta - mu) / (a e^{it}) = 1 + u: linear (default) Im u / omega and Re u; `exact`: arg(w) / omega, log|w|."""
    u = e / (a_rows * np.exp(1j * t_h * W))
    if exact:
        return np.angle(1 + u) / W, np.log(np.abs(1 + u))
    return u.imag / W, u.real


def wide(values, donor, tissue) -> pd.DataFrame:
    return pd.DataFrame({"d": donor, "k": tissue, "v": values}).pivot(index="d", columns="k", values="v")


# ---- factor analysis with missing entries (local: not in clocksync) --------------------------------------------------

def _posterior(Y0, M, B, psi):
    P = np.einsum("nj,jq,jp->nqp", M / psi, B, B)
    S = np.linalg.inv(np.eye(B.shape[1]) + P)
    m = np.einsum("nqp,np->nq", S, (Y0 * M / psi) @ B)
    return m, S


def loglik(Y, B, psi) -> float:
    M = ~np.isnan(Y)
    ll = 0.0
    for i in range(len(Y)):
        o = M[i]
        C = B[o] @ B[o].T + np.diag(psi[o])
        y = Y[i, o]
        ll -= 0.5 * (np.linalg.slogdet(2 * np.pi * C)[1] + y @ np.linalg.solve(C, y))
    return ll


def fa(Y: np.ndarray, q: int, n_iter: int = 2000, tol: float = 1e-7):
    """Y = B g + e, g ~ N(0, I_q), e ~ N(0, diag psi); NaN entries missing, columns centred by the caller. Returns B
    (p x q), psi (p), posterior means m (n x q) and covariances S (n x q x q)."""
    M = ~np.isnan(Y)
    Y0 = np.where(M, Y, 0.0)
    nk = M.sum(0)
    psi = (Y0 ** 2).sum(0) / nk
    if q == 0:
        return np.zeros((Y.shape[1], 0)), psi, np.zeros((len(Y), 0)), np.zeros((len(Y), 0, 0))
    C = np.nan_to_num(pd.DataFrame(Y).cov(min_periods=10).to_numpy())
    w, U = np.linalg.eigh(C)
    B = U[:, ::-1][:, :q] * np.sqrt(np.maximum(w[::-1][:q], 1e-3 * psi.mean()))
    psi = np.maximum(psi - (B ** 2).sum(1), 0.05 * psi)
    ll_old = -np.inf
    for it in range(n_iter):
        m, S = _posterior(Y0, M, B, psi)
        mm = S + np.einsum("nq,np->nqp", m, m)
        for j in range(Y.shape[1]):
            r = M[:, j]
            B[j] = np.linalg.solve(mm[r].sum(0), (Y0[r, j][:, None] * m[r]).sum(0))
            fit = m[r] @ B[j]
            psi[j] = max(np.mean(Y0[r, j] ** 2 - 2 * Y0[r, j] * fit + np.einsum("q,nqp,p->n", B[j], mm[r], B[j])),
                         1e-6 * psi.mean())
        if it % 5 == 4:  # the likelihood is checked every 5 EM steps (the legacy checked every step)
            ll = loglik(Y, B, psi)
            if abs(ll - ll_old) < tol * abs(ll):
                break
            ll_old = ll
    m, S = _posterior(Y0, M, B, psi)
    return B, psi, m, S


def heldout(Y, Yh, B, psi) -> np.ndarray:
    """Predictive log-density of each held-out entry (Yh non-NaN) given the donor's training entries (Y)."""
    M = ~np.isnan(Y)
    m, S = _posterior(np.where(M, Y, 0.0), M, B, psi)
    i, j = np.nonzero(~np.isnan(Yh))
    mean = np.einsum("nq,nq->n", m[i], B[j])
    var = psi[j] + np.einsum("nq,nqp,np->n", B[j], S[i], B[j])
    return -0.5 * (np.log(2 * np.pi * var) + (Yh[i, j] - mean) ** 2 / var)


def rank_cv(beta, V, t_h, donor, tissue, qs=(0, 1, 2, 3, 4), n_folds=10, permute=False, seed=0,
            exact: bool = False) -> pd.DataFrame:
    """Held-out log-density per entry of the phase, amplitude and joint departures for q shared factors (10-fold
    cross-validation of donor x tissue entries; the timed fit is refitted on the training rows of every fold).
    `permute`: donors shuffled within each tissue (no sharing between tissues)."""
    rng = np.random.default_rng(seed)
    donor = donor.copy()
    if permute:
        for k in np.unique(tissue):
            i = np.flatnonzero(tissue == k)
            donor[i] = donor[rng.permutation(i)]
    fold = rng.integers(0, n_folds, len(beta))
    rows = []
    for f in range(n_folds):
        train = fold != f
        lab, _, a, _, e = timed_fit(beta, V, t_h, tissue, train)
        ph, am = departures(e, a[np.searchsorted(lab, tissue)], t_h, exact)
        for name, cols in (("phase", [ph]), ("amplitude", [am]), ("joint", [ph, am])):
            Ys, Hs = [], []
            for c in cols:
                Wd = wide(c, donor, tissue)
                Wt = wide(np.where(train, c, np.nan), donor, tissue).reindex_like(Wd)
                centre, scale = Wt.mean(), Wt.std()
                Ys.append(((Wt - centre) / scale).to_numpy())
                Hs.append(((Wd.where(Wt.isna()) - centre) / scale).to_numpy())
            Y, H = np.hstack(Ys), np.hstack(Hs)
            for q in qs:
                B, psi, _, _ = fa(Y, q)
                rows.append(dict(analysis=name, q=q, fold=f, ll=heldout(Y, H, B, psi).sum(), n=int((~np.isnan(H)).sum())))
    out = pd.DataFrame(rows).groupby(["analysis", "q"])[["ll", "n"]].sum().reset_index()
    out["ll_per_entry"] = out.ll / out.n
    out["gain_vs_q0"] = out.ll_per_entry - out.groupby("analysis").ll_per_entry.transform("first")
    out["permuted"] = permute
    return out


def varimax(L: np.ndarray, n_iter: int = 200) -> np.ndarray:
    p, q = L.shape
    R = np.eye(q)
    for _ in range(n_iter):
        Lr = L @ R
        u, _, vt = np.linalg.svd(L.T @ (Lr ** 3 - Lr * (Lr ** 2).sum(0) / p))
        R = u @ vt
    return L @ R


# ---- brain vs periphery of the shared departure (local: not in clocksync) -----------------------------------------

def _score(X: pd.DataFrame) -> pd.Series:
    return X.mean(axis=1).where(X.notna().sum(axis=1) >= 2)


def reliability(X: pd.DataFrame, rng, n_split: int) -> float:
    """Split-half (tissues split at random) reliability of the block score, Spearman-Brown, mean over splits."""
    r = []
    for _ in range(n_split):
        c = rng.permutation(X.columns)
        h = pd.concat([_score(X[c[: len(c) // 2]]), _score(X[c[len(c) // 2:]])], axis=1).dropna()
        rr = h.corr().iloc[0, 1]
        r.append(2 * rr / (1 + rr))
    return float(np.mean(r))


def zwide(values, donor, tissue) -> pd.DataFrame:
    Wd = wide(values, donor, tissue)
    return (Wd - Wd.mean()) / Wd.std()


def block_stat(Z: pd.DataFrame, rng, n_split: int) -> dict:
    """Correlation over donors (>= 2 tissues in each block) between the mean standardised departure of the brain
    regions and of the other tissues, and its disattenuated value r* = r / sqrt(rel_periphery rel_brain)."""
    isb = Z.columns.str.startswith("Brain")
    both = pd.concat([_score(Z.loc[:, ~isb]), _score(Z.loc[:, isb])], axis=1).dropna()
    r = both.corr().iloc[0, 1]
    rp, rb = reliability(Z.loc[both.index, ~isb], rng, n_split), reliability(Z.loc[both.index, isb], rng, n_split)
    return dict(donors=len(both), r=r, rel_periph=rp, rel_brain=rb,
                r_star=r / np.sqrt(rp * rb) if rp > 0 and rb > 0 else np.nan)


def block_boot(Z: pd.DataFrame, rng, n_boot: int, n_split: int = 20) -> tuple[float, float]:
    """95 % donor-bootstrap interval of r* (reliabilities recomputed with `n_split` splits in each replicate)."""
    isb = Z.columns.str.startswith("Brain")
    el = Z.index[(Z.loc[:, isb].notna().sum(1) >= 2) & (Z.loc[:, ~isb].notna().sum(1) >= 2)]
    Ze = Z.loc[el]
    bs = [block_stat(Ze.iloc[rng.integers(0, len(Ze), len(Ze))].reset_index(drop=True), rng, n_split)["r_star"]
          for _ in range(n_boot)]
    bs = np.where(np.isfinite(bs), bs, np.inf)
    return tuple(np.quantile(bs, [0.025, 0.975])) if n_boot else (np.nan, np.nan)
