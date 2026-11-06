"""A225 helpers: A97 coherence from a (pairwise-complete) covariance, batched over many arms / permutations.

The A97 coherence (clocksync.decompose.components) regresses gene g on the 2-d clock state estimated from the other genes,
y_g ~ (y_-g X_-g (X_-g' X_-g)^-1) a; with gene-centred y this is a function of the covariance S only:
    a_g = (M_g' S_-g,-g M_g)^-1 M_g' S_-g,g,   M_g = X_-g (X_-g' X_-g)^-1,   phase_g = atan2(a_g[1], a_g[0]),
    coherence = |sum_g amp_g exp(i (phase_g - psi_g))| / sum_g amp_g.
`arm_cov` builds S for many row subsets at once with missing values (pairwise-complete: each gene pair uses the rows where both
are observed, means taken over those rows); with no missing value it equals the usual covariance, and `coherence_from_cov`
equals `components(...)['coherence']` (checked in `selfcheck`)."""
import numpy as np


def arm_cov(Y, M, W):
    """Pairwise-complete covariances of rows selected by each column of W.
    Y: n x G (values at missing entries ignored), M: n x G observed mask (0/1), W: n x P row indicators (0/1).
    Returns (P, G, G) covariance (NaN where fewer than 3 complete rows)."""
    n, G = Y.shape
    Y0 = np.where(M > 0, Y, 0.0); M = M.astype(float)
    ij = lambda A, B: (A[:, :, None] * B[:, None, :]).reshape(n, G * G)
    N = (W.T @ ij(M, M)).reshape(-1, G, G)
    A = (W.T @ ij(Y0, M)).reshape(-1, G, G)          # A[p,i,j] = sum over complete rows of y_i
    Sxy = (W.T @ ij(Y0, Y0)).reshape(-1, G, G)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = (Sxy - A * A.transpose(0, 2, 1) / N) / (N - 1)
    S[N < 3] = np.nan
    return S


def _design(wf):
    X, amp, psi = wf.X, np.asarray(wf.amp, float), np.asarray(wf.psi, float); G = len(amp); out = []
    for g in range(G):
        keep = np.arange(G) != g; Xk = X[keep]; out.append((keep, Xk @ np.linalg.inv(Xk.T @ Xk)))
    return out, amp, psi


def coherence_from_cov(S, wf, _cache={}):
    """A97 coherence for covariances S (..., G, G) of the genes of waveform `wf` (same order). NaN covariances give NaN."""
    key = id(wf)
    if key not in _cache:
        _cache[key] = (_design(wf), wf)
    (des, amp, psi), _ = _cache[key]
    S = np.asarray(S, float); lead = S.shape[:-2]; S = S.reshape((-1,) + S.shape[-2:]); ph = np.empty((len(S), len(amp)))
    for g, (keep, Mg) in enumerate(des):
        Skk = S[:, keep][:, :, keep]; Skg = S[:, keep, g]
        A = np.einsum("ik,pij,jl->pkl", Mg, Skk, Mg); b = np.einsum("ik,pi->pk", Mg, Skg)
        with np.errstate(all="ignore"):
            a = np.linalg.solve(A, b[..., None])[..., 0]
        ph[:, g] = np.arctan2(a[:, 1], a[:, 0])
    R = np.abs(np.sum(amp * np.exp(1j * (ph - psi)), axis=1)) / amp.sum()
    return R.reshape(lead)


def min_pair_n(M, idx, lab):
    """Smallest number of rows with both genes observed over all gene pairs, in the smaller-supported arm of a cell (arms lab == 1 and lab == 0 of rows `idx`)."""
    out = []
    for a in (1, 0):
        Ma = M[idx[lab == a]]; out.append((Ma.T @ Ma).min())
    return int(min(out))


def coh_masked(Y, M, idx, wf):
    """Coherence of rows `idx` of Y (n x G) with observed mask M (pairwise-complete)."""
    W = np.zeros((len(Y), 1)); W[idx, 0] = 1
    return float(coherence_from_cov(arm_cov(Y, M, W), wf)[0])


def selfcheck(Y, idx, wf, coh_ref):
    """Max abs difference between coherence_from_cov (no mask) and `coh_ref(Y, idx, wf)` (the A97 `components`)."""
    return abs(coh_masked(Y, np.ones_like(Y), idx, wf) - coh_ref(Y, idx, wf))


# ======================================================================================================================================
# semi-synthetic simulation (A225 part A)
# ======================================================================================================================================
LN2 = np.log(2)
EST = ("raw", "adjusted", "masked")


def centre_cols(X):
    return X - X.mean(0)


def leave_out_phases(Y, cw):
    """Per-tumour phase (rad) of the projection of the cohort-centred Y (n x G) on the complex waveform weights cw, gene j left out (column j)."""
    Yc = Y - Y.mean(0); tot = Yc @ cw
    return np.angle(tot[:, None] - Yc * cw)


def fit_cn(y, cn, ph):
    """OLS y ~ 1 + beta cn + cos ph + sin ph on rows with finite cn (A224 `fit_gene` model, one cohort). Returns (beta, residual SD, n), (0, SD of y, n) if cn has no variance."""
    ok = np.isfinite(cn) & np.isfinite(y) & np.isfinite(ph)
    if ok.sum() < 30 or cn[ok].std() == 0:                    # no CN call or no variation: not adjusted; residual SD = SD of y (phase terms not removed)
        return 0.0, float(np.nanstd(y, ddof=1)), int(ok.sum())
    y, cn, ph = y[ok], cn[ok], ph[ok]
    X = np.column_stack([cn, np.cos(ph), np.sin(ph)]); Xc = X - X.mean(0); yc = y - y.mean()
    b = np.linalg.lstsq(Xc, yc, rcond=None)[0]; res = yc - Xc @ b
    return float(b[0]), float(np.sqrt((res ** 2).sum() / (len(y) - 4))), int(ok.sum())


def adjust(Y, cn, ph):
    """Y - beta_g CN_g with beta_g fitted per gene on (Y, cn, ph) columns (ph: n x G phase covariate of each gene); missing CN counts as 0 in the offset."""
    b = np.array([fit_cn(Y[:, j], cn[:, j], ph[:, j])[0] for j in range(Y.shape[1])])
    return Y - b * np.nan_to_num(cn), b


def views(Y, cn, ph):
    """The three estimator inputs {name: (Y, observed mask)}: raw, dosage-adjusted (betas re-estimated), masked (GISTIC != 0 set missing)."""
    Ya, _ = adjust(Y, cn, ph)
    return {"raw": (Y, np.ones_like(Y)), "adjusted": (Ya, np.ones_like(Y)), "masked": (Y, (np.nan_to_num(cn) == 0).astype(float))}


def cell_stats(Y, M, cells, wf, nperm, rng, obs=True):
    """Per cell: observed D = coh(label 1) - coh(label 0) (if `obs`) and, if nperm, nperm label-permuted D (labels permuted within the cell).
    Returns (obs: (ncell,) or None, null: (nperm, ncell) or None)."""
    o, nul = [], []
    for c in cells:
        idx, lab = c["idx"], c["lab"]; Yc, Mc = Y[idx], M[idx]
        lp = rng.permuted(np.tile(lab, (nperm, 1)), axis=1).T if nperm else np.empty((len(lab), 0), int)
        a, b = ([lab == 1], [lab == 0]) if obs else ([], [])
        W = np.column_stack(a + [lp == 1] + b + [lp == 0]).astype(float)
        r = coherence_from_cov(arm_cov(Yc, Mc, W), wf); k = int(obs) + nperm
        if obs:
            o.append(r[0] - r[k])
        if nperm:
            nul.append(r[int(obs):k] - r[k + int(obs):])
    return (np.array(o) if obs else None), (np.column_stack(nul) if nperm else None)


def block_null(blocks, nperm, seed):
    """Null (nperm,) of the mean over all cells of D, labels permuted within cell; blocks: [(Y, M, cells, wf)]. NaN draws kept as NaN (counted by the caller)."""
    rng = np.random.default_rng(seed)
    nl = np.column_stack([cell_stats(Y, M, cells, wf, nperm, rng, obs=False)[1] for Y, M, cells, wf in blocks if cells])
    return np.nanmean(nl, axis=1)          # a cell whose permuted arm has a gene pair with < 3 complete rows is left out of that draw's mean


def block_obs(blocks):
    """Observed per-cell D of all blocks, concatenated."""
    return np.concatenate([cell_stats(Y, M, cells, wf, 0, None)[0] for Y, M, cells, wf in blocks if cells])


def gen_clock(par, lam, damp, rng):
    """Simulated ln expression (n x G) of the clock genes: internal phase uniform; m_g (0) + lam_i amp_g cos(phi - psi_g) + CN (beta + CN x phase terms) + Gaussian noise (residual SD)."""
    n = par["n"]; phi = rng.uniform(0, 2 * np.pi, n); lam_i = lam * np.where(par["is_other"], damp, 1.0)
    clock = lam_i[:, None] * par["amp"] * np.cos(phi[:, None] - par["psi"])
    dose = par["cn0"] * (par["beta"] + par["ic"] * np.cos(phi)[:, None] + par["is"] * np.sin(phi)[:, None])
    return clock + dose + rng.normal(size=(n, len(par["amp"]))) * par["sd"], phi


def gen_pseudo(ps, n, rng):
    """Pseudo genes: no clock term; own CN and beta, Gaussian noise with their own residual SD."""
    return ps["cn0"] * ps["beta"] + rng.normal(size=(n, ps["cn0"].shape[1])) * ps["sd"]


def luma_coh(Y, cells, wf):
    """Weighted (by LumA count) mean over cells of the raw coherence of the label-1 tumours (A219 `descriptive_per_group`)."""
    v = np.array([coherence_from_cov(arm_cov(Y[c["idx"]], np.ones((len(c["idx"]), Y.shape[1])), (c["lab"] == 1)[:, None].astype(float)), wf)[0] for c in cells])
    w = np.array([(c["lab"] == 1).sum() for c in cells], float); return float((v * w).sum() / w.sum())


def calibrate(par, cells, wf, target, seed, nrep=20, tol=0.003):
    """lam such that the mean raw LumA coherence over `nrep` simulations (common random numbers) equals `target` (bisection on [0.01, 5])."""
    def f(lam):
        return np.mean([luma_coh(gen_clock(par, lam, 1.0, np.random.default_rng([seed, i]))[0], cells, wf) for i in range(nrep)])
    lo, hi = 0.01, 5.0; flo, fhi = f(lo), f(hi)
    if not flo <= target <= fhi:
        return np.nan, flo, fhi
    for _ in range(40):
        mid = np.sqrt(lo * hi); fm = f(mid)
        if abs(fm - target) < tol:
            return mid, fm, None
        lo, hi = (mid, hi) if fm < target else (lo, mid)
    return mid, fm, None


def true_D(par, cells, lam, damp, wf, nrep, seed):
    """Planted coherence difference per cell: mean over `nrep` copy-number-free simulations (clock term + Gaussian noise at the real arm sizes, no CN terms, raw estimator)
    of coh(label 1) - coh(label 0). The population coherence of a rank-2 clock is ~1 at any amplitude (phase is consistent), so the planted damping only shows at finite n:
    the truth is therefore taken at the real cell sizes."""
    p0 = dict(par, beta=np.zeros_like(par["beta"]), ic=np.zeros_like(par["ic"]), **{"is": np.zeros_like(par["is"])})
    v = [cell_stats(gen_clock(p0, lam, damp, np.random.default_rng([seed, i]))[0], np.ones((par["n"], len(par["amp"]))), cells, wf, 0, None)[0] for i in range(nrep)]
    return np.mean(v, axis=0)


def replicate(r, ctx, lams, damp, nperm, seed):
    """One replicate, scenarios N (damp 1) and D (damp `damp`), all estimators, both cohorts. Rows: rep, scenario, estimator, level, D, p, pseudo_D."""
    from clocksync import stats as cstats
    rows = []
    for scen, dmp in (("N", 1.0), ("D", damp)):
        obs, nul, pobs = {e: {} for e in EST}, {e: {} for e in EST}, {e: {} for e in EST}
        for k, (name, c) in enumerate(ctx.items()):
            par = c["par"]; rg = np.random.default_rng([seed, r, k, 0])           # same clock draws for N and D
            Y, phi = gen_clock(par, lam=lams[name], damp=dmp, rng=rg)
            v = views(Y, par["cn"], leave_out_phases(Y, par["cw"]))
            rp = np.random.default_rng([seed, r, k, 1 + (scen == "D")])
            ps = c["pseudo"][r % len(c["pseudo"])]; Yp = gen_pseudo(ps, par["n"], np.random.default_rng([seed, r, k, 3]))
            ph14 = np.repeat(np.angle((Y - Y.mean(0)) @ par["cw"])[:, None], Yp.shape[1], 1)
            vp = views(Yp, ps["cn"], ph14)
            for e in EST:
                if c["cells_by"][e]:
                    obs[e][name], nul[e][name] = cell_stats(*v[e], c["cells_by"][e], par["wf"], nperm, rp)
                    pobs[e][name], _ = cell_stats(*vp[e], c["cells_by"][e], par["wf"], 0, rp)
        for e in EST:
            for lev in list(ctx) + ["pooled"]:
                nm = [x for x in (list(ctx) if lev == "pooled" else [lev]) if x in obs[e]]
                if not nm:
                    continue
                o = float(np.concatenate([obs[e][x] for x in nm]).mean()); nl = np.nanmean(np.column_stack([nul[e][x] for x in nm]), axis=1)
                rows.append(dict(rep=r, scenario=scen, estimator=e, level=lev, D=o, p=cstats.perm_p(nl, o, "greater", tie=0.0),
                                 pseudo_D=float(np.nanmean(np.concatenate([pobs[e][x] for x in nm])))))
    return rows


def replicate_chunk(rs, ctx, lams, damp, nperm, seed):
    return [row for r in rs for row in replicate(r, ctx, lams, damp, nperm, seed)]
