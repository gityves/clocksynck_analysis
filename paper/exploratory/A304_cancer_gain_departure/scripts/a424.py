"""A424: is the TCGA proliferation gain loss tumour purity? Per-pair gain contribution on proliferation and ABSOLUTE purity.

Registered in paper/docs/registered_experiments.md ("A424"; builds on A304). Data: A304's 517 TCGA tumour-normal pairs (10
cancers), projector, centring and Poisson correction as `lib_a304`. Purity: ABSOLUTE (`tcga_purity_absolute/mastercalls.txt`,
call status == called, tumour sample of the pair matched on the first 15 characters of its barcode); the A127 expression
proxy for all pairs when fewer than 80 % of pairs have it. Pairs without purity leave the regressions.

Per-pair contribution  u_i = Re(z_T,i conj(z_N,i) exp(-i arg C_c)) / V_N,c   (states centred per cancer, C_c and V_N,c over all
pairs of the cancer, V_N Poisson-corrected; the cancer mean of u is G up to a factor (n-1)/n).
Model: u ~ proliferation_z + purity_z (z within cancer), cancer intercepts, pair weights 1.
  H424.1  proliferation slope < 0 at fixed purity; Freedman-Lane (proliferation residualised on purity, residuals permuted within
          cancer; 10,000; seed [424, 1]); specificity: same slope on the 200 matched pseudo-clock sets of A304.
  H424.2  purity slope < 0 at fixed proliferation; roles exchanged (seed [424, 2]); Holm over the two.
Descriptive (score): (1) unadjusted proliferation slope; (2) leave-one-cancer-out of H424.1 and of A304 Q_G; (3) BRCA slope
within PAM50 strata; (4) per-cancer slopes, correlation of proliferation with purity, Cochran Q / I^2; (5) the model on CPTAC-3.

  power  simulation at the real design (A304 calibration of S, g, psi per cancer; observed proliferation and purity vectors,
         no clock slope): null, gain falling with proliferation (k_G across the within-cancer IQR), gain proportional to stroma
         share 1 - purity; rejection of H424.1 / H424.2 at 0.05 and 0.025, H424.1 without the purity adjustment (confound size);
         ABSOLUTE coverage and proliferation-purity correlation      -> output/a424/power/{power,design,coverage,corr}.csv
  score  refuses to run without power.csv      -> output/a424/score/{tests,unadjusted,per_cancer,heterogeneity,loo,pam50,cptac,
                                                  companions,coverage}.csv
    python a424.py --stage power|score [--jobs 4] [--reps 1000] [--perms 1000]   (score: --perms 10000 --sets 200)
    python a424.py --check            (2 cancers, simulated states only; prints no observed slope)
Long runs: RUN_CAP_GB=20 paper/common/run.sh a424.py --stage power
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib_a304 as L  # noqa: E402
import power_a304 as PW  # noqa: E402

OUT = L.OUT / "a424"
ABS = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tcga_purity_absolute/mastercalls.txt")
COVERAGE_RULE = 0.80
MIN_CANCER = 10      # pairs with purity for a per-cancer fit
SEEDS = {"H424.1": [424, 1], "H424.2": [424, 2]}
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


# ---- design -----------------------------------------------------------------------------------------------------------

def load_absolute() -> pd.Series:
    """ABSOLUTE purity by sample barcode (first 15 characters): called rows with a value; first aliquot if several."""
    t = pd.read_csv(ABS, sep="\t")
    t = t[(t["call status"] == "called") & t.purity.notna()]
    t = t.assign(key=t["sample"].str[:15]).sort_values("sample").drop_duplicates("key")
    return t.set_index("key").purity


def zwithin(x: np.ndarray, g: np.ndarray) -> np.ndarray:
    """z-score within groups g (NaN stays NaN)."""
    out = np.full(len(x), np.nan)
    for k in np.unique(g):
        i = (g == k) & np.isfinite(x)
        if i.sum() > 1:
            out[i] = (x[i] - x[i].mean()) / x[i].std(ddof=1)
    return out


class Design:
    """Pairs of the cancers kept, in cell order: cancer index, proliferation, purity (chosen source), ok mask, subtype."""

    def __init__(self, cells: list[dict], samples, cancers: list[str]):
        ab = load_absolute()
        rows = []
        for c in cells:
            case, _ = L.sh.pair_frames(samples, c["cancer"])
            case = case.loc[c["ann"].index]
            rows.append(pd.DataFrame(dict(cancer=c["cancer"], prolif=c["ann"].prolif.to_numpy(float),
                                          proxy=c["ann"].purity.to_numpy(float), subtype=c["ann"].subtype.to_numpy(),
                                          absolute=case.barcode.str[:15].map(ab).to_numpy(float)), index=c["ann"].index))
        t = pd.concat(rows)
        self.cov = t.groupby("cancer", sort=False).absolute.agg(n="size", n_absolute=lambda s: int(s.notna().sum()))
        self.coverage = float(t.absolute.notna().mean())
        self.source = "ABSOLUTE" if self.coverage >= COVERAGE_RULE else "proxy"
        t["purity"] = t.absolute if self.source == "ABSOLUTE" else t.proxy
        self.t = t.reset_index(drop=True)
        self.g = np.concatenate([np.full(len(r), i) for i, r in enumerate(rows)])
        self.names = [c["cancer"] for c in cells]
        self.ok = self.t.purity.notna().to_numpy() & np.isfinite(self.t.prolif.to_numpy())
        self.slices = np.cumsum([0] + [len(r) for r in rows])
        # keep cancers with enough pairs; z within cancer among the kept pairs
        keep = np.array([self.ok[self.g == i].sum() >= MIN_CANCER for i in range(len(cells))])
        self.ok &= keep[self.g]
        self.px = zwithin(np.where(self.ok, self.t.prolif, np.nan), self.g)
        self.pp = zwithin(np.where(self.ok, self.t.purity, np.nan), self.g)
        self.gk = self.g[self.ok]
        self.xk, self.pk = self.px[self.ok], self.pp[self.ok]

    def split(self, u_all: np.ndarray) -> list[np.ndarray]:
        return [u_all[a:b] for a, b in zip(self.slices[:-1], self.slices[1:])]


def contributions(zt, zn, vn) -> np.ndarray:
    """u_i of one cancer (pair order) from raw complex states and the normal's Poisson variances."""
    ztc, znc = zt - zt.mean(), zn - zn.mean()
    n = len(zn)
    VN = np.sum(np.abs(znc) ** 2) / (n - 1) - np.mean(vn)
    C = np.sum(ztc * np.conj(znc)) / (n - 1)
    return np.real(ztc * np.conj(znc) * np.exp(-1j * np.angle(C))) / VN


def contrib_all(cells: list[dict]) -> np.ndarray:
    return np.concatenate([contributions(c["zt"], c["zn"], c["vn"]) for c in cells])


# ---- regression -------------------------------------------------------------------------------------------------------

def demean(x: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Within-group centring along the last axis (cancer intercepts)."""
    x = np.asarray(x, float)
    out = x.copy()
    for k in np.unique(g):
        i = g == k
        out[..., i] = x[..., i] - x[..., i].mean(axis=-1, keepdims=True)
    return out


def perm_within(x: np.ndarray, g: np.ndarray, n: int, rng) -> np.ndarray:
    """n permutations of x within groups g (n x len(x))."""
    P = np.empty((n, len(x)))
    for k in np.unique(g):
        i = np.flatnonzero(g == k)
        P[:, i] = x[i][np.argsort(rng.random((n, len(i))), axis=1)]
    return P


def two_slopes(y, x, p):
    """OLS slopes of y on (x, p), no intercept (inputs centred); x may be (P, n). Returns slope of x and slope of p."""
    sxx, spp, sxp = (x * x).sum(-1), (p * p).sum(-1), (x * p).sum(-1)
    sxy, spy = (x * y).sum(-1), (p * y).sum(-1)
    det = sxx * spp - sxp ** 2
    return (spp * sxy - sxp * spy) / det, (sxx * spy - sxp * sxy) / det


def one_slope(y, x):
    return (x * y).sum(-1) / (x * x).sum(-1)


def fl_test(y, x, p, g, n_perm, rng, adjust=True):
    """Slope of x on y (given p) with Freedman-Lane permutation within g. Inputs are centred within g.
    Returns (slope, p one-sided less)."""
    if not adjust:
        obs = float(one_slope(y, x))
        null = one_slope(y, perm_within(x, g, n_perm, rng))
    else:
        obs = float(two_slopes(y, x, p)[0])
        xr = x - (x * p).sum() / (p * p).sum() * p
        null = two_slopes(y, perm_within(xr, g, n_perm, rng), p)[0]
    return obs, L.p_less(null, obs)


def models(u_all: np.ndarray, D: Design, n_perm: int, rngs: dict | None, adjust_h1: bool = True) -> dict:
    """Pooled slopes (and permutation p when rngs is given) of u on proliferation / purity within cancer."""
    y = demean(u_all[D.ok], D.gk)
    x, p = demean(D.xk, D.gk), demean(D.pk, D.gk)
    res = dict(slope_prolif=float(two_slopes(y, x, p)[0]), slope_purity=float(two_slopes(y, x, p)[1]),
               slope_prolif_unadj=float(one_slope(y, x)), slope_purity_unadj=float(one_slope(y, p)))
    if rngs is not None:
        res["p_prolif"] = fl_test(y, x, p, D.gk, n_perm, rngs["H424.1"])[1]
        res["p_purity"] = fl_test(y, p, x, D.gk, n_perm, rngs["H424.2"])[1]
        res["p_prolif_unadj"] = fl_test(y, x, p, D.gk, n_perm, rngs["unadj"], adjust=False)[1]
    return res


def holm(p: dict) -> dict:
    k = sorted(p, key=p.get)
    out, run = {}, 0.0
    for i, key in enumerate(k):
        run = max(run, min(1.0, (len(k) - i) * p[key]))
        out[key] = run
    return out


# ---- power stage ------------------------------------------------------------------------------------------------------

def gain_vector(D: Design, i: int, g0: float, scenario: str, k: float) -> np.ndarray:
    """Per-pair gain of cancer i: calibrated g, falling with proliferation by k across the IQR, or proportional to the
    stroma share 1 - purity (cancer mean fixed), or constant."""
    a, b = D.slices[i], D.slices[i + 1]
    n = b - a
    if scenario == "prolif":
        x = D.t.prolif.to_numpy(float)[a:b]
        q1, q3 = np.nanpercentile(x, [25, 75])
        return g0 * np.exp(np.log(k) * (x - np.nanmedian(x)) / (q3 - q1))
    if scenario == "purity":
        s = 1.0 - D.t.purity.to_numpy(float)[a:b] if D.source == "ABSOLUTE" else \
            np.clip(0.5 - 0.25 * D.pp_all[a:b], 0.05, None)      # proxy: z-scale mapped to a stroma share (see docs)
        s = np.where(np.isfinite(s), s, np.nanmean(s))
        return g0 * s / s.mean()
    return np.full(n, g0)


def one_rep(cells, D: Design, cal, scenario: str, k: float, n_perm: int, seed) -> dict:
    rng = np.random.default_rng(seed)
    us = []
    for i, cell in enumerate(cells):
        kc = cal.loc[cell["cancer"]]
        n = cell["n"]
        h = PW.cn(rng, kc.S, n)
        zn = h + PW.cn(rng, kc.psi_N2, n) + PW.cn(rng, cell["vn"], n)
        zt = gain_vector(D, i, kc.g, scenario, k) * h + PW.cn(rng, kc.psi_T2, n) + PW.cn(rng, cell["vt"], n)
        us.append(contributions(zt, zn, cell["vn"]))
    rngs = {"H424.1": rng, "H424.2": rng, "unadj": rng}
    return models(np.concatenate(us), D, n_perm, rngs)


SCEN = [("null", "null", 1.0)] + [(f"prolif x{k}", "prolif", k) for k in (0.85, 0.7, 0.5, 0.35)] + [("purity only", "purity", 1.0)]


def design_tables(D: Design):
    cov = D.cov.assign(coverage=D.cov.n_absolute / D.cov.n).reset_index()
    cov.loc[len(cov)] = ["all", int(D.cov.n.sum()), int(D.cov.n_absolute.sum()), D.coverage]
    cov["source_used"] = D.source
    rows = []
    for i, c in enumerate(D.names):
        m = (D.g == i) & D.ok
        r = np.corrcoef(D.t.prolif[m], D.t.purity[m])[0, 1] if m.sum() > 2 else np.nan
        mp = (D.g == i) & np.isfinite(D.t.proxy) & np.isfinite(D.t.absolute)
        rows.append(dict(cancer=c, n_pairs=int((D.g == i).sum()), n_in_model=int(m.sum()),
                         r_prolif_purity=r, r_prolif_absolute=float(np.corrcoef(D.t.prolif[mp], D.t.absolute[mp])[0, 1]),
                         r_prolif_proxy=float(np.corrcoef(D.t.prolif[(D.g == i)], D.t.proxy[(D.g == i)])[0, 1]),
                         r_absolute_proxy=float(np.corrcoef(D.t.absolute[mp], D.t.proxy[mp])[0, 1])))
    rows.append(dict(cancer="pooled within cancer", n_pairs=len(D.g), n_in_model=int(D.ok.sum()),
                     r_prolif_purity=float(np.corrcoef(D.xk, D.pk)[0, 1])))
    return cov, pd.DataFrame(rows)


def setup(samples, states, ann, cancers):
    cells = [c for c in L.design(states, ann, "H304.1") if c["cancer"] in cancers]
    D = Design(cells, samples, cancers)
    D.pp_all = zwithin(D.t.purity.to_numpy(float), D.g)
    return cells, D


def stage_power(a, cancers):
    out = OUT / ("power_check" if a.check else "power")
    out.mkdir(parents=True, exist_ok=True)
    samples = L.sh.load_tcga()
    states = L.q.tcga_states()
    ann = L.annotations(samples)
    cells, D = setup(samples, states, ann, cancers)
    cal = PW.calibrate({c: states[c] for c in cancers}, L.reliability())
    cov, corr = design_tables(D)
    cov.to_csv(out / "coverage.csv", index=False)
    corr.to_csv(out / "corr.csv", index=False)
    log(f"purity source {D.source}, ABSOLUTE coverage {D.coverage:.3f}")
    pd.DataFrame(dict(cancer=D.names, n=[c["n"] for c in cells],
                      n_in_model=[int((D.gk == i).sum()) for i in range(len(cells))])).to_csv(out / "design.csv", index=False)
    rows = []
    for si, (name, scen, k) in enumerate(SCEN):
        res = pd.DataFrame(Parallel(n_jobs=a.jobs)(
            delayed(one_rep)(cells, D, cal, scen, k, a.perms, [424, 3, si, r]) for r in range(a.reps)))
        row = dict(scenario=name, k=k, reps=len(res), mean_slope_prolif=res.slope_prolif.mean(),
                   mean_slope_purity=res.slope_purity.mean(), mean_slope_prolif_unadj=res.slope_prolif_unadj.mean(),
                   sd_slope_prolif=res.slope_prolif.std(), sd_slope_purity=res.slope_purity.std())
        for key, lab in (("p_prolif", "H424.1"), ("p_purity", "H424.2"), ("p_prolif_unadj", "H424.1_unadjusted")):
            row[f"rej_{lab}_05"] = float((res[key] < 0.05).mean())
            row[f"rej_{lab}_025"] = float((res[key] < 0.025).mean())
        rows.append(row)
        log(name, {k_: round(v, 3) for k_, v in row.items() if k_.startswith("rej")})
    P = pd.DataFrame(rows)
    rule = P[P.scenario == "prolif x0.5"].iloc[0]
    pur = P[P.scenario == "purity only"].iloc[0]
    P["H424.1_scored"] = bool(rule["rej_H424.1_025"] >= 0.8)
    P["H424.2_scored"] = bool(pur["rej_H424.2_025"] >= 0.8)
    P.to_csv(out / "power.csv", index=False)
    log("power\n", P.round(3).to_string())


# ---- score stage ------------------------------------------------------------------------------------------------------

def companion_u(c: str, ann_c, n_sets: int, samples) -> list[np.ndarray]:
    """u of cancer c on each of the A304 matched pseudo-clock sets (seed [304, 7, cancer index])."""
    from clocksync.config import limit_threads
    from clocksync.paired import cohort_states, pseudo_waveforms
    limit_threads(1)
    case, _ = L.sh.pair_frames(samples, c)
    ct, cn = L.sh.tcga_full_counts(c, case.index)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([304, 7, L.CANCERS.index(c)]), L.sh.WF, kind="matched")
    out = []
    for w in ws:
        g = list(w.genes)
        cell = L.cell_from_raw(c, cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w), ann_c, "H304.1")
        out.append(contributions(cell["zt"], cell["zn"], cell["vn"]))
    return out


def ols_se(y, X):
    """Slopes and SEs of y on columns of X (centred within cancer already)."""
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    s2 = (r @ r) / (len(y) - X.shape[1] - 1)
    return b, np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))


def per_cancer_table(u_all, D: Design) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for i, c in enumerate(D.names):
        m = (D.g == i) & D.ok
        if m.sum() < MIN_CANCER:
            continue
        y, x, p = (v - v.mean() for v in (u_all[m], D.px[m], D.pp[m]))
        b, se = ols_se(y, np.c_[x, p])
        b1, se1 = ols_se(y, x[:, None])
        rows.append(dict(cancer=c, n=int(m.sum()), slope_prolif=b[0], se_prolif=se[0], slope_purity=b[1], se_purity=se[1],
                         slope_prolif_unadj=b1[0], se_prolif_unadj=se1[0], r_prolif_purity=float(np.corrcoef(x, p)[0, 1])))
    T = pd.DataFrame(rows)
    het = []
    for key in ("prolif", "purity", "prolif_unadj"):
        w = 1 / T[f"se_{key}"] ** 2
        mu = (w * T[f"slope_{key}"]).sum() / w.sum()
        Q = float((w * (T[f"slope_{key}"] - mu) ** 2).sum())
        df = len(T) - 1
        from scipy.stats import chi2
        het.append(dict(slope=key, k=len(T), fixed_effect=mu, Q=Q, df=df, p=float(chi2.sf(Q, df)),
                        I2=max(0.0, (Q - df) / Q) if Q > 0 else 0.0))
    return T, pd.DataFrame(het)


def leave_one_out(u_all, D: Design, cells, n_perm) -> pd.DataFrame:
    rows = []
    vals = [L.cell_values(c, c["lab"][None, :]) for c in cells]
    for i, c in enumerate(D.names + ["none"]):
        keep = [j for j in range(len(D.names)) if D.names[j] != c]
        m = np.isin(D.g, keep) & D.ok
        gk = D.g[m]
        y, x, p = demean(u_all[m], gk), demean(D.px[m], gk), demean(D.pp[m], gk)
        rng = np.random.default_rng([424, 5, i])
        obs, pv = fl_test(y, x, p, gk, n_perm, rng)
        Q = L.pooled([cells[j] for j in keep], [vals[j] for j in keep])["Q_G"][0]
        rows.append(dict(left_out=c, slope_prolif=obs, p=pv, A304_Q_G=float(Q)))
    return pd.DataFrame(rows)


def pam50(u_all, D: Design, n_perm) -> pd.DataFrame:
    """BRCA proliferation slope with PAM50 stratum intercepts (strata >= 5 pairs), with and without the purity adjustment."""
    i = D.names.index("BRCA")
    m = (D.g == i) & D.ok & pd.notna(D.t.subtype).to_numpy()
    sub = D.t.subtype.to_numpy()
    cnt = pd.Series(sub[m]).value_counts()
    m &= np.isin(sub, cnt.index[cnt >= 5])
    st = pd.factorize(sub[m])[0]
    x, p = zwithin(D.t.prolif.to_numpy(float)[m], np.zeros(m.sum())), zwithin(D.t.purity.to_numpy(float)[m], np.zeros(m.sum()))
    y, x, p = demean(u_all[m], st), demean(x, st), demean(p, st)
    rng = np.random.default_rng([424, 3])
    rows = []
    for name, adj in (("PAM50 strata + purity", True), ("PAM50 strata", False)):
        obs, pv = fl_test(y, x, p, st, n_perm, rng, adjust=adj)
        rows.append(dict(model=name, n=int(m.sum()), strata="; ".join(f"{k} {v}" for k, v in pd.Series(sub[m]).value_counts().items()),
                         slope_prolif=obs, p=pv))
    return pd.DataFrame(rows)


def cptac(a) -> pd.DataFrame:
    """The model on the CPTAC-3 pairs of A304r (proliferation = G2M; purity proxy = minus the mean of the stromal and immune
    scores, z within type)."""
    import a304r as R
    K = R.K
    st, ann = R.data()
    cells = R.cells_of(st, ann)
    us, g, x, p, ok = [], [], [], [], []
    for i, c in enumerate(cells):
        cases, ct, _ = K.load_type(c["cancer"])
        ix = [cases.index(v) for v in c["ann"].index]
        s = (R.gene_score(ct, R.STROMA).to_numpy()[ix] + R.gene_score(ct, R.IMMUNE).to_numpy()[ix]) / 2
        us.append(contributions(c["zt"], c["zn"], c["vn"]))
        g.append(np.full(c["n"], i)); x.append(c["ann"].prolif.to_numpy(float)); p.append(-s)
    u, g, x, p = np.concatenate(us), np.concatenate(g), np.concatenate(x), np.concatenate(p)
    xz, pz = zwithin(x, g), zwithin(p, g)
    y, xd, pd_ = demean(u, g), demean(xz, g), demean(pz, g)
    rng = np.random.default_rng([424, 6])
    rows = []
    for name, adj in (("adjusted for purity proxy", True), ("unadjusted", False)):
        obs, pv = fl_test(y, xd, pd_, g, a.perms, rng, adjust=adj)
        rows.append(dict(model=name, types=len(cells), pairs=len(u), slope_prolif=obs, p=pv))
    obs, pv = fl_test(y, pd_, xd, g, a.perms, rng)
    rows.append(dict(model="purity proxy at fixed proliferation", types=len(cells), pairs=len(u), slope_prolif=obs, p=pv))
    return pd.DataFrame(rows)


def stage_score(a, cancers):
    out = OUT / ("score_check" if a.check else "score")
    out.mkdir(parents=True, exist_ok=True)
    pw = pd.read_csv(OUT / "power" / "power.csv")
    scored = {"H424.1": bool(pw["H424.1_scored"].iloc[0]), "H424.2": bool(pw["H424.2_scored"].iloc[0])}
    samples = L.sh.load_tcga()
    states = L.q.tcga_states()
    ann = L.annotations(samples)
    cells, D = setup(samples, states, ann, cancers)
    cov, corr = design_tables(D)
    cov.to_csv(out / "coverage.csv", index=False)
    corr.to_csv(out / "corr.csv", index=False)
    u_all = contrib_all(cells)
    rngs = {"H424.1": np.random.default_rng(SEEDS["H424.1"]), "H424.2": np.random.default_rng(SEEDS["H424.2"]),
            "unadj": np.random.default_rng([424, 4])}
    obs = models(u_all, D, a.perms, rngs)
    per_cancer = dict(zip(cancers, Parallel(n_jobs=a.jobs)(
        delayed(companion_u)(c, ann[c], a.sets, samples) for c in D.names)))
    sets = []
    for j in range(a.sets):
        u = np.concatenate([per_cancer[c][j] for c in D.names])
        sets.append(models(u, D, 0, None))
    comp = pd.DataFrame(sets)
    comp.insert(0, "set", range(len(comp)))
    comp.to_csv(out / "companions.csv", index=False)
    rows = []
    for h, key, pkey in (("H424.1", "slope_prolif", "p_prolif"), ("H424.2", "slope_purity", "p_purity")):
        s = comp[key].dropna()
        rows.append(dict(hypothesis=h, slope=obs[key], p=obs[pkey], matched_median=float(s.median()),
                         matched_q05=float(s.quantile(.05)), matched_q95=float(s.quantile(.95)),
                         matched_pct=float((s < obs[key]).mean() * 100), matched_n=len(s), purity_source=D.source,
                         pairs=int(D.ok.sum()), power_scored=scored[h]))
    T = pd.DataFrame(rows)
    ph = holm(dict(zip(T.hypothesis, T.p)))
    T["p_holm"] = T.hypothesis.map(ph)
    T["clock_specific"] = T.slope < T.matched_q05
    T["detected"] = (T.p_holm < 0.05) & T.clock_specific & T.power_scored
    T.to_csv(out / "tests.csv", index=False)
    pd.DataFrame([dict(model="unadjusted proliferation", slope=obs["slope_prolif_unadj"], p=obs["p_prolif_unadj"],
                       n=int(D.ok.sum()))]).to_csv(out / "unadjusted.csv", index=False)
    pcan, het = per_cancer_table(u_all, D)
    pcan.to_csv(out / "per_cancer.csv", index=False)
    het.to_csv(out / "heterogeneity.csv", index=False)
    leave_one_out(u_all, D, cells, a.perms).to_csv(out / "loo.csv", index=False)
    pam50(u_all, D, a.perms).to_csv(out / "pam50.csv", index=False)
    cptac(a).to_csv(out / "cptac.csv", index=False)
    log("tests\n", T.round(4).to_string())


# ---- check ------------------------------------------------------------------------------------------------------------

def check(a):
    """Code check on 2 cancers: power path on 2 replicates, and every score-stage function on simulated contributions.
    Prints counts and finiteness only."""
    cancers = ["BRCA", "THCA"]
    a.reps, a.perms, a.jobs = 4, 50, 1
    stage_power(a, cancers)
    samples = L.sh.load_tcga()
    states = L.q.tcga_states()
    ann = L.annotations(samples)
    cells, D = setup(samples, states, ann, cancers)
    cal = PW.calibrate({c: states[c] for c in cancers}, L.reliability())
    rng = np.random.default_rng(0)
    us = []
    for i, cell in enumerate(cells):
        kc = cal.loc[cell["cancer"]]
        h = PW.cn(rng, kc.S, cell["n"])
        zn = h + PW.cn(rng, kc.psi_N2, cell["n"])
        us.append(contributions(kc.g * h + PW.cn(rng, kc.psi_T2, cell["n"]), zn, cell["vn"]))
    u = np.concatenate(us)
    rngs = {"H424.1": np.random.default_rng(1), "H424.2": np.random.default_rng(2), "unadj": np.random.default_rng(3)}
    res = models(u, D, 50, rngs)
    T, het = per_cancer_table(u, D)
    loo = leave_one_out(u, D, cells, 50)
    pm = pam50(u, D, 50)
    comp = companion_u("THCA", ann["THCA"], 1, samples)
    finite = all(np.isfinite(list(res.values()))) and np.isfinite(T.select_dtypes(float).to_numpy()).all() \
        and np.isfinite(het.select_dtypes(float).to_numpy()).all() and np.isfinite(loo.slope_prolif).all() \
        and np.isfinite(pm.slope_prolif).all() and np.isfinite(comp[0]).all() and len(comp[0]) == len(ann["THCA"])
    log(f"check: {len(cells)} cancers, {int(D.ok.sum())} pairs in model, score-stage functions finite on simulated u: {finite}")
    assert finite


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--perms", type=int, default=None)
    ap.add_argument("--sets", type=int, default=200)
    a = ap.parse_args()
    from clocksync.config import limit_threads
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.check:
        check(a)
    elif a.stage == "power":
        a.perms = a.perms or 1000
        stage_power(a, L.CANCERS)
    elif a.stage == "score":
        a.perms = a.perms or 10_000
        stage_score(a, L.CANCERS)
    else:
        ap.error("--stage or --check required")
    log("done")


if __name__ == "__main__":
    main()
