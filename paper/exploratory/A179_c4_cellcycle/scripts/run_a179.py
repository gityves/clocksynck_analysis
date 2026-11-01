"""A179 (C4, bulk): is the proliferation coherence loss carried by clock-gene variation that tracks the cell cycle?
Five proliferation contrasts (H127.1, H171.1, H171.2, H173.1, H173.3) with the registered strata, arms and centring (A176 definitions,
imported from A176's run_a176.py). Within each arm x stratum every clock gene is regressed on the G1/S (HALLMARK_E2F_TARGETS
without clock genes) and G2M (A97) scores; coherence is recomputed on the residuals. S = dR_conditioned - dR. Null: 200 random score pairs
(A168 expression-matched pseudo sets for the G1/S and the G2M gene lists, pool without clock-like, G2M, purity and E2F genes).
    python run_a179.py [--nrand 200] [--jobs 4]"""
import argparse
import sys
import time
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A176_c3_plane_variance/scripts"))
ap = argparse.ArgumentParser()
ap.add_argument("--nrand", type=int, default=200); ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
import run_a176 as D  # noqa: E402
D.args.jobs = args.jobs
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

L, A166, say, EXPR, components = D.L, D.A166, D.say, D.EXPR, D.components
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
SEED = zlib.crc32(b"A179") % (2 ** 31)
NAMES = ["H127.1", "H171.1", "H171.2", "H173.1", "H173.3"]
K = norm.ppf(0.95) + norm.ppf(0.80)  # 2.49
t0 = time.time()

R = {n: r for n, r in D.build_all().items() if n in NAMES}
for n, r in R.items():
    r["wf"] = L.wf_of(r["cg"]); r["Y"] = D.Yof(r, r["cg"])
reg = D.registered_dR()
for n in NAMES:  # first: reproduce the registered dR (1e-6)
    v = R[n]["c"].dR_registered(R[n]["Y"], R[n]["wf"]); assert abs(v - reg[n]) < 1e-6, (n, v, reg[n])
say("registered dR reproduced to 1e-6:", {n: round(reg[n], 6) for n in NAMES})


def score(X, genes, center=None):
    s = L.sh._zmean(X, genes)
    if center is not None:  # GSE68465: scores z-scored within site, as the registered G2M score
        s = s.groupby(center).transform(lambda v: (v - v.mean()) / v.std())
    return s.to_numpy(float)


def parts(c, Y, wf, S=None):
    """Per-stratum arm components; with S (n x 2 score matrix) every arm's Y is first residualised on [1, S].
    Returns mean dR, dlog_plane, dlog_off over strata and per-arm mean R^2 (n_genes, for case and reference arm)."""
    d, r2 = [], {c.case: [], c.ref: []}
    for i, l in zip(c.G, c.obs):
        cc = {}
        for v in (c.case, c.ref):
            j = i[l == v]; y = Y[j]; y = y - y.mean(0)
            if S is not None:
                C = np.column_stack([np.ones(len(j)), S[j]])
                res = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
                r2[v].append(1 - res.var(0) / y.var(0)); y = res
            cc[v] = components(y, np.zeros(y.shape[1]), 1, wf)
        d.append([cc[c.case][k] - cc[c.ref][k] for k in ("coherence", "log_plane", "log_off")])
    return np.mean(d, 0), {v: np.mean(x, 0) if x else None for v, x in r2.items()}


def cond_dR(r, S):
    return parts(r["c"], r["Y"], r["wf"], S)[0][0]


rows, desc = [], []
for ci, n in enumerate(NAMES):
    r = R[n]; X = r["X"]; cg = r["cg"]
    e2f = [g for g in L.hallmark("HALLMARK_E2F_TARGETS") if g in X.columns]
    g2m = [g for g in L.sh.hallmark_g2m() if g in X.columns]
    cen = r["center"]
    g1s_sc, g2m_sc = score(X, e2f, cen), score(X, g2m, cen)
    S_obs = np.column_stack([g1s_sc, g2m_sc])
    rng = np.random.default_rng(np.random.SeedSequence([SEED, ci]))
    pool = L.pool_genes(X, extra_exclude=e2f)
    sets1 = L.pseudo_sets(X, pool, e2f, args.nrand, rng); sets2 = L.pseudo_sets(X, pool, g2m, args.nrand, rng)
    Srand = [np.column_stack([score(X, a, cen), score(X, b, cen)]) for a, b in zip(sets1, sets2)]
    dR = reg[n]
    (u, _), (cd, r2) = parts(r["c"], r["Y"], r["wf"]), parts(r["c"], r["Y"], r["wf"], S_obs)
    assert abs(u[0] - dR) < 1e-6
    Sobs = cd[0] - dR
    rs = np.array(Parallel(n_jobs=args.jobs)(delayed(cond_dR)(r, s) for s in Srand)) - dR
    p = D.cs.stats.perm_p(rs, Sobs, "greater", tie=0.0)
    rows.append(dict(contrast=n, n_genes=len(cg), n_e2f=len(e2f), n_g2m=len(g2m), n_strata=len(r["c"].G), dR=dR, dR_conditioned=cd[0], S=Sobs,
                     rand_S_mean=rs.mean(), rand_S_sd=rs.std(ddof=1), rand_S_q95=np.quantile(rs, 0.95), p=p, n_rand=len(rs),
                     MDE_S=K * rs.std(ddof=1), shrinkage_share=Sobs / abs(dR), MDE_below_half_abs_dR=bool(K * rs.std(ddof=1) < 0.5 * abs(dR))))
    for v, arm in ((r["c"].case, "high"), (r["c"].ref, "low")):
        desc += [dict(contrast=n, quantity="r2_two_scores", arm=arm, gene=g, value=x) for g, x in zip(cg, r2[v])]
    desc += [dict(contrast=n, quantity=q, arm="high-low", gene="all", value=val) for q, val in
             (("dlog_plane_unconditioned", u[1]), ("dlog_plane_conditioned", cd[1]), ("dlog_off_unconditioned", u[2]), ("dlog_off_conditioned", cd[2]),
              ("dlog_plane_minus_off_unconditioned", u[1] - u[2]), ("dlog_plane_minus_off_conditioned", cd[1] - cd[2]),
              ("mean_r2_high", r2[r["c"].case].mean()), ("mean_r2_low", r2[r["c"].ref].mean()))]
    say(n, f"dR {dR:.4f} -> {cd[0]:.4f}, S {Sobs:.4f}, random S {rs.mean():.4f} +- {rs.std(ddof=1):.4f}, p {p:.4f}", f"{time.time() - t0:.0f}s")

res = pd.DataFrame(rows)
o = np.argsort(res.p.to_numpy())
holm = np.minimum(1, np.maximum.accumulate(res.p.to_numpy()[o] * (len(res) - np.arange(len(res)))))
res["p_holm"] = holm[np.argsort(o)]
res["pass"] = (res.p_holm < 0.05) & (res.S > 0)
res.to_csv(OUT / "tests.csv", index=False)
pd.DataFrame(desc).to_csv(OUT / "descriptive.csv", index=False)
npass, med = int(res["pass"].sum()), float(res.shrinkage_share.median())
supported = bool(npass >= 3 and med >= 0.5)
nmde = int(res.MDE_below_half_abs_dR.sum())
verdict = dict(n_pass=npass, median_shrinkage_share=med, H179_1_supported=supported, n_contrasts_MDE_S_below_half_abs_dR=nmde,
               median_MDE_S=float(res.MDE_S.median()), median_half_abs_dR=float((0.5 * res.dR.abs()).median()),
               not_explained_by_linear_cell_cycle_composition=bool((not supported) and nmde >= 3))
pd.Series(verdict).to_csv(OUT / "verdict.csv", header=False)
print(res.round(4).T.to_string()); print(verdict)
say("done A179", f"{time.time() - t0:.0f}s")
