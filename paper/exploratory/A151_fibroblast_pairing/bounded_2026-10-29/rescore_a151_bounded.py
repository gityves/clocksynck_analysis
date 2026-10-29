"""A151 re-scored with the bounded statistic (clocksync `stats.ccorr`, `nan_on_flag`). DESCRIPTIVE / POST HOC (2026-10-29,
roadmap 5.5): registered verdicts are not changed. Fibroblast pairing (H151.1, H151.2) and the endothelial companion.

Same cohorts, units, null seeds and pseudo-set seeds as `compute_a151.py --stage score`; a unit (or a pseudo set within a
cohort) whose corrected variance is <= 2 SE is NaN and drops out of the pooled mean (weights renormalised). The registered
output is not touched (`F.units` is not called: it rewrites units.csv).
    python rescore_a151_bounded.py [--jobs J] [--n-perm 10000] [--n-pseudo 200]
Outputs: output/{per_cohort,pooled,tests,pseudo_summary}_bounded.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER / "exploratory/A149_paired_tumour_normal/bounded_2026-10-29"))
sys.path.insert(0, str(HERE.parents[0] / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import bounded_lib as B  # noqa: E402
import compute_a151 as F  # noqa: E402
from common import cli  # noqa: E402

L = F.L
OUT = HERE / "output"


def pseudo_one(study, kind, n_pseudo):
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    a, b = F.PAIRS[kind]
    pb, _ = L.table(study)
    u = L.pair_states(pb, a, b)
    G = F.genome_pair(study, (a, b), u["pts"])
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
    rng = np.random.default_rng([F.SEED, 44, len(kind), sum(map(ord, study))])  # the A151 stream
    new, old, flag = [], [], []
    for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF):
        p = L.pair_states(G, a, b, wf=w, pts=u["pts"])
        if p:
            s = B.stat(p["za"], p["vA"], p["zb"], p["vB"])
            new.append(s["new"]); old.append(s["old"]); flag.append(s["flag"])
        else:
            new.append(np.nan); old.append(np.nan); flag.append(True)
    return study, kind, np.array(new), np.array(old), np.array(flag)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    res = [(s, u) for s, u, e in Parallel(n_jobs=args.jobs)(delayed(F._units)(s) for s in L.cohorts()) if u]
    per, pooled, tests, psum = [], [], [], []
    oldP = pd.read_csv(F.OUT / "pooled.csv")
    for kind in F.PAIRS:
        xs_all = []
        for s, u in sorted(res):
            if kind in u:
                x = u[kind]
                nl = B.perm_null(x, args.n_perm, np.random.default_rng([F.SEED, len(kind), sum(map(ord, s))]))
                st = B.stat(x["za"], x["vA"], x["zb"], x["vB"])
                xs_all.append((s, L.cancer_of(s), x, nl, B.estimable(x, nl)))
                per.append(dict(kind=kind, cohort=s, cancer=L.cancer_of(s), n=x["n"], old_statistic=st["old"], old_gt1=st["old"] > 1,
                                flagged=st["flag"], bounded_statistic=st["new"], estimable_bounded=xs_all[-1][4],
                                estimable_registered=bool(L.estimable(x, B.q.perm_null(x['za'], x['vA'], x['zb'], x['vB'], args.n_perm, np.random.default_rng([F.SEED, len(kind), sum(map(ord, s))])))),
                                p_bounded=B.perm_p(nl, st["new"]) if np.isfinite(st["new"]) else np.nan))
        for name, xs in F.scopes(xs_all):
            if not xs:
                continue
            w = [x[2]["n"] for x in xs]
            nul = L.pool(np.array([x[3] for x in xs]), w)
            obs = float(L.pool(np.array([[B.new_stat(x[2])] for x in xs]), w)[0])
            o = oldP[(oldP.kind == kind) & (oldP.scope == name)]
            pooled.append(dict(kind=kind, scope=name, n=sum(w), n_cohorts=len(xs), statistic=obs, p=B.perm_p(nul, obs), **L.summarise(nul),
                               old_n=o.n.iloc[0] if len(o) else np.nan, old_n_cohorts=o.n_cohorts.iloc[0] if len(o) else np.nan,
                               old_statistic=o.statistic.iloc[0] if len(o) else np.nan, old_p=o.p.iloc[0] if len(o) else np.nan,
                               old_null_mean=o.null_mean.iloc[0] if len(o) else np.nan))
    P, PER = pd.DataFrame(pooled), pd.DataFrame(per)
    # pseudo-clock sets for the estimable cohorts of each kind
    jobs = [(s, k) for k in F.PAIRS for s in PER[(PER.kind == k) & PER.estimable_bounded].cohort]
    ps = Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(s, k, args.n_pseudo) for s, k in jobs)
    d = {(s, k): (a, b, c) for s, k, a, b, c in ps}
    for kind in F.PAIRS:
        sel = PER[(PER.kind == kind) & PER.estimable_bounded]
        w = sel.n.to_numpy(float)
        pse = L.pool(np.array([d[(s, kind)][0] for s in sel.cohort]), w)
        pse_old = L.pool(np.array([d[(s, kind)][1] for s in sel.cohort]), w)
        clock = float(P[(P.kind == kind) & (P.scope == "pooled_all")].statistic.iloc[0])
        ok = np.isfinite(pse)
        fl = np.concatenate([d[(s, kind)][2] for s in sel.cohort])
        og = np.concatenate([d[(s, kind)][1] for s in sel.cohort])
        psum.append(dict(kind=kind, n_cohorts=len(sel), n_sets=int(ok.sum()), clock_bounded=clock, pseudo_median=float(np.median(pse[ok])),
                         pseudo_q95=float(np.quantile(pse[ok], .95)), p_pseudo=float((1 + np.sum(pse[ok] >= clock)) / (1 + ok.sum())),
                         share_cohort_set_flagged=float(fl.mean()), share_cohort_set_old_gt1=float(np.nanmean(og > 1)),
                         pseudo_median_unbounded=float(np.nanmedian(pse_old)), pseudo_q95_unbounded=float(np.nanquantile(pse_old, .95))))
    PS = pd.DataFrame(psum)
    # Holm over H151.1 / H151.2 (fibroblast) as registered
    f = P[(P.kind == "fibroblast") & (P.scope == "pooled_all")].iloc[0]
    pf = PS[PS.kind == "fibroblast"].iloc[0]
    ps_ = np.array([f.p, pf.p_pseudo]); o = np.argsort(ps_); h = np.empty(2)
    h[o[0]] = min(1, 2 * ps_[o[0]]); h[o[1]] = max(h[o[0]], min(1, ps_[o[1]]))
    T = pd.DataFrame([dict(hypothesis="H151.1", statistic=f.statistic, p=f.p, p_holm=h[0]),
                      dict(hypothesis="H151.2", statistic=f.statistic, pseudo_median=pf.pseudo_median, pseudo_q95=pf.pseudo_q95,
                           n_sets=pf.n_sets, p=pf.p_pseudo, p_holm=h[1])])
    T["verdict_bounded_descriptive"] = np.where(T.p_holm < 0.05, "pass", "fail")
    old = pd.read_csv(F.OUT / "tests.csv")
    T = T.merge(old[["hypothesis", "statistic", "pseudo_median", "pseudo_q95", "p", "p_holm", "verdict"]].rename(
        columns=lambda c: c if c == "hypothesis" else "old_" + c), on="hypothesis")
    PER.to_csv(OUT / "per_cohort_bounded.csv", index=False)
    P.to_csv(OUT / "pooled_bounded.csv", index=False)
    PS.to_csv(OUT / "pseudo_summary_bounded.csv", index=False)
    T.to_csv(OUT / "tests_bounded.csv", index=False)
    summ = PER.groupby("kind").agg(units=("cohort", "size"), flagged=("flagged", "sum"), old_gt1=("old_gt1", "sum"),
                                   estimable_registered=("estimable_registered", "sum"), estimable_bounded=("estimable_bounded", "sum"))
    summ.to_csv(OUT / "breakdown_bounded.csv")
    print(summ, flush=True); print(T.round(4).T.to_string(), flush=True); print(PS.round(4).T.to_string(), flush=True)
    print(P.drop(columns=[]).round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
