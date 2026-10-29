"""A149 re-scored with the bounded statistic (clocksync `stats.ccorr`, `nan_on_flag`). DESCRIPTIVE / POST HOC (2026-10-29,
roadmap 5.5): registered verdicts are not changed.

Same units (`output/units.pkl`), pooling, scopes, null seeds and pseudo-set seeds as `compute_a149.py --stage score`;
the only change is that a unit whose corrected variance is <= 2 SE is NaN (value and null) and drops out of the pooled
mean. Bootstrap intervals, D149.1 interval and D149.2 heterogeneity are not recomputed (point estimates only).
    python rescore_a149_bounded.py [--jobs J] [--n-perm 10000] [--n-pseudo 200]
Outputs: output/{per_cohort,pooled,tests,pseudo_companion,d149_1}_bounded.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SC = HERE.parents[0] / "scripts"
sys.path.insert(0, str(SC))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import bounded_lib as B  # noqa: E402
import compute_a149 as C  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE / "output"
KINDS = C.KINDS


def nulls(res, n_perm):
    """As `compute_a149.nulls` (same streams), bounded nulls and estimability."""
    out = {k: [] for k in KINDS}
    for r in sorted(res, key=lambda r: r["study"]):
        for ki, k in enumerate(KINDS):
            u = r.get(k)
            if not u:
                continue
            rng = np.random.default_rng([L.SEED, ki, sum(map(ord, r["study"]))])
            nl = B.perm_null(u, n_perm, rng)
            out[k].append((r["study"], r["cancer"], u, nl, B.estimable(u, nl)))
    return out


def old_nulls(res, n_perm):
    return C.nulls(res, n_perm)


def pseudo_one(r, n_pseudo):
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    pb, mal = L.table(r["study"])
    u = L.pair_states(pb, mal, L.NORM)
    G = L.genome_epithelium(r["study"], mal, u["pts"])
    nrm = G[G.compartment == L.NORM]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = nrm[genes].div(nrm.depth, axis=0)
    wg = list(L.q.WF.genes)
    lm = np.log(frac.mean() + 1e-12)[(nrm[genes].mean() >= 1) | frac.columns.isin(wg)]
    rng = np.random.default_rng([L.SEED, 44, sum(map(ord, r["study"]))])  # the A149 stream
    new, old, flag = [], [], []
    for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF):
        p = L.pair_states(G, mal, L.NORM, wf=w, pts=u["pts"])
        if p:
            s = B.stat(p["za"], p["vA"], p["zb"], p["vB"])
            new.append(s["new"]); old.append(s["old"]); flag.append(s["flag"])
        else:
            new.append(np.nan); old.append(np.nan); flag.append(True)
    return r["study"], np.array(new), np.array(old), np.array(flag)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    res = C.load_units()
    nl = nulls(res, args.n_perm)
    per, pooled = [], []
    for k in KINDS:
        for x in nl[k]:
            s = B.stat(x[2]["za"], x[2]["vA"], x[2]["zb"], x[2]["vB"])
            per.append(dict(kind=k, cohort=x[0], cancer=x[1], n=x[2]["n"], old_statistic=s["old"], flagged=s["flag"],
                            bounded_statistic=s["new"], estimable_bounded=x[4], old_gt1=s["old"] > 1,
                            null_mean=float(np.nanmean(x[3])) if np.isfinite(x[3]).any() else np.nan,
                            p_bounded=B.perm_p(x[3], s["new"]) if np.isfinite(s["new"]) else np.nan))
        for name, xs in C.scopes(nl[k]):
            if xs:
                w = [x[2]["n"] for x in xs]
                nul = L.pool(np.array([x[3] for x in xs]), w)
                obs = float(L.pool(np.array([[B.new_stat(x[2])] for x in xs]), w)[0])
                pooled.append(dict(kind=k, hypothesis=C.HYP[k], scope=name, n=sum(w), n_cohorts=len(xs), statistic=obs,
                                   p=B.perm_p(nul, obs), **L.summarise(nul), excess=obs - float(nul.mean()),
                                   cohorts=";".join(x[0] for x in xs)))
    per, P = pd.DataFrame(per), pd.DataFrame(pooled)
    prim = P[(P.scope == "pooled_all") & P.kind.isin(["s1", "s2", "s3"])].copy()
    prim["p_holm"] = C.holm(prim.p.to_numpy())
    prim["verdict_bounded_descriptive"] = np.where(prim.p_holm < 0.05, "pass", "fail")

    # old pooled from the saved registered output, for the side-by-side
    oldP = pd.read_csv(L.OUT / "score" / "pooled.csv")
    oldT = pd.read_csv(L.OUT / "score" / "tests.csv")
    P = P.merge(oldP[["kind", "scope", "n", "n_cohorts", "statistic", "p", "null_mean"]].rename(
        columns={c: "old_" + c for c in ("n", "n_cohorts", "statistic", "p", "null_mean")}), on=["kind", "scope"], how="left")
    prim = prim.merge(oldT[["hypothesis", "n", "n_cohorts", "statistic", "p", "p_holm", "verdict"]].rename(
        columns={c: "old_" + c for c in ("n", "n_cohorts", "statistic", "p", "p_holm", "verdict")}), on="hypothesis", how="left")

    # pseudo-clock companion of H149.1 (S1, estimable cohorts), flagged sets dropped
    s1 = [x for x in nl["s1"] if x[4]]
    ps = Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(next(r for r in res if r["study"] == x[0]), args.n_pseudo) for x in s1)
    d = {s: (a, b, c) for s, a, b, c in ps}
    w = [x[2]["n"] for x in s1]
    pse = L.pool(np.array([d[x[0]][0] for x in s1]), w)
    pse_old = L.pool(np.array([d[x[0]][1] for x in s1]), w)
    flag_share = float(np.mean(np.concatenate([d[x[0]][2] for x in s1])))
    gt1_old = float(np.mean(np.concatenate([d[x[0]][1] for x in s1]) > 1))
    clock = float(prim.loc[prim.kind == "s1", "statistic"].iloc[0])
    ok = np.isfinite(pse)
    comp = dict(n_cohorts=len(s1), n_sets=int(ok.sum()), pseudo_median=float(np.median(pse[ok])), pseudo_q95=float(np.quantile(pse[ok], .95)),
                clock_statistic_bounded=clock, p_pseudo=float((1 + np.sum(pse[ok] >= clock)) / (1 + ok.sum())),
                clock_specific=bool(clock > np.quantile(pse[ok], .95)), share_cohort_set_flagged=flag_share,
                share_cohort_set_old_gt1=gt1_old, pseudo_median_unbounded=float(np.nanmedian(pse_old)),
                pseudo_q95_unbounded=float(np.nanquantile(pse_old, .95)))
    # D149.1 point estimate: S2 retained share (S2 - null2) / (S1' - null1') on the S2 patients, bounded
    sub = {x[0]: x for x in nl["s2"] if x[4]}
    s1sub = [(r["study"], dict(za=r["s3"]["za"], vA=r["s3"]["vA"], zb=r["s3n"]["za"], vB=r["s3n"]["vA"], n=r["s3"]["n"]))
             for r in res if r["study"] in sub and r.get("s3")]
    ww = [u["n"] for _, u in s1sub]
    n1 = L.pool(np.array([B.perm_null(u, args.n_perm, np.random.default_rng([L.SEED, 9, i])) for i, (_, u) in enumerate(s1sub)]), ww)
    o1 = float(L.pool(np.array([[B.new_stat(u)] for _, u in s1sub]), ww)[0])
    p2 = P[(P.kind == "s2") & (P.scope == "pooled_all")].iloc[0]
    d1 = dict(n=sum(ww), n_cohorts=len(s1sub), s1_same_patients=o1, null1=float(np.nanmean(n1)), s2=p2.statistic, null2=p2.null_mean,
              retained=(p2.statistic - p2.null_mean) / (o1 - np.nanmean(n1)))
    per.to_csv(OUT / "per_cohort_bounded.csv", index=False)
    P.to_csv(OUT / "pooled_bounded.csv", index=False)
    prim.to_csv(OUT / "tests_bounded.csv", index=False)
    pd.DataFrame([comp]).to_csv(OUT / "pseudo_companion_bounded.csv", index=False)
    pd.DataFrame([d1]).to_csv(OUT / "d149_1_bounded.csv", index=False)
    summ = per.groupby("kind").agg(units=("cohort", "size"), flagged=("flagged", "sum"), old_gt1=("old_gt1", "sum"),
                                   not_estimable_bounded=("estimable_bounded", lambda s: int((~s).sum())))
    summ.to_csv(OUT / "breakdown_bounded.csv")
    print(summ, flush=True)
    print(prim.round(4).T.to_string(), flush=True)
    print(comp, d1, flush=True)


if __name__ == "__main__":
    main()
