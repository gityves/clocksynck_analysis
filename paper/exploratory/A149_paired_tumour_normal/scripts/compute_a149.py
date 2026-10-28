"""A149: does the malignant clock follow the patient's own normal epithelium in single cells, and through the stroma?

Registered in paper/docs/registered_experiments.md (A149, notes 1-2). Stages:
  tables  build (new cohorts) or load the cohort tables, list the units of every cohort        -> output/units.csv
  power   permutation nulls only (no observed statistic): null mean / SD / q95 / MDE           -> output/power/*.csv
  score   observed S1-S3, p, Holm, pseudo-clock companion, D149.1 retained share, D149.2 per cancer -> output/score/*.csv
The score stage regenerates the nulls from the same seeds and checks them against output/power.

    python compute_a149.py --stage tables|power|score [--jobs J] [--n-perm 10000] [--n-boot 2000] [--n-pseudo 200]
"""

from __future__ import annotations

import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
from common import cli  # noqa: E402

KINDS = ("s1", "s2", "s3", "s3n")
HYP = {"s1": "H149.1", "s2": "H149.2", "s3": "H149.3", "s3n": "H149.3 companion (normal epithelium vs host)"}
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f} s]", *a, flush=True)


def _units(study):
    from clocksync.config import limit_threads
    limit_threads(1)
    try:
        return L.units(study)
    except Exception as e:  # noqa: BLE001  (a cohort that cannot be built is reported, not fatal)
        return dict(study=study, cancer=L.cancer_of(study), error=repr(e)[:200])


def stage_tables(args):
    res = Parallel(n_jobs=args.jobs)(delayed(_units)(s) for s in L.cohorts())
    with open(L.OUT / "units.pkl", "wb") as f:
        pickle.dump(res, f)
    rows = [dict(study=r["study"], cancer=r["cancer"], new=r["study"] in L.NEW, malignant=r.get("malignant", ""),
                 error=r.get("error", ""), **{f"n_{k}": (r[k]["n"] if r.get(k) else 0) for k in KINDS}) for r in res]
    t = pd.DataFrame(rows)
    t.to_csv(L.OUT / "units.csv", index=False)
    log(t.to_string())
    log("patients S1:", int(t.n_s1.sum()), "S2:", int(t.n_s2.sum()), "cohorts S1:", int((t.n_s1 > 0).sum()))


def load_units():
    with open(L.OUT / "units.pkl", "rb") as f:
        return [r for r in pickle.load(f) if "error" not in r]


def nulls(res, n_perm):
    """Per kind: list of (study, cancer, unit, null, estimable); each (cohort, kind) has its own seed."""
    out = {k: [] for k in KINDS}
    for ci, r in enumerate(sorted(res, key=lambda r: r["study"])):
        for ki, k in enumerate(KINDS):
            u = r.get(k)
            if not u:
                continue
            rng = np.random.default_rng([L.SEED, ki, sum(map(ord, r["study"]))])
            nl = L.null_of(u, n_perm, rng)
            out[k].append((r["study"], r["cancer"], u, nl, L.estimable(u, nl)))
    return out


def scopes(lst):
    """Pooled scopes: all, each cancer type with >= 12 patients, and the two classes of note 1."""
    est = [x for x in lst if x[4]]
    yield "pooled_all", est
    by = {}
    for x in est:
        by.setdefault(x[1], []).append(x)
    for c, xs in sorted(by.items()):
        if sum(x[2]["n"] for x in xs) >= 12:
            yield f"cancer:{c}", xs
    yield "class:keeping", [x for x in est if x[1] in L.KEEPING]
    yield "class:at_null", [x for x in est if x[1] in L.AT_NULL]


def pooled_null(xs):
    return L.pool(np.array([x[3] for x in xs]), [x[2]["n"] for x in xs])


def stage_power(args):
    res = load_units()
    nl = nulls(res, args.n_perm)
    rows = []
    for k in KINDS:
        for x in nl[k]:
            rows.append(dict(kind=k, scope=x[0], cancer=x[1], n=x[2]["n"], estimable=x[4], **L.summarise(x[3])))
        for name, xs in scopes(nl[k]):
            if xs:
                rows.append(dict(kind=k, scope=name, cancer="", n=sum(x[2]["n"] for x in xs), estimable=True,
                                 n_cohorts=len(xs), **L.summarise(pooled_null(xs))))
    t = pd.DataFrame(rows)
    (L.OUT / "power").mkdir(parents=True, exist_ok=True)
    t.to_csv(L.OUT / "power" / "nulls.csv", index=False)
    log(t[t.scope.str.contains(":|pooled")].round(4).to_string())


def holm(ps):
    order, adj, run = np.argsort(ps), np.empty(len(ps)), 0.0
    for i, j in enumerate(order):
        run = max(run, min(1.0, (len(ps) - i) * ps[j]))
        adj[j] = run
    return adj


def perm_p(null, obs):
    x = null[np.isfinite(null)]
    return (1 + np.sum(x >= obs)) / (1 + len(x))


def pseudo_one(r, n_pseudo):
    """Pooled-ready pseudo-clock statistics of one cohort's S1 unit (note 2 iv)."""
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
    lm = np.log(frac.mean() + 1e-12)[(nrm[genes].mean() >= 1) | frac.columns.isin(wg)]  # pool filter; clock genes kept
    rng = np.random.default_rng([L.SEED, 44, sum(map(ord, r["study"]))])
    out = []
    for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF):
        p = L.pair_states(G, mal, L.NORM, wf=w, pts=u["pts"])
        out.append(L.obs_of(p) if p else np.nan)
    return r["study"], np.array(out)


def stage_score(args):
    res = load_units()
    nl = nulls(res, args.n_perm)
    ref = pd.read_csv(L.OUT / "power" / "nulls.csv")
    out = L.OUT / "score"
    out.mkdir(parents=True, exist_ok=True)
    per, pooled = [], []
    boots = {}
    for k in KINDS:
        for ci, x in enumerate(nl[k]):
            rr = ref[(ref.kind == k) & (ref.scope == x[0])].iloc[0]
            s = L.summarise(x[3])
            assert np.isclose(s["null_mean"], rr.null_mean, equal_nan=True), f"null drift {k} {x[0]}"
            obs = L.obs_of(x[2])
            b = L.boot_of(x[2], args.n_boot, np.random.default_rng([L.SEED, 7, ci, len(k)]))
            boots[(k, x[0])] = b
            per.append(dict(kind=k, cohort=x[0], cancer=x[1], n=x[2]["n"], estimable=x[4], statistic=obs,
                            null_mean=s["null_mean"], null_sd=s["null_sd"], p=perm_p(x[3], obs) if np.isfinite(obs) else np.nan))
        for name, xs in scopes(nl[k]):
            if not xs:
                continue
            w = [x[2]["n"] for x in xs]
            obs = float(L.pool(np.array([[L.obs_of(x[2])] for x in xs]), w)[0])
            nul = pooled_null(xs)
            bt = L.pool(np.array([boots[(k, x[0])] for x in xs]), w)
            ex = bt - nul.mean()
            pooled.append(dict(kind=k, hypothesis=HYP[k], scope=name, n=sum(w), n_cohorts=len(xs), statistic=obs,
                               p=perm_p(nul, obs), **L.summarise(nul), excess=obs - nul.mean(),
                               excess_lo=np.nanquantile(ex, .025), excess_hi=np.nanquantile(ex, .975),
                               excess_se=np.nanstd(ex, ddof=1), cohorts=";".join(x[0] for x in xs)))
    P = pd.DataFrame(pooled)
    prim = P[(P.scope == "pooled_all") & P.kind.isin(["s1", "s2", "s3"])].copy()
    prim["p_holm"] = holm(prim.p.to_numpy())
    prim["verdict"] = np.where(prim.p_holm < 0.05, "pass", "fail")

    # pseudo-clock companion of H149.1
    s1 = [x for x in nl["s1"] if x[4]]
    ps = dict(Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(next(r for r in res if r["study"] == x[0]), args.n_pseudo)
                                          for x in s1))
    pse = L.pool(np.array([ps[x[0]] for x in s1]), [x[2]["n"] for x in s1])
    comp = dict(n_sets=int(np.isfinite(pse).sum()), median=float(np.nanmedian(pse)), q95=float(np.nanquantile(pse, .95)),
                clock_statistic=float(prim.loc[prim.kind == "s1", "statistic"].iloc[0]))
    comp["clock_specific"] = comp["clock_statistic"] > comp["q95"]

    # D149.1 retained share on the S2 patients: (S2 - null2) / (S1' - null1') with S1' on the same patients
    sub = {x[0]: x for x in nl["s2"] if x[4]}
    s1sub = []
    for r in res:
        if r["study"] in sub and r.get("s3"):
            u = dict(za=r["s3"]["za"], vA=r["s3"]["vA"], zb=r["s3n"]["za"], vB=r["s3n"]["vA"], n=r["s3"]["n"])
            s1sub.append((r["study"], u))
    w = [u["n"] for _, u in s1sub]
    n1 = L.pool(np.array([L.null_of(u, args.n_perm, np.random.default_rng([L.SEED, 9, i])) for i, (_, u) in enumerate(s1sub)]), w)
    o1 = float(L.pool(np.array([[L.obs_of(u)] for _, u in s1sub]), w)[0])
    b1 = L.pool(np.array([L.boot_of(u, args.n_boot, np.random.default_rng([L.SEED, 7, i, 99])) for i, (_, u) in enumerate(s1sub)]), w)
    p2 = P[(P.kind == "s2") & (P.scope == "pooled_all")].iloc[0]
    b2 = L.pool(np.array([boots[("s2", s)] for s, _ in s1sub]), w)
    ratio_b = (b2 - p2.null_mean) / (b1 - n1.mean())
    d1 = dict(n=sum(w), s1_same_patients=o1, null1=n1.mean(), s2=p2.statistic, null2=p2.null_mean,
              retained=(p2.statistic - p2.null_mean) / (o1 - n1.mean()),
              lo=np.nanquantile(ratio_b, .025), hi=np.nanquantile(ratio_b, .975))

    # D149.2 heterogeneity across cancers: I^2 of the excess
    nmap = {k: {x[0]: x[2]["n"] for x in nl[k]} for k in KINDS}

    def boot_excess(k, row):
        cs_ = row.cohorts.split(";")
        return L.pool(np.array([boots[(k, s)] for s in cs_]), [nmap[k][s] for s in cs_]) - row.null_mean

    het = []
    for k in ("s1", "s2", "s3"):
        h = dict(kind=k)
        c = P[(P.kind == k) & P.scope.str.startswith("cancer:")]
        if len(c) >= 2:
            wv = 1 / c.excess_se ** 2
            m = np.sum(wv * c.excess) / wv.sum()
            Q = float(np.sum(wv * (c.excess - m) ** 2))
            h.update(n_cancers=len(c), Q=Q, df=len(c) - 1, I2=max(0.0, (Q - (len(c) - 1)) / Q) if Q > 0 else 0.0)
        kp, an = (P[(P.kind == k) & (P.scope == s)] for s in ("class:keeping", "class:at_null"))
        if len(kp) and len(an):
            d = boot_excess(k, kp.iloc[0]) - boot_excess(k, an.iloc[0])
            h.update(keeping_minus_at_null=float(kp.excess.iloc[0] - an.excess.iloc[0]),
                     contrast_lo=float(np.nanquantile(d, .025)), contrast_hi=float(np.nanquantile(d, .975)))
        het.append(h)

    pd.DataFrame(per).to_csv(out / "per_cohort.csv", index=False)
    P.to_csv(out / "pooled.csv", index=False)
    prim.to_csv(out / "tests.csv", index=False)
    pd.DataFrame([comp]).to_csv(out / "pseudo_companion.csv", index=False)
    pd.DataFrame([d1]).to_csv(out / "d149_1.csv", index=False)
    pd.DataFrame(het).to_csv(out / "d149_2_heterogeneity.csv", index=False)
    log(prim[["hypothesis", "n", "n_cohorts", "statistic", "null_mean", "null_q95", "mde", "p", "p_holm", "verdict"]].round(4).to_string())
    log(comp, d1)
    log(P[P.scope != "pooled_all"][["kind", "scope", "n", "statistic", "null_mean", "excess", "excess_lo", "excess_hi", "p"]].round(3).to_string())
    log(pd.DataFrame(het).round(3).to_string())


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["tables", "power", "score"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-boot", type=int, default=2_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    L.OUT.mkdir(parents=True, exist_ok=True)
    {"tables": stage_tables, "power": stage_power, "score": stage_score}[args.stage](args)


if __name__ == "__main__":
    main()
