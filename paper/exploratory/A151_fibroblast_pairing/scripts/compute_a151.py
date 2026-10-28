"""A151: tumour vs own adjacent-normal fibroblast pairing, with a pseudo-clock specificity test (endothelial companion).

Registered in paper/docs/registered_experiments.md (A151). Cohorts, tables and statistics of A149 (`a149_lib`).
  power  units and permutation nulls only                                   -> output/power_nulls.csv, output/units.csv
  score  observed pairing, p, pseudo-clock sets, Holm, per cancer           -> output/{tests,pooled,per_cohort,pseudo}.csv

    python compute_a151.py --stage power|score [--jobs J] [--n-perm 10000] [--n-pseudo 200]
"""

from __future__ import annotations

import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER / "exploratory" / "A149_paired_tumour_normal" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
from common import cli  # noqa: E402

OUT = PAPER / "exploratory" / "A151_fibroblast_pairing" / "output"
SEED = 151
PAIRS = {"fibroblast": ("tumour:fibroblast", "adjacent_normal:fibroblast"),
         "endothelial": ("tumour:endothelial", "adjacent_normal:endothelial")}


def genome_pair(study: str, comps: tuple, pts: list) -> pd.DataFrame:
    """Genome-wide pseudobulks of the two compartments for `pts` (labels and cells as `cancer_sc`; no copy-number split
    is needed for stromal classes)."""
    us, _ = L.ST.units(study)
    us = [u for u in L.ST.majority_chemistry(us) if u[1] in set(pts)]
    sums, depth, strat, genes = {}, {}, {}, None
    for sample, patient, kind, load, chem, _ in us:
        X, _g, umi, comp, barcodes, sym = load()
        lab = L.B._labels(comp, kind, sample, barcodes, None)
        genes = pd.Index(sym) if genes is None else genes
        for c in comps:
            i = np.flatnonzero(lab == c)
            if len(i):
                k = (patient, c)
                sums[k] = sums.get(k, 0) + np.asarray(X[i].sum(0)).ravel()
                depth[k] = depth.get(k, 0) + float(umi[i].sum())
                strat[k] = chem
    keys = sorted(sums)
    G = pd.DataFrame(np.vstack([sums[k] for k in keys]), columns=genes)
    G.columns = [L.canonical(str(g)) for g in G.columns]
    G = G.T.groupby(level=0).sum().T
    meta = pd.DataFrame({"donor": [k[0] for k in keys], "compartment": [k[1] for k in keys],
                         "stratum": [strat[k] for k in keys], "depth": [depth[k] for k in keys]})
    return pd.concat([meta, G.reset_index(drop=True)], axis=1)


def _units(study):
    from clocksync.config import limit_threads
    limit_threads(1)
    try:
        pb, _ = L.table(study)
    except Exception as e:  # noqa: BLE001
        return study, {}, repr(e)[:200]
    out = {k: L.pair_states(pb, a, b) for k, (a, b) in PAIRS.items()
           if a in set(pb.compartment) and b in set(pb.compartment)}
    return study, {k: v for k, v in out.items() if v}, ""


def units(jobs):
    res = Parallel(n_jobs=jobs)(delayed(_units)(s) for s in L.cohorts())
    pd.DataFrame([dict(study=s, cancer=L.cancer_of(s), error=e, **{f"n_{k}": (u[k]["n"] if k in u else 0) for k in PAIRS})
                  for s, u, e in res]).to_csv(OUT / "units.csv", index=False)
    return [(s, u) for s, u, e in res if u]


def nulls(res, kind, n_perm):
    out = []
    for s, u in sorted(res):
        if kind in u:
            x = u[kind]
            nl = L.null_of(x, n_perm, np.random.default_rng([SEED, len(kind), sum(map(ord, s))]))
            out.append((s, L.cancer_of(s), x, nl, L.estimable(x, nl)))
    return out


def scopes(lst):
    est = [x for x in lst if x[4]]
    yield "pooled_all", est
    by = {}
    for x in est:
        by.setdefault(x[1], []).append(x)
    for c, xs in sorted(by.items()):
        if sum(x[2]["n"] for x in xs) >= 12:
            yield f"cancer:{c}", xs


def pseudo_one(study, kind, n_pseudo):
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    a, b = PAIRS[kind]
    pb, _ = L.table(study)
    u = L.pair_states(pb, a, b)
    G = genome_pair(study, (a, b), u["pts"])
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    wg = list(L.q.WF.genes)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(wg)]
    rng = np.random.default_rng([SEED, 44, len(kind), sum(map(ord, study))])
    vals = []
    for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF):
        p = L.pair_states(G, a, b, wf=w, pts=u["pts"])
        vals.append(L.obs_of(p) if p else np.nan)
    return study, np.array(vals)


def perm_p(null, obs):
    x = null[np.isfinite(null)]
    return (1 + np.sum(x >= obs)) / (1 + len(x))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    OUT.mkdir(parents=True, exist_ok=True)
    res = units(args.jobs)
    rows, pooled, per = [], [], []
    for kind in PAIRS:
        nl = nulls(res, kind, args.n_perm)
        for name, xs in scopes(nl):
            if not xs:
                continue
            w = [x[2]["n"] for x in xs]
            nul = L.pool(np.array([x[3] for x in xs]), w)
            row = dict(kind=kind, scope=name, n=sum(w), n_cohorts=len(xs), cohorts=";".join(x[0] for x in xs), **L.summarise(nul))
            if args.stage == "score":
                obs = float(L.pool(np.array([[L.obs_of(x[2])] for x in xs]), w)[0])
                row.update(statistic=obs, p=perm_p(nul, obs))
            pooled.append(row)
        if args.stage == "score":
            for x in nl:
                obs = L.obs_of(x[2])
                per.append(dict(kind=kind, cohort=x[0], cancer=x[1], n=x[2]["n"], estimable=x[4], statistic=obs,
                                null_mean=float(np.nanmean(x[3])), p=perm_p(x[3], obs) if np.isfinite(obs) else np.nan))
    P = pd.DataFrame(pooled)
    if args.stage == "power":
        P.to_csv(OUT / "power_nulls.csv", index=False)
        print(P.drop(columns="cohorts").round(4).to_string(), flush=True)
        return
    ref = pd.read_csv(OUT / "power_nulls.csv")
    m = P.merge(ref[["kind", "scope", "null_mean"]], on=["kind", "scope"], suffixes=("", "_power"))
    assert np.allclose(m.null_mean, m.null_mean_power), "null differs from the power step"
    fib = [x for x in nulls(res, "fibroblast", 1) if x[4]]
    ps = dict(Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(x[0], "fibroblast", args.n_pseudo) for x in fib))
    pse = L.pool(np.array([ps[x[0]] for x in fib]), [x[2]["n"] for x in fib])
    clock = float(P[(P.kind == "fibroblast") & (P.scope == "pooled_all")].statistic.iloc[0])
    p2 = (1 + np.sum(pse[np.isfinite(pse)] >= clock)) / (1 + np.isfinite(pse).sum())
    p1 = float(P[(P.kind == "fibroblast") & (P.scope == "pooled_all")].p.iloc[0])
    ps_ = np.array([p1, p2])
    o = np.argsort(ps_)
    holm = np.empty(2)
    holm[o[0]] = min(1, 2 * ps_[o[0]])
    holm[o[1]] = max(holm[o[0]], min(1, ps_[o[1]]))
    T = pd.DataFrame([dict(hypothesis="H151.1", statistic=clock, p=p1, p_holm=holm[0]),
                      dict(hypothesis="H151.2", statistic=clock, pseudo_median=float(np.nanmedian(pse)),
                           pseudo_q95=float(np.nanquantile(pse, .95)), n_sets=int(np.isfinite(pse).sum()), p=p2, p_holm=holm[1])])
    T["verdict"] = np.where(T.p_holm < 0.05, "pass", "fail")
    T.to_csv(OUT / "tests.csv", index=False)
    P.to_csv(OUT / "pooled.csv", index=False)
    pd.DataFrame(per).to_csv(OUT / "per_cohort.csv", index=False)
    pd.DataFrame({x[0]: ps[x[0]] for x in fib}).to_csv(OUT / "pseudo.csv", index=False)
    print(T.round(4).to_string(), flush=True)
    print(P.drop(columns="cohorts").round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
