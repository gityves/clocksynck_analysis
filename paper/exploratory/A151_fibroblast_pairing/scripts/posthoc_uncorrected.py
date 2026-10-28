"""A149 note 5 / A151 note 2 (post hoc): clock vs pseudo-clock pairing with the uncorrected |complex correlation|.

Same cohorts, patients, genome tables and pseudo-set seeds as the registered companions (A149 `pseudo_one`, A151
`pseudo_one`); per-cohort values saved. Outputs: paper/exploratory/A151_fibroblast_pairing/output/posthoc_*.csv
    python posthoc_uncorrected.py [--jobs J] [--n-pseudo 200]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import compute_a151 as F  # noqa: E402
from common import cli  # noqa: E402

L = F.L
EPI = ("epithelium", None)  # A149 S1: malignant compartment of the cohort vs adjacent-normal epithelium


def raw(u) -> float:
    za, zb = u["za"] - u["za"].mean(), u["zb"] - u["zb"].mean()
    return float(abs(np.sum(za * np.conj(zb))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2)))


def one(study, kind, n_pseudo):
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    pb, mal = L.table(study)
    if kind == "epithelium":
        a, b = mal, L.NORM
        seed = [L.SEED, 44, sum(map(ord, study))]  # A149 pseudo_one stream
        u = L.pair_states(pb, a, b)
        G = L.genome_epithelium(study, mal, u["pts"])
    else:
        a, b = F.PAIRS[kind]
        seed = [F.SEED, 44, len(kind), sum(map(ord, study))]  # A151 pseudo_one stream
        u = L.pair_states(pb, a, b)
        G = F.genome_pair(study, (a, b), u["pts"])
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
    rng = np.random.default_rng(seed)
    rv, cv = [], []
    for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF):
        p = L.pair_states(G, a, b, wf=w, pts=u["pts"])
        rv.append(raw(p) if p else np.nan)
        cv.append(L.obs_of(p) if p else np.nan)
    return dict(kind=kind, cohort=study, n=u["n"], clock_raw=raw(u), clock_corrected=L.obs_of(u),
                pseudo_raw=np.array(rv), pseudo_corrected=np.array(cv))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    a149 = pd.read_csv(L.OUT / "score" / "per_cohort.csv")
    a151 = pd.read_csv(F.OUT / "per_cohort.csv")
    jobs = [(r.cohort, "epithelium") for r in a149[(a149.kind == "s1") & a149.estimable].itertuples()]
    jobs += [(r.cohort, r.kind) for r in a151[a151.estimable].itertuples()]
    res = Parallel(n_jobs=args.jobs)(delayed(one)(s, k, args.n_pseudo) for s, k in jobs)
    per, pooled = [], []
    for kind in ("epithelium", "fibroblast", "endothelial"):
        rs = [r for r in res if r["kind"] == kind]
        w = np.array([r["n"] for r in rs], float)
        clock = float(np.sum(w * [r["clock_raw"] for r in rs]) / w.sum())
        P = np.array([r["pseudo_raw"] for r in rs])
        Pw = np.where(np.isfinite(P), w[:, None], 0)
        pse = np.nansum(P * Pw, 0) / Pw.sum(0)
        ok = np.isfinite(pse)
        pooled.append(dict(kind=kind, n=int(w.sum()), n_cohorts=len(rs), clock_raw=clock, pseudo_median=float(np.median(pse[ok])),
                           pseudo_q95=float(np.quantile(pse[ok], .95)), n_sets=int(ok.sum()),
                           p=(1 + np.sum(pse[ok] >= clock)) / (1 + ok.sum()),
                           share_corrected_pseudo_gt1=float(np.nanmean(np.concatenate([r["pseudo_corrected"] for r in rs]) > 1))))
        for r in rs:
            per.append(dict(kind=kind, cohort=r["cohort"], n=r["n"], clock_raw=r["clock_raw"], clock_corrected=r["clock_corrected"],
                            pseudo_raw_median=float(np.nanmedian(r["pseudo_raw"])), pseudo_raw_q95=float(np.nanquantile(r["pseudo_raw"], .95)),
                            pseudo_corrected_median=float(np.nanmedian(r["pseudo_corrected"])),
                            share_corrected_pseudo_gt1=float(np.nanmean(r["pseudo_corrected"] > 1))))
    out = F.OUT
    pd.DataFrame(pooled).to_csv(out / "posthoc_pooled.csv", index=False)
    pd.DataFrame(per).to_csv(out / "posthoc_per_cohort.csv", index=False)
    print(pd.DataFrame(pooled).round(4).to_string(), flush=True)
    print(pd.DataFrame(per).round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
