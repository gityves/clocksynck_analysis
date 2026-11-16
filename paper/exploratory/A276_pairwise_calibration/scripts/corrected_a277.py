"""A277: the pairwise test with each gene set read against its own no-pairing null (after A276).

Registered in paper/docs/registered_experiments.md (A277). For each gene set j (clock, 200 matched, 200 responsive per
atlas, exactly as A272 / A274) and statistic (T_SV, D of A272; D, T_SV, T_CX of A274): excess_j = T_j - mu_j, mu_j the
mean of T_j over placebo draws (donors permuted within compartment x stratum).

  calibration  mu from N1 (500 draws, seeds [277, 1, d]); in each of N2 (500 draws, seeds [277, 2, d]) the clock's excess
               against the q95 of the companions' excess. V277.1: rate <= 0.075 per combination; V277.2: mean <= 0.065.
  restated     only if both pass: observed excess with mu from N1 and N2 pooled; H277.1-5 = H272.1, H272.2, H274.1-3.

    python corrected_a277.py [--draws 500] [--jobs 6]
Outputs (../output/a277/): calibration.csv, null_draws.npz (statistics per draw and set), tests.csv, per_atlas.csv, mu.csv.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import calibrate_a276 as C  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

SEED = 277
TESTS = {"H277.1": ("H272.1 T_SV", "both"), "H277.2": ("H272.2 D", "both"), "H277.3": ("H274.1 D", "both"),
         "H277.4": ("H274.2 T_SV", "matched"), "H277.5": ("H274.3 T_CX", "both")}


def load_all(n_pseudo):
    data, groups = {}, {}
    for study, names in (("A272", C.A272_ATLASES), ("A274", C.A274_ATLASES)):
        data[study] = {}
        for n in names:
            key, S = C.load(study, n, n_pseudo)
            pairs = C.A272.pair_index(key)
            meta = pd.DataFrame([dict(atlas=n, a=p["a"], b=p["b"], family_a=C.A150.FAM.family(p["a"]),
                                      family_b=C.A150.FAM.family(p["b"])) for p in pairs])
            meta["klass"] = [C.klass(study, n, fa, fb) for fa, fb in zip(meta.family_a, meta.family_b)]
            data[study][n] = (key, S, meta)
            C.log(study, n, "rows", len(key), "pairs", len(pairs))
        cnt = pd.concat([m for _, _, m in data[study].values()]).groupby(["atlas", "klass"]).size().unstack(fill_value=0)
        col = lambda c: cnt[c] if c in cnt else pd.Series(0, index=cnt.index)  # noqa: E731
        if study == "A272":
            groups[study] = dict(sv=names, ly=[n for n in names if col("LY_X")[n] > 0])
        else:
            body = [n for n in names if n in C.A274.BODY]
            sv = [n for n in body if col("SV_X")[n] >= C.A274.MIN_PAIRS]
            groups[study] = dict(sv=sv, both=[n for n in sv if col("LY_X")[n] >= C.A274.MIN_PAIRS],
                                 cx=[n for n in names if n in C.A274.CORTEX and col("CX_diff")[n] >= C.A274.MIN_PAIRS])
    return data, groups


def draw(seed, data, groups):
    """(statistics {test: vector over sets}, class medians {(study, atlas, class): vector}); seed None = real pairing."""
    warnings.filterwarnings("ignore")
    rng = None if seed is None else np.random.default_rng(seed)
    cm = {}
    for study, atl in data.items():
        rows, metas = [], []
        for n, (key, S, meta) in atl.items():
            pairs = C.A272.pair_index(key if rng is None else C.A272.placebo(key, rng))
            rows.append(C.pair_abs(S, pairs))
            metas.append(meta)
        cm[study] = C.class_medians(np.concatenate(rows, 1), pd.concat(metas, ignore_index=True))
    flat = {(s, n, k): v for s in cm for (n, k), v in cm[s].items()}
    return C.statistics(cm, groups), flat


def run(seeds, data, groups, jobs):
    res = Parallel(n_jobs=jobs, verbose=0)(delayed(draw)(s, data, groups) for s in seeds)
    st = {t: np.array([r[0][t] for r in res]) for t in res[0][0]}
    cm = {k: np.array([r[1][k] for r in res]) for k in res[0][1]}
    return st, cm


def compare(ex, n_pseudo):
    """Clock excess against each family's excess: (q95, p, above) per family."""
    out = {}
    for kind, sl in (("matched", slice(1, 1 + n_pseudo)), ("responsive", slice(1 + n_pseudo, 1 + 2 * n_pseudo))):
        comp = ex[sl][np.isfinite(ex[sl])]
        q = float(np.quantile(comp, 0.95))
        out[kind] = dict(q95=q, median=float(np.median(comp)), p=float((1 + np.sum(comp >= ex[0])) / (1 + len(comp))),
                         above=bool(ex[0] > q))
    return out


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=500)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(C.AMENDMENT) / "a277"
    out.mkdir(parents=True, exist_ok=True)
    C.A150.ATLASES.update(C.A274.ALL)
    data, groups = load_all(args.n_pseudo)
    n1, cm1 = run([[SEED, 1, d] for d in range(args.draws)], data, groups, args.jobs)
    C.log("N1 done")
    n2, cm2 = run([[SEED, 2, d] for d in range(args.draws)], data, groups, args.jobs)
    C.log("N2 done")
    np.savez_compressed(out / "null_draws.npz", **{f"N1 {t}": v for t, v in n1.items()}, **{f"N2 {t}": v for t, v in n2.items()})

    cal = []
    for t in n1:
        mu1 = np.nanmean(n1[t], 0)
        rej = {"matched": [], "responsive": []}
        for v in n2[t]:
            for kind, r in compare(v - mu1, args.n_pseudo).items():
                rej[kind].append(r["above"])
        for kind in rej:
            cal.append(dict(test=t, family=kind, draws=len(rej[kind]), rejection_rate=float(np.mean(rej[kind]))))
    cal = pd.DataFrame(cal)
    cal["V277.1_pass"] = cal.rejection_rate <= 0.075
    cal.to_csv(out / "calibration.csv", index=False)
    v1, v2 = bool(cal["V277.1_pass"].all()), bool(cal.rejection_rate.mean() <= 0.065)
    print(cal.round(3).to_string())
    print("V277.1", "PASS" if v1 else "FAIL", "| V277.2 mean", round(cal.rejection_rate.mean(), 4), "PASS" if v2 else "FAIL")
    if not (v1 and v2):
        C.log("calibration failed: restated tests not read")
        return

    obs, cmo = draw(None, data, groups)
    mu = {t: np.nanmean(np.concatenate([n1[t], n2[t]]), 0) for t in n1}
    pd.DataFrame({t: m for t, m in mu.items()}).to_csv(out / "mu.csv", index_label="set")
    tests = []
    for h, (t, rule) in TESTS.items():
        ex = obs[t] - mu[t]
        r = compare(ex, args.n_pseudo)
        row = dict(test=h, restates=t, clock=float(obs[t][0]), clock_null_mu=float(mu[t][0]), clock_excess=float(ex[0]))
        for kind in r:
            row.update({f"{kind}_{k}": v for k, v in r[kind].items()})
        row["passes"] = r["matched"]["above"] and (rule == "matched" or r["responsive"]["above"])
        tests.append(row)
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "tests.csv", index=False)
    print(tests.round(3).T.to_string())

    pa = []
    for k, v in cmo.items():
        m = np.nanmean(np.concatenate([cm1[k], cm2[k]]), 0)
        r = compare(v - m, args.n_pseudo)
        row = dict(study=k[0], atlas=k[1], klass=k[2], clock=float(v[0]), clock_excess=float(v[0] - m[0]))
        for kind in r:
            row.update({f"{kind}_{x}": y for x, y in r[kind].items()})
        pa.append(row)
    pd.DataFrame(pa).to_csv(out / "per_atlas.csv", index=False)
    C.log("done")


if __name__ == "__main__":
    main()
