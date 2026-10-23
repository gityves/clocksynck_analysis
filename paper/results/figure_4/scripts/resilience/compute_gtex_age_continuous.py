"""A94b: H94.1 with age as a continuous variable (post hoc amendment of A94; docs/PLAN_internal_desync.md (archive)).

Estimand: per tissue and cohort, the OLS slope of the per-donor squared amplitude r2 = |z|^2 - v on age (per decade)
divided by the tissue's mean r2; median over tissues; primary = mean of the fast-death and ventilated medians.
Null: age permuted among donors within Hardy class x sex (10,000; the same permutation for all tissues of a cohort),
two-sided p. Sensitivities: S1 rank version, S2 Hardy 2 only, S3 per cohort, S4 leave-one-tissue-out; plus the tie-break
spread of the tertile H94.2 and H94.3 (descriptive). Outputs: a94b_continuous.csv, a94b_per_tissue.csv, a94b_null.csv,
a94b_tiebreak_h942_h943.csv (aggregates only). Legacy script: none.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import compute_gtex_age as A  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402


def slopes(c, age: pd.Series, rank: bool = False, keep=None) -> dict:
    """tissue -> relative slope per decade (rank version: slope of the within-tissue rank of r2 on age, per decade)."""
    out = {}
    for t, d in c.t.items():
        if t not in A.BODY:
            continue
        a = age.reindex(d["subj"]).to_numpy(float)
        r = d["r2"]
        m = np.isfinite(a) & (np.ones(len(r), bool) if keep is None else np.isin(d["subj"], keep))
        if m.sum() < A.MIN_N:
            continue
        a, r = a[m] / 10, r[m]
        if rank:
            r = pd.Series(r).rank().to_numpy()
        a = a - a.mean()
        b = float((a * (r - r.mean())).sum() / (a * a).sum())
        if rank:
            out[t] = b
        elif r.mean() > 0:
            out[t] = b / r.mean()
    return out


def permute(age: pd.Series, strata: pd.Series, rng) -> pd.Series:
    out = age.copy()
    for _, idx in age.groupby(strata.reindex(age.index)).groups.items():
        out.loc[idx] = rng.permutation(age.loc[idx].to_numpy())
    return out


def stat(sl: dict, drop=None) -> float:
    v = [x for t, x in sl.items() if t != drop]
    return float(np.median(v)) if v else np.nan


def main(argv=None):
    p = cli.parser(__doc__)
    p.add_argument("--n-perm", type=int, default=None)
    p.add_argument("--draws", type=int, default=30, help="tie-break draws for the tertile H94.2 / H94.3 spread")
    args = p.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(A.FOLDER, run["quick"])
    n_perm = args.n_perm or (200 if run["quick"] else 10000)
    rng = np.random.default_rng(940 + args.seed)
    s = A.load_states(A.BODY + ["Pancreas"] + A._brain_tissues(), args.subsample, args.seed).dropna(subset=["tod"])
    brain = sorted(t for t in s.tissue.unique() if t.startswith("Brain"))
    C = {"fast": A.Cohort(s[s.hardy.isin([1, 2])], True), "vent": A.Cohort(s[s.hardy == 0], False)}
    age, strata = {}, {}
    for k, c in C.items():
        dn = c.donors.dropna(subset=["age"])
        age[k] = dn.age.astype(float)
        strata[k] = dn.hardy.astype(str) + "|" + dn.sex.astype(str)
    hard2 = C["fast"].donors.index[C["fast"].donors.hardy == 2].to_numpy()
    rows, per_t, nulls = [], [], {}

    def run_one(name, fn, cohorts):
        """fn(cohort, age) -> {tissue: slope}; statistic per cohort and combined; permutation null."""
        obs = {k: fn(C[k], age[k]) for k in cohorts}
        for k in cohorts:
            for t, x in obs[k].items():
                per_t.append(dict(analysis=name, cohort=k, tissue=t, slope=x))
        o = {k: stat(obs[k]) for k in cohorts}
        null = {k: np.empty(n_perm) for k in cohorts}
        for i in range(n_perm):
            for k in cohorts:
                null[k][i] = stat(fn(C[k], permute(age[k], strata[k], rng)))
        comb = np.nanmean([o[k] for k in cohorts])
        ncomb = np.nanmean([null[k] for k in cohorts], axis=0)
        for k in cohorts:
            rows.append(dict(analysis=name, scope=k, tissues=len(obs[k]), stat=o[k], p_two=A.p_two(o[k], null[k]),
                             null_sd=float(null[k].std())))
        rows.append(dict(analysis=name, scope="combined", tissues=sum(len(v) for v in obs.values()), stat=comb,
                         p_two=A.p_two(comb, ncomb), null_sd=float(ncomb.std())))
        nulls[name] = ncomb
        return obs

    obs = run_one("H94b.1", lambda c, a: slopes(c, a), ["fast", "vent"])
    run_one("S1_rank", lambda c, a: slopes(c, a, rank=True), ["fast", "vent"])
    # S2: Hardy 2 only (fast deaths); permutation within sex among those donors
    c2, a2 = C["fast"], age["fast"][age["fast"].index.isin(hard2)]
    st2 = strata["fast"].reindex(a2.index).str.split("|").str[1]
    o2 = stat(slopes(c2, a2, keep=hard2))
    n2 = np.array([stat(slopes(c2, permute(a2, st2, rng), keep=hard2)) for _ in range(n_perm)])
    rows.append(dict(analysis="S2_hardy2", scope="fast", tissues=len(slopes(c2, a2, keep=hard2)), stat=o2,
                     p_two=A.p_two(o2, n2), null_sd=float(n2.std())))
    # S4: leave-one-tissue-out range of the primary
    allt = sorted(set(obs["fast"]) | set(obs["vent"]))
    loo = [np.nanmean([stat(obs["fast"], drop=t), stat(obs["vent"], drop=t)]) for t in allt]
    rows.append(dict(analysis="S4_leave_one_tissue_out", scope="combined", tissues=len(allt), stat=float(np.median(loo)),
                     p_two=np.nan, null_sd=np.nan, lo=float(np.min(loo)), hi=float(np.max(loo))))
    pd.DataFrame(rows).to_csv(out / "a94b_continuous.csv", index=False)
    pd.DataFrame(per_t).to_csv(out / "a94b_per_tissue.csv", index=False)
    pd.DataFrame({k: v for k, v in nulls.items()}).to_csv(out / "a94b_null.csv", index=False)

    # tie-break spread of the tertile H94.2 (wall-time alignment) and H94.3 (brain-body coherence), fast deaths
    c, tb = C["fast"], []
    for d in range(args.draws):
        lab = A.tertiles(age["fast"], np.random.default_rng(2000 + d))
        tb.append(dict(draw=d,
                       H94_2=A.contrast(c, lab, 2, 0, A.walltime_R, A.BODY, 94),
                       H94_3=A.contrast(c, lab, 2, 0, lambda cc, dd: A.coherence(cc, dd, brain, A.BODY), None, 94,
                                        per_tissue=False)))
    tb = pd.DataFrame(tb)
    tb.to_csv(out / "a94b_tiebreak_h942_h943.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string())
    print(tb.describe().round(4).to_string())


if __name__ == "__main__":
    main()
