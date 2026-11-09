"""A235: are metabolic organs phase-shifted relative to the body clock in T2D and obesity? GTEx, A164 matching, A233 estimation.

Registered in paper/docs/registered_experiments.md (A235). Mouse prediction (high-fat diet): metabolic organs advanced relative
to the host in cases (offset_diff_h > 0, `Dataset.contrast`: case minus control, > 0 advanced). Labels and matching as A164:
T2D vs non-diabetic (up to 3 controls), obesity BMI >= 30 vs 18.5 <= BMI < 25 without T2D (1:1); exact sex and Hardy class,
age caliper 5 y, nearest age then ischaemic time. Machinery: ../../A231_age_coupling/scripts/matched_label.py. Two stages:
  power  permutation null only (no observed statistic): MDE_80 of the module offset and d* per cohort and pooled, module-minus-
         companion offset (fast), for the A164 host and the sensitivity host; matching and balance ->
         output/{power,power_mc,matching,balance,balance_summary,eligible_summary}.csv.
  score  (run only after the power note is dated) H235.1 (fast module offset, one-sided p, Holm over conditions, pseudo-clock
         and companion), H235.2 (vent module offset, descriptive), descriptives and the sensitivity host ->
         output/score_*.csv. `--placebo`: random within-set case labels (code test, output/placebo/).

    python compute_a235.py --stage power [--subsample N] [--conditions t2d obese] [--jobs 8] [--n-perm 1000]
Controlled access: conditions, BMI, exact age and ischaemic time stay in memory; only aggregates are written.

Choices not fixed by the registration
- Module pooling: inverse-null-variance weights over the module tissues of a cohort (A233); the pooled (vent + fast) value
  weights cohorts by inverse null variance (descriptive only).
- One-sided p (H235.1 and module minus companion): add-one, upper tail of the fast-cohort permutation null of the same pooled
  statistic. The two-sided p is A233's (`pooled_tests`).
- Pseudo-clock companion: 50 sets on the fast cohort only (H235.1 is a fast test); q95 of |pooled module offset|, pooled with
  the clock null weights (A233 `pseudo_summary`).
- Module log z2 ratio: unweighted mean over module tissues of log(z2_case / z2_control), per cohort (descriptive, no null).
- Sensitivity host: A164 anchors minus "Muscle - Skeletal" and the two adiposes (the module tissues among the anchors); same
  matched sets, its own permutations and bootstraps (suffix `_noMET`).
"""

from __future__ import annotations

import sys
import time
import warnings
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A231_age_coupling/scripts"))
import matched_label as ML  # noqa: E402
from matched_label import A, B, S, say  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
MET, COMP = list(A.MET), list(A.COMP)
NOMET = tuple(t for t in S.TISSUES if t not in MET)  # sensitivity host
HOSTS = {"": None, "_noMET": NOMET}
CONDITIONS = ("t2d", "obese")


def design_of(name: str) -> ML.Design:
    if name == "t2d":
        pools, match = (lambda q: (q.t2d == 1, q.t2d == 0)), ML.match_ratio
    else:
        pools = lambda q: ((q.bmi >= 30) & (q.t2d == 0), (q.bmi >= 18.5) & (q.bmi < 25) & (q.t2d == 0))
        match = partial(ML.match_ratio, ratio=1)
    return ML.Design(name, pools, match, min_arm=0, min_cases=10, balance=("age", "ischaemic_h", "bmi", "t2d", "male"),
                     modules={"module": MET, "companion": COMP})


def cohort_pool(arr, tw):
    """Pooled module value of one cohort (rows: replicates or a single observation)."""
    return B.wmean(np.atleast_2d(arr), tw)


def module_minus_companion(obs, nulls, msel, csel, c="fast"):
    """(observed, null draws) of module minus companion offset in cohort c, offset weights from the clock null."""
    if c not in msel or c not in csel:
        return np.nan, np.array([np.nan])
    twm = B.tissue_weights(nulls[c][:, msel[c], 1])
    twc = B.tissue_weights(nulls[c][:, csel[c], 1])
    null = cohort_pool(nulls[c][:, msel[c], 1], twm) - cohort_pool(nulls[c][:, csel[c], 1], twc)
    o = float(cohort_pool(obs[c][msel[c], 1], twm)[0] - cohort_pool(obs[c][csel[c], 1], twc)[0])
    return o, null


def p_upper(est, null):
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    return float((1 + np.sum(null >= est)) / (1 + len(null))) if np.isfinite(est) and len(null) else np.nan


def power(args, out, quick, cohorts):
    prows, mrows, brows, crows = [], [], [], []
    for name in args.conditions:
        d = design_of(name)
        for suf, anc in HOSTS.items():
            plans, dss, anchors, pbs, M, Bt = ML.build(d, cohorts, quick=quick, tag="gtex_a235", anchor_set=anc)
            if not suf:
                mrows.append(M)
                brows.append(Bt)
                say(f"{name}: " + "; ".join(f"{c} {len(p)} eligible tissues, {int(M[(M.cohort == c) & M.eligible].n_sets.sum())}"
                                            f" sets" for c, p in plans.items()))
            if not plans:
                continue
            nulls = ML.nulls_of(name, plans, dss, out, quick, args.n_perm, args.jobs, args.seed, suffix=suf)
            msel, csel = ML.selection(plans, MET), ML.selection(plans, COMP)
            if msel:
                prows += [dict(r, host=suf or "A164") for r in ML.power_rows(nulls, plans, M, msel, "module", name, args.seed)]
            _, null = module_minus_companion({c: np.zeros((len(plans[c]), 2)) for c in plans}, nulls, msel, csel)
            sd = float(np.nanstd(null, ddof=1))
            crows.append(dict(design=name, host=suf or "A164", readout="module_minus_companion_offset_h", set="fast",
                              null_sd=sd, mde_80=ML.MDE_K * sd))
            pd.DataFrame(prows).to_csv(out / "power.csv", index=False)  # partial results survive a failure
    M, Bt = pd.concat(mrows, ignore_index=True), pd.concat(brows, ignore_index=True)
    M.to_csv(out / "matching.csv", index=False)
    Bt.to_csv(out / "balance.csv", index=False)
    ML.balance_summary(Bt, M).to_csv(out / "balance_summary.csv", index=False)
    ML.eligibility_summary(M).to_csv(out / "eligible_summary.csv", index=False)
    pd.DataFrame(crows).to_csv(out / "power_mc.csv", index=False)
    P = pd.DataFrame(prows)
    P.to_csv(out / "power.csv", index=False)
    pd.set_option("display.width", 250)
    say("power\n" + P[P.readout != "Q"].drop(columns=["bound"], errors="ignore").round(4).to_string(index=False))
    say("module minus companion\n" + pd.DataFrame(crows).round(4).to_string(index=False))


def score(args, out, quick, cohorts):
    n_pseudo = 5 if quick else 50
    pooled, tissues, tests, z2rows = [], [], [], []
    for name in args.conditions:
        d = design_of(name)
        for suf, anc in HOSTS.items():
            host = suf or "A164"
            plans, dss, anchors, pbs, M, _ = ML.build(d, cohorts, quick=quick, tag="gtex_a235", anchor_set=anc)
            if not plans:
                continue
            if args.placebo:
                plans = ML.placebo_plans(plans, args.seed + 7)
            nulls = ML.nulls_of(name, plans, dss, out, quick, args.n_perm, args.jobs, args.seed, suffix=suf)
            boots = ML.boots_of(name, pbs, anchors, plans, out, quick, args.n_boot, args.jobs, args.seed, suffix=suf)
            obs, T = ML.observed(dss, plans)
            T.insert(0, "host", host)
            T.insert(0, "design", name)
            tissues.append(T)
            msel, csel = ML.selection(plans, MET), ML.selection(plans, COMP)
            if not msel:
                continue
            P = ML.pooled_tests(obs, nulls, boots, msel, "module", name)
            P.insert(0, "host", host)
            pooled.append(P)
            r = dict(condition=name, host=host)
            if "fast" in msel:
                tw = B.tissue_weights(nulls["fast"][:, msel["fast"], 1])
                fo = P[(P.set == "fast") & (P.readout == "offset_h")].iloc[0]
                r.update(fast_offset_h=fo.estimate, fast_lo90=fo.lo90, fast_hi90=fo.hi90, fast_mde_80=fo.mde_80,
                         fast_p_one_sided=p_upper(fo.estimate, cohort_pool(nulls["fast"][:, msel["fast"], 1], tw)),
                         fast_p_two_sided=fo.p_perm, fast_equivalent_1h=bool(fo.lo90 > -1 and fo.hi90 < 1))
                o, null = module_minus_companion(obs, nulls, msel, csel)
                r.update(fast_module_minus_companion_h=o, fast_module_minus_companion_p=p_upper(o, null))
                if not suf:
                    fplans = {"fast": plans["fast"]}
                    ps = ML.pseudo_clock(fplans, {"fast": anchors["fast"]}, n_pseudo, args.jobs, args.seed)
                    spec, off_v, _ = ML.pseudo_summary(ps, {"fast": obs["fast"]}, {"fast": nulls["fast"]},
                                                       {"fast": msel["fast"]}, {"fast": msel["fast"]})
                    pd.DataFrame(dict(offset_h=off_v)).to_csv(out / f"score_pseudo_{name}.csv", index=False)
                    r.update(offset_pseudo_abs_q95=spec["offset_pseudo_abs_q95"],
                             offset_clock_specific=spec["offset_clock_specific"])
            if "vent" in msel:
                vo = P[(P.set == "vent") & (P.readout == "offset_h")].iloc[0]
                r.update(vent_offset_h=vo.estimate, vent_lo90=vo.lo90, vent_hi90=vo.hi90, vent_mde_80=vo.mde_80,
                         vent_p_two_sided=vo.p_perm)
            po = P[(P.set == "pooled") & (P.readout == "offset_h")]
            if len(po):
                r.update(pooled_offset_h=po.iloc[0].estimate, pooled_lo90=po.iloc[0].lo90, pooled_hi90=po.iloc[0].hi90)
            tests.append(r)
            if not suf:
                mt = T[T.tissue.isin(MET)]
                for c, g in mt.groupby("cohort"):
                    lz = np.log(g.z2_case.astype(float) / g.z2_control.astype(float))
                    z2rows.append(dict(condition=name, cohort=c, n_tissues=len(g), module_log_z2_ratio=float(np.nanmean(lz))))
            say(f"{name} host {host} scored")
    pd.concat(pooled, ignore_index=True).to_csv(out / "score_pooled.csv", index=False)
    pd.concat(tissues, ignore_index=True).to_csv(out / "score_tissues.csv", index=False)
    pd.DataFrame(z2rows).to_csv(out / "score_z2.csv", index=False)
    Tt = pd.DataFrame(tests)
    main = Tt.host == "A164"
    if "fast_p_one_sided" in Tt:
        Tt.loc[main, "fast_p_holm"] = ML.holm(Tt.loc[main, "fast_p_one_sided"].to_numpy())
    Tt.to_csv(out / "score_tests.csv", index=False)
    pd.set_option("display.width", 250)
    say("tests\n" + Tt.round(4).T.to_string())


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--conditions", nargs="*", default=list(CONDITIONS), choices=list(CONDITIONS))
    ap.add_argument("--placebo", action="store_true", help="score stage code test on random within-set labels")
    args = ap.parse_args(argv)
    quick = cli.setup(args)["quick"]
    out = OUT / ("subsample" if quick else ".")
    if args.placebo:
        out = OUT / "placebo"
    (out / ML.DRAW_SUB).mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", message="below the paper settings")
    t0 = time.time()
    cohorts = ML.load_cohorts(args.subsample, args.seed)
    (power if args.stage == "power" else score)(args, out, quick, cohorts)
    say(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
