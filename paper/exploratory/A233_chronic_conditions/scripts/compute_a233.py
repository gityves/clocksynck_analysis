"""A233: do chronic conditions other than T2D and obesity change organ coupling or phase relations? GTEx, A164 procedure.

Registered in paper/docs/registered_experiments.md (A233). Conditions (case = 1, control = 0, unknown excluded): hypertension
(MHHTN), ischaemic heart disease (MHHRTDIS), COPD (MHCOPD), major depression (MHDPRSSN), smoking (MHSMKSTS Yes vs No). Matching
blind to clock data: A164 T2D rule (up to 3 controls, exact sex and Hardy class, age caliper 5 y, nearest age then ischaemic
time). Machinery shared with A231: ../../A231_age_coupling/scripts/matched_label.py. Two stages:
  power  permutation null only (no observed statistic): per condition MDE_80 of the body-wide pooled d* and offset, module d*
         and offset, MDE of Q (body-wide and module); matching, balance (SMD of age, ischaemic time, BMI, T2D; sex) ->
         output/{power,power_tissues,matching,balance,balance_summary,eligible_summary}.csv.
  score  (run only after the power note is dated) per condition H233.c.1-3 with companions; Holm over conditions per
         hypothesis -> output/score_*.csv. `--placebo`: random within-set case labels (code test, output/placebo/).

    python compute_a233.py --stage power [--subsample N] [--conditions htn ..] [--jobs 8] [--n-perm 1000]
Controlled access: conditions, exact age and ischaemic time stay in memory; only aggregates are written.

Choices not fixed by the registration
- All non-cultured tissues enter the body-wide readouts (no tissue excluded, sex-specific organs included: matching is exact on
  sex); a tissue enters with >= 10 matched cases in the cohort (>= 1 case and 1 control in the pool to be matched at all).
- Module d* and offset: inverse-null-variance pooling over the module tissues of a cohort, then over cohorts (A208), weights
  from the same permutations (the registration says "A164 / A208 estimation"). Module tissue names: `MODULES` below
  (depression: every "Brain - ..." tissue). Companion tissues of the module-minus-companion contrast: A164 COMP (tibial nerve and
  the two skins), pooled the same way.
- Body-wide Q: all eligible tissues of both cohorts in one sum (as A231); a module Q is also reported (power only).
- Permutation: the case label is swapped with a random member of its matched set (A164 `labels`), shared across tissues.
- Holm (score): per hypothesis (c.2 Q, c.3 module d*, c.3 module offset) over the five conditions.
"""

from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A231_age_coupling/scripts"))
import matched_label as ML  # noqa: E402
from matched_label import A, B, say  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
HEART = ["Heart - Left Ventricle", "Heart - Atrial Appendage"]
MODULES = {
    "htn": HEART + ["Artery - Aorta", "Artery - Coronary", "Artery - Tibial", "Kidney - Cortex"],
    "ihd": HEART + ["Artery - Coronary", "Artery - Aorta"],
    "copd": ["Lung"],
    "smoking": ["Lung"],
    "depression": lambda t: t.startswith("Brain - "),
}
COLUMN = {"htn": "htn", "ihd": "ihd", "copd": "copd", "depression": "depression", "smoking": "smoking"}


def design_of(name: str) -> ML.Design:
    col = COLUMN[name]
    pools = lambda q: (q[col] == 1, q[col] == 0)
    return ML.Design(name, pools, ML.match_ratio, min_arm=0, min_cases=10, balance=("age", "ischaemic_h", "bmi", "t2d", "male"),
                     modules={"module": MODULES[name], "companion": A.COMP})


def power(args, out, quick, cohorts):
    prows, trows, mrows, brows = [], [], [], []
    for name in args.conditions:
        d = design_of(name)
        plans, dss, anchors, pbs, M, Bt = ML.build(d, cohorts, quick=quick, tissues=args.tissues, tag="gtex_cond")
        mrows.append(M)
        brows.append(Bt)
        say(f"{name}: " + "; ".join(f"{c} {len(p)} eligible tissues, {int(M[(M.cohort == c) & M.eligible].n_sets.sum())} sets"
                                    for c, p in plans.items()))
        if not plans:
            continue
        nulls = ML.nulls_of(name, plans, dss, out, quick, args.n_perm, args.jobs, args.seed)
        for label, names in (("bodywide", lambda t: True), ("module", d.modules["module"])):
            sel = ML.selection(plans, names)
            if sel:
                prows += ML.power_rows(nulls, plans, M, sel, label, name, args.seed)
            else:
                prows.append(dict(design=name, selection=label, readout="none", set="no eligible module tissue"))
        trows += ML.tissue_rows(nulls, plans, M, name)
        pd.DataFrame(prows).to_csv(out / "power.csv", index=False)  # partial results survive a failure
        pd.DataFrame(trows).to_csv(out / "power_tissues.csv", index=False)
    M, Bt = pd.concat(mrows, ignore_index=True), pd.concat(brows, ignore_index=True)
    M.to_csv(out / "matching.csv", index=False)
    Bt.to_csv(out / "balance.csv", index=False)
    ML.balance_summary(Bt, M).to_csv(out / "balance_summary.csv", index=False)
    ML.eligibility_summary(M).to_csv(out / "eligible_summary.csv", index=False)
    P = pd.DataFrame(prows)
    P.to_csv(out / "power.csv", index=False)
    pd.set_option("display.width", 250)
    say("power\n" + P.drop(columns=["bound"], errors="ignore").round(4).to_string(index=False))


def score(args, out, quick, cohorts):
    n_pseudo = 5 if quick else 50
    rows, tt, comps = [], [], []
    keep = {}
    for name in args.conditions:
        d = design_of(name)
        plans, dss, anchors, pbs, M, _ = ML.build(d, cohorts, quick=quick, tissues=args.tissues, tag="gtex_cond")
        if not plans:
            continue
        if args.placebo:
            plans = ML.placebo_plans(plans, args.seed + 7)
        nulls = ML.nulls_of(name, plans, dss, out, quick, args.n_perm, args.jobs, args.seed)
        boots = ML.boots_of(name, pbs, anchors, plans, out, quick, args.n_boot, args.jobs, args.seed)
        obs, T = ML.observed(dss, plans)
        T.insert(0, "design", name)
        tt.append(T)
        all_sel = ML.selection(plans, lambda t: True)
        mod_sel = ML.selection(plans, d.modules["module"])
        cmp_sel = ML.selection(plans, d.modules["companion"])
        P = pd.concat([ML.pooled_tests(obs, nulls, boots, all_sel, "bodywide", name)]
                      + ([ML.pooled_tests(obs, nulls, boots, mod_sel, "module", name)] if mod_sel else []), ignore_index=True)
        rows.append(P)
        Q = ML.q_test(obs, nulls, all_sel)
        # module minus companion d* on the same permutations (organ-specific), per cohort pooled with the module weights
        mc = dict(design=name, module_minus_companion_d_star=np.nan, p_module_minus_companion=np.nan)
        if mod_sel and cmp_sel:
            cc = [c for c in mod_sel if c in cmp_sel]
            sm, sc = {c: mod_sel[c] for c in cc}, {c: cmp_sel[c] for c in cc}
            _, twm, cwm, _ = B.pooled_null(ML.sub(nulls, sm), 0)
            _, twc, _, _ = B.pooled_null(ML.sub(nulls, sc), 0)
            wm = lambda a, tw, s: sum(cwm[c] * B.wmean(a[c][:, s[c]] if a[c].ndim == 2 else a[c][s[c]][None], tw[c]) for c in cwm)
            nm = {c: nulls[c][:, :, 0] for c in cwm}
            on = {c: obs[c][:, 0] for c in cwm}
            null_diff = wm(nm, twm, sm) - wm(nm, twc, sc)
            o_diff = np.asarray(wm(on, twm, sm) - wm(on, twc, sc)).item()  # (1,) array -> scalar
            mc.update(module_minus_companion_d_star=o_diff, p_module_minus_companion=A.p2(o_diff, null_diff))
        comps.append(mc)
        ps = ML.pseudo_clock(plans, {c: anchors[c] for c in plans}, n_pseudo, args.jobs, args.seed)
        spec, off_v, q_v = ML.pseudo_summary(ps, obs, nulls, mod_sel or all_sel, all_sel)
        pd.DataFrame(dict(offset_h=off_v, Q=q_v)).to_csv(out / f"score_pseudo_{name}.csv", index=False)
        pm = P[(P.set == "pooled")].set_index(["selection", "readout"])
        keep[name] = dict(Q=Q, spec=spec, pm=pm, mc=mc, has_mod=bool(mod_sel))
        say(f"{name} scored")
    pd.concat(rows, ignore_index=True).to_csv(out / "score_pooled.csv", index=False)
    pd.concat(tt, ignore_index=True).to_csv(out / "score_tissues.csv", index=False)
    pd.DataFrame(comps).to_csv(out / "score_companion_module_minus_companion.csv", index=False)
    names = list(keep)
    hq = ML.holm([keep[n]["Q"]["p_perm"] for n in names])
    mods = [n for n in names if keep[n]["has_mod"]]
    hd = dict(zip(mods, ML.holm([keep[n]["pm"].loc[("module", "d_star"), "p_perm"] for n in mods]))) if mods else {}
    ho = dict(zip(mods, ML.holm([keep[n]["pm"].loc[("module", "offset_h"), "p_perm"] for n in mods]))) if mods else {}
    tests = []
    for i, n in enumerate(names):
        k, pm = keep[n], keep[n]["pm"]
        bd = pm.loc[("bodywide", "d_star")]
        r = dict(condition=n, bodywide_d_star=bd.estimate, bodywide_lo90=bd.lo90, bodywide_hi90=bd.hi90,
                 bodywide_mde_80=bd.mde_80, H233_1=(("equivalent" if bd.equivalent else "not shown") if bd.mde_80 <= ML.D_BOUND
                                                    else "not powered (descriptive)"),
                 Q=k["Q"]["Q"], Q_p=k["Q"]["p_perm"], Q_p_holm=hq[i], Q_null_q95=k["Q"]["null_q95"],
                 Q_pseudo_q95=k["spec"]["Q_pseudo_q95"], Q_clock_specific=k["spec"]["Q_clock_specific"])
        if k["has_mod"]:
            md, mo = pm.loc[("module", "d_star")], pm.loc[("module", "offset_h")]
            r.update(module_d_star=md.estimate, module_d_star_p=md.p_perm, module_d_star_p_holm=hd[n],
                     module_minus_companion_d_star=k["mc"]["module_minus_companion_d_star"],
                     module_minus_companion_p=k["mc"]["p_module_minus_companion"],
                     module_offset_h=mo.estimate, module_offset_p=mo.p_perm, module_offset_p_holm=ho[n],
                     offset_pseudo_abs_q95=k["spec"]["offset_pseudo_abs_q95"],
                     offset_clock_specific=k["spec"]["offset_clock_specific"])
        tests.append(r)
    T = pd.DataFrame(tests)
    T.to_csv(out / "score_tests.csv", index=False)
    say("tests\n" + T.round(4).to_string(index=False))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--tissues", nargs="*", default=None, help="restrict to these tissues (check)")
    ap.add_argument("--conditions", nargs="*", default=list(COLUMN), choices=list(COLUMN))
    ap.add_argument("--placebo", action="store_true", help="score stage code test on random within-set labels")
    args = ap.parse_args(argv)
    quick = cli.setup(args)["quick"] or bool(args.tissues)
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
