"""A306: do chronic conditions change the gain (log z2 ratio) or the departure (log psi ratio) of the organ they affect?
GTEx, A233 matched procedure; readout = module minus companion on the same permutations.

Registered in paper/docs/registered_experiments.md (A306). Conditions: T2D (module A164 MET), hypertension (arteries,
heart), smoking (lung); companion A164 COMP (tibial nerve, both skins). Matching, anchors, permutations: A233
(`A231_age_coupling/scripts/matched_label.py`); only module and companion tissues are fitted.

  power  permutation nulls only (no observed statistic): per condition and readout the null SD of
         Delta = module pooled - companion pooled, MDE_80 = 2.49 SD, scored flag (MDE_80 <= log 2)
         -> output/{power,power_tissues,matching,balance_summary}.csv
  score  (only after the dated power note) observed per-tissue readouts, Delta, two-sided permutation p, Holm over the
         scored readouts; arm-prior sensitivity on the observed values -> output/score_{tests,tissues,sensitivity}.csv

    python compute_a306.py --stage power|score [--jobs 4] [--n-perm 1000] [--conditions t2d htn smoking]
Controlled access: labels, ages and times stay in memory; only aggregates are written.
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
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
HEART = ["Heart - Left Ventricle", "Heart - Atrial Appendage"]
ARTERY = ["Artery - Aorta", "Artery - Coronary", "Artery - Tibial"]
MODULES = {"t2d": list(A.MET), "htn": ARTERY + HEART, "smoking": ["Lung"]}
READOUTS = ("log_z2_ratio", "log_psi_ratio")
NAMES = {"log_z2_ratio": "gain", "log_psi_ratio": "departure"}
BOUND = float(np.log(2.0))


def design_of(name: str) -> ML.Design:
    pools = lambda q: (q[name] == 1, q[name] == 0)
    return ML.Design(name, pools, ML.match_ratio, min_arm=0, min_cases=10, balance=("age", "ischaemic_h", "bmi", "t2d", "male"),
                     modules={"module": MODULES[name], "companion": list(A.COMP)})


def stat4(ds, lab, t, host_prior="pooled") -> list:
    try:
        tab = ds.contrast(lab, targets=[t], n_perm=0, n_boot=0, blocks="none", host_prior=host_prior).table
        if t in tab.index and "log_z2_ratio" in tab.columns:
            return [float(tab.loc[t, k]) for k in READOUTS]
    except Exception:
        pass
    return [np.nan] * len(READOUTS)


def perm_chunk(ds, plan, seeds):
    """(replicates, tissues, 2) of log z2 ratio and log psi ratio, case label swapped within matched sets (A164 labels)."""
    limit_threads(1)
    donors = sorted({d for s in plan.values() for x in s for d in x})
    out = []
    for sd in seeds:
        u = dict(zip(donors, np.random.default_rng(sd).random(len(donors))))
        out.append([stat4(ds, A.labels(plan[t], key=u), t) for t in plan])
    return np.array(out)


def nulls_of(name, plans, dss, out, quick, n_perm, jobs, seed):
    n_perm = n_perm or (20 if quick else 1000)
    return {c: A.cached(out / ML.DRAW_SUB / f"null_{name}_{c}.npy", quick, perm_chunk, n_perm, jobs,
                        A.seed_of("A306", name, c, seed), f"null {name} {c}", dss[c], plans[c]) for c in plans}


def delta(arrs: dict, mod: dict, cmp_: dict, j: int, nulls: dict):
    """Module pooled minus companion pooled for readout j: weights from the nulls (tissue inverse variance within cohort,
    cohort inverse variance of the module null), applied to `arrs` (cohort -> (rows, tissues, readouts))."""
    cc = [c for c in mod if c in cmp_]
    sm, sc = {c: mod[c] for c in cc}, {c: cmp_[c] for c in cc}
    _, twm, cwm, _ = B.pooled_null(ML.sub(nulls, sm), j)
    _, twc, _, _ = B.pooled_null(ML.sub(nulls, sc), j)
    f = lambda tw, s: sum(cwm[c] * B.wmean(arrs[c][:, s[c], j], tw[c]) for c in cwm)
    return f(twm, sm) - f(twc, sc), f(twm, sm), f(twc, sc)


def build(name, cohorts, quick):
    d = design_of(name)
    tissues = d.modules["module"] + d.modules["companion"]
    plans, dss, anchors, pbs, M, Bt = ML.build(d, cohorts, quick=quick, tissues=tissues, tag="gtex_a306")
    return d, plans, dss, M, Bt


def power(args, out, quick, cohorts):
    prows, trows, mrows, brows = [], [], [], []
    for name in args.conditions:
        d, plans, dss, M, Bt = build(name, cohorts, quick)
        mrows.append(M); brows.append(Bt)
        say(f"{name}: " + "; ".join(f"{c} {list(p)}" for c, p in plans.items()))
        nulls = nulls_of(name, plans, dss, out, quick, args.n_perm, args.jobs, args.seed)
        mod, cmp_ = ML.selection(plans, d.modules["module"]), ML.selection(plans, d.modules["companion"])
        for j, ro in enumerate(READOUTS):
            if not (mod and cmp_):
                prows.append(dict(condition=name, readout=NAMES[ro], note="no module or companion tissue"))
                continue
            nd, nm, nc = delta(nulls, mod, cmp_, j, nulls)
            sd = float(np.nanstd(nd, ddof=1))
            prows.append(dict(condition=name, readout=NAMES[ro], cohorts=";".join(c for c in mod if c in cmp_),
                              module_tissues=sum(len(v) for v in mod.values()),
                              companion_tissues=sum(len(v) for v in cmp_.values()),
                              module_sets=int(M[M.eligible & M.tissue.isin(d.modules["module"])].n_sets.sum()),
                              null_mean=float(np.nanmean(nd)), null_sd=sd, mde_80=ML.MDE_K * sd,
                              null_sd_module=float(np.nanstd(nm, ddof=1)), null_sd_companion=float(np.nanstd(nc, ddof=1)),
                              finite=float(np.isfinite(nd).mean()), scored=ML.MDE_K * sd <= BOUND))
            for c in nulls:
                for k, t in enumerate(plans[c]):
                    trows.append(dict(condition=name, readout=NAMES[ro], cohort=c, tissue=t,
                                      n_sets=len(plans[c][t]), null_sd=float(np.nanstd(nulls[c][:, k, j], ddof=1)),
                                      finite=float(np.isfinite(nulls[c][:, k, j]).mean())))
        pd.DataFrame(prows).to_csv(out / "power.csv", index=False)
        pd.DataFrame(trows).to_csv(out / "power_tissues.csv", index=False)
    M, Bt = pd.concat(mrows, ignore_index=True), pd.concat(brows, ignore_index=True)
    M.to_csv(out / "matching.csv", index=False)
    ML.balance_summary(Bt, M).to_csv(out / "balance_summary.csv", index=False)
    pd.set_option("display.width", 250)
    say("power\n" + pd.DataFrame(prows).round(4).to_string(index=False))


def score(args, out, quick, cohorts):
    P = pd.read_csv(out / "power.csv")
    scored = {(r.condition, r.readout) for r in P.itertuples() if r.scored}
    tests, tis, sens = [], [], []
    for name in args.conditions:
        d, plans, dss, M, _ = build(name, cohorts, quick)
        nulls = nulls_of(name, plans, dss, out, quick, args.n_perm, args.jobs, args.seed)
        mod, cmp_ = ML.selection(plans, d.modules["module"]), ML.selection(plans, d.modules["companion"])
        obs = {c: np.array([[stat4(dss[c], A.labels(plans[c][t]), t) for t in plans[c]]]) for c in plans}
        arm = {c: np.array([[stat4(dss[c], A.labels(plans[c][t]), t, "arm") for t in plans[c]]]) for c in plans}
        for c in plans:
            for k, t in enumerate(plans[c]):
                tis.append(dict(condition=name, cohort=c, tissue=t, n_sets=len(plans[c][t]),
                                role="module" if t in d.modules["module"] else "companion",
                                log_z2_ratio=obs[c][0, k, 0], log_psi_ratio=obs[c][0, k, 1],
                                log_z2_ratio_arm=arm[c][0, k, 0], log_psi_ratio_arm=arm[c][0, k, 1]))
        for j, ro in enumerate(READOUTS):
            nd, _, _ = delta(nulls, mod, cmp_, j, nulls)
            o, om, oc = (float(np.asarray(x).item()) for x in delta(obs, mod, cmp_, j, nulls))
            a, _, _ = delta(arm, mod, cmp_, j, nulls)
            tests.append(dict(condition=name, readout=NAMES[ro], delta=o, module=om, companion=oc,
                              p=A.p2(o, nd), null_sd=float(np.nanstd(nd, ddof=1)), mde_80=ML.MDE_K * float(np.nanstd(nd, ddof=1)),
                              scored=(name, NAMES[ro]) in scored))
            sens.append(dict(condition=name, readout=NAMES[ro], delta_arm_prior=float(np.asarray(a).item()),
                             p_arm_prior=A.p2(float(np.asarray(a).item()), nd)))
        say(f"{name} scored")
    T = pd.DataFrame(tests)
    m = T.scored
    if m.any():
        T.loc[m, "p_holm"] = ML.holm(T.p[m].to_numpy())
    T["verdict"] = np.where(~T.scored, "descriptive (not powered)",
                            np.where(T.p_holm < 0.05, "detected", "not detected"))
    T.to_csv(out / "score_tests.csv", index=False)
    pd.DataFrame(tis).to_csv(out / "score_tissues.csv", index=False)
    pd.DataFrame(sens).to_csv(out / "score_sensitivity.csv", index=False)
    pd.set_option("display.width", 250)
    say("tests\n" + T.round(4).to_string(index=False) + "\n" + pd.DataFrame(sens).round(4).to_string(index=False))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--conditions", nargs="*", default=list(MODULES), choices=list(MODULES))
    args = ap.parse_args(argv)
    quick = cli.setup(args)["quick"]
    out = OUT / ("subsample" if quick else ".")
    (out / ML.DRAW_SUB).mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", message="below the paper settings")
    t0 = time.time()
    cohorts = ML.load_cohorts(args.subsample, args.seed)
    (power if args.stage == "power" else score)(args, out, quick, cohorts)
    say(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
