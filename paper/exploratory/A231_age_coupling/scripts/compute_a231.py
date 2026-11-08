"""A231: do organ coupling or organ phase relations change with age? GTEx, A208 procedure with age as the label.

Registered in paper/docs/registered_experiments.md (A231). Old (60-70 y) = case, young (21-49 y) = control; 50-59 y unused.
Two stages:
  power  permutation null only (no observed statistic is computed or printed): null SD and MDE_80 (2.49 x null SD) of the
         pooled d* and offset per cohort and overall, MDE of Q (SD s in h of per-tissue offset shifts detected with 80 % at
         the null q95), matching and balance -> output/{power,power_tissues,matching,balance,balance_summary}.csv.
  score  (run only after the power note is dated) observed d* / offset per tissue, matched-set bootstrap, pooled tests
         H231.1-3, companions (3 lowest-coupling tissues, 50 pseudo-clock sets for offset and Q), D231.4 timed cosinor
         comparator (fast), D231.5 fast donors who died 08-18 h only, D231.6 host spread -> output/score_*.csv.
         `--placebo`: the case label is a random within-set draw (code test, writes under output/placebo/).

    python compute_a231.py --stage power [--subsample N] [--tissues T ..] [--jobs 8] [--n-perm 1000]
Controlled access: exact age and time of death stay in memory; only aggregates are written.

Choices not fixed by the registration
- Old = 60 <= age <= 70, young = 21 <= age <= 49 (exact AGE, integers 20-70).
- Matching `matched_label.match_pairs`: per tissue among donors with data in it; the smaller arm is processed first, fewest
  candidates first (ties by id); candidates = same sex and Hardy class, nearest ischaemic time (missing = tissue median).
- Eligible tissue: >= 20 donors per arm among donors with data in the tissue (before matching) and >= 10 matched pairs; no
  tissue excluded for being sex-specific (matching is exact on sex).
- Estimation and pooling exactly as A208 (imported); permutations: one uniform key per donor, shared across tissues.
- Q: all eligible tissues of both cohorts in one sum (also reported per cohort); w_t = 1 / null variance of the tissue offset
  from the 1,000 permutations (as the pooling weights); non-finite tissue draws dropped and Q rescaled to the full weight.
  MDE_Q: s found by bisection with one fixed N(0,1) matrix on the 1,000 null draws, power = share of shifted draws above the
  null q95 of Q (the registered 'SD of independent per-tissue shifts').
- Holm over H231.1-3 uses p of d* (two-sided), of the pooled offset and of Q.
- D231.4: per arm, Dataset.timed_loadings on the arm's matched donors (all clock genes, GTEx waveform, ML own noise, hours =
  time of death); tissue phase = clock_offset_h, minus the arm's circular mean over tissues; old minus young per tissue,
  pooled with the clocksync offset weights; SE from the same matched-set bootstrap (replicates resampled jointly over
  tissues); agreement = Jammalamadaka-SenGupta circular correlation and Pearson r of the per-tissue differences.
- D231.5: fast donors with 8 <= time of death < 18 h, then everything re-matched, re-permuted and re-bootstrapped (eligibility
  as above on that subset); no pseudo-clock companion; agreement with the full-day values = per-tissue correlation over
  common tissues and the difference of pooled estimates.
- D231.6: per cohort, median over eligible tissues of the host spread V of the old and young arm (contrast table V_case,
  V_control) and of log(V_old / V_young).
"""

from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matched_label as ML  # noqa: E402
from matched_label import A, B, S, cs, say  # noqa: E402
from clocksync import stats as CS  # noqa: E402
from common import cli  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
W = 2 * np.pi / 24
OLD, YOUNG = (60, 70), (21, 49)


def pools(q):
    return (q.age >= OLD[0]) & (q.age <= OLD[1]), (q.age >= YOUNG[0]) & (q.age <= YOUNG[1])


DESIGN = ML.Design("age", pools, ML.match_pairs, min_arm=20, min_cases=10, balance=("ischaemic_h", "bmi", "t2d", "male"))


# ---- D231.4 timed comparator -------------------------------------------------------------------------------------------
def timed_diff(frames: dict, tod: pd.Series, plan_arms: dict) -> np.ndarray:
    """Per tissue, old minus young tissue phase (h) relative to each arm's circular mean over tissues.
    frames: tissue -> pseudobulk rows; plan_arms: tissue -> (old donor ids, young donor ids) as in `frames`."""
    tissues = list(plan_arms)
    rel = []
    for a in (0, 1):
        rows = pd.concat([frames[t][frames[t].donor.isin(set(plan_arms[t][a]))] for t in tissues]).reset_index(drop=True)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tab = cs.Dataset(rows, anchors=None, name="timed").timed_loadings(tod)
            ph = tab.clock_offset_h.reindex(tissues).to_numpy(float)
        except Exception:
            ph = np.full(len(tissues), np.nan)
        ok = np.isfinite(ph)
        m = CS.circmean(ph[ok]) if ok.any() else np.nan
        rel.append(CS.wrap(ph - m))
    return np.asarray(CS.wrap(rel[0] - rel[1]), float)


def timed_obs(pb, tod, plan):
    frames = {t: pb[pb.compartment == t] for t in plan}
    return timed_diff(frames, tod, {t: ([s[0] for s in plan[t]], [s[1] for s in plan[t]]) for t in plan})


def timed_boot_chunk(pb, tod, plan, seeds):
    """Matched-set bootstrap of `timed_obs` (cases resampled jointly over tissues, as `compute_a208.boot_chunk`)."""
    A.limit_threads(1)
    frames = {t: pb[pb.compartment == t] for t in plan}
    cases = sorted({s[0] for v in plan.values() for s in v})
    by_case = {t: {s[0]: s for s in plan[t]} for t in plan}
    out = []
    for sd in seeds:
        draw = np.random.default_rng(sd).choice(len(cases), len(cases))
        fr, arms, hours = {}, {}, {}
        for t in plan:
            pairs = [(by_case[t][cases[i]], j) for j, i in enumerate(draw) if cases[i] in by_case[t]]
            parts = []
            for a in (0, 1):
                m = pd.DataFrame([(s[a], f"{s[a]}@{j}") for s, j in pairs], columns=["orig", "copy"])
                f = frames[t].merge(m, left_on="donor", right_on="orig")
                f["donor"] = f["copy"]
                parts.append(f.drop(columns=["orig", "copy"]))
                hours.update({c: tod[o] for o, c in zip(m["orig"], m["copy"])})
            fr[t] = pd.concat(parts)
            arms[t] = (list(parts[0].donor), list(parts[1].donor))
        out.append(timed_diff(fr, pd.Series(hours), arms))
    return np.array(out)


# ---- one pass of H231.1-3 on a cohort set ------------------------------------------------------------------------------
def readouts(args, out, quick, cohorts, suffix, placebo, pseudo: bool):
    """Estimates, tests and timed comparator on `cohorts` (full: vent + fast; D231.5: fast, 08-18 h). Returns a dict of tables."""
    plans, dss, anchors, pbs, M, Bt = ML.build(DESIGN, cohorts, quick=quick, tissues=args.tissues, tag="gtex_age")
    if placebo:
        plans = ML.placebo_plans(plans, args.seed + 7)
    nulls = ML.nulls_of("age", plans, dss, out, quick, args.n_perm, args.jobs, args.seed, suffix)
    boots = ML.boots_of("age", pbs, anchors, plans, out, quick, args.n_boot, args.jobs, args.seed, suffix)
    obs, T = ML.observed(dss, plans)
    for c in plans:
        T.loc[T.cohort == c, "null_sd_d_star"] = [float(np.nanstd(nulls[c][:, k, 0], ddof=1)) for k, t in enumerate(plans[c])
                                                   if t in set(T[T.cohort == c].tissue)]
    sel = ML.selection(plans, lambda t: True)
    P = ML.pooled_tests(obs, nulls, boots, sel, "all", "age")
    qt = ML.q_test(obs, nulls, sel)
    res = dict(plans=plans, M=M, T=T, P=P, Q=qt, obs=obs, nulls=nulls, boots=boots, sel=sel, pbs=pbs)
    # D231.4: timed comparator, fast cohort only (needs time of death)
    if "fast" in plans:
        tod = cohorts["fast"][2].tod
        tpath = out / ML.DRAW_SUB / f"timed_age{suffix}_fast"
        t_obs = timed_obs(pbs["fast"], tod, plans["fast"])
        t_boot = A.cached(Path(str(tpath) + "_boot.npy"), quick, timed_boot_chunk, args.n_boot or (10 if quick else 200),
                          args.jobs, A.seed_of("A231timed", suffix, args.seed), "timed boot", pbs["fast"], tod, plans["fast"])
        res.update(timed_obs=t_obs, timed_boot=t_boot)
    return res


def timed_table(res, label):
    """Timed pooled offset (fast), its bootstrap SE, clocksync's value and agreement per tissue."""
    if "timed_obs" not in res:
        return pd.DataFrame()
    j = list(res["plans"]).index("fast")
    tw = B.tissue_weights(res["nulls"]["fast"][:, :, 1])
    obs_t, boot_t = res["timed_obs"], res["timed_boot"]
    est_t = float(B.wmean(obs_t[None], tw)[0])
    se_t = float(np.nanstd(B.wmean(boot_t, tw), ddof=1))
    o_cs = res["obs"]["fast"][:, 1]
    est_c = float(B.wmean(o_cs[None], tw)[0])
    pr = res["P"][(res["P"].readout == "offset_h") & (res["P"].set == "fast")].iloc[0]
    ok = np.isfinite(obs_t) & np.isfinite(o_cs)
    return pd.DataFrame([dict(sample=label, n_tissues=int(ok.sum()), timed_offset_h=est_t, timed_boot_se=se_t,
                              clocksync_offset_h=est_c, clocksync_boot_se=pr.boot_se, se_ratio_timed_over_clocksync=se_t / pr.boot_se,
                              circ_corr_js=float(CS.js_corr(obs_t[ok] * W, o_cs[ok] * W)),
                              pearson_r=float(np.corrcoef(obs_t[ok], o_cs[ok])[0, 1]))])


def spread_table(res, label):
    """D231.6: host spread V of the old and young arm, median over tissues per cohort."""
    T = res["T"]
    rows = []
    for c, g in T.groupby("cohort"):
        rows.append(dict(sample=label, cohort=c, n_tissues=len(g), median_V_old=g.V_case.median(),
                         median_V_young=g.V_control.median(), median_log_V_old_over_young=float(np.log(g.V_case / g.V_control).median()),
                         share_tissues_old_lower=float((g.V_case < g.V_control).mean())))
    return pd.DataFrame(rows)


# ---- stages ------------------------------------------------------------------------------------------------------------
def power(args, out, quick):
    cohorts = ML.load_cohorts(args.subsample, args.seed)
    plans, dss, anchors, pbs, M, Bt = ML.build(DESIGN, cohorts, quick=quick, tissues=args.tissues, tag="gtex_age")
    M.to_csv(out / "matching.csv", index=False)
    Bt.to_csv(out / "balance.csv", index=False)
    ML.balance_summary(Bt, M).to_csv(out / "balance_summary.csv", index=False)
    ML.eligibility_summary(M).to_csv(out / "eligible_summary.csv", index=False)
    for c in plans:
        e = M[(M.cohort == c) & M.eligible]
        say(f"{c}: {len(plans[c])} eligible tissues, {int(e.n_sets.sum())} tissue pairs, "
            f"{sum(len(q) for q in [cohorts[c][1]])} donors with a known age")
    nulls = ML.nulls_of("age", plans, dss, out, quick, args.n_perm, args.jobs, args.seed)
    sel = ML.selection(plans, lambda t: True)
    P = pd.DataFrame(ML.power_rows(nulls, plans, M, sel, "all", "age", args.seed))
    P.to_csv(out / "power.csv", index=False)
    pd.DataFrame(ML.tissue_rows(nulls, plans, M, "age")).to_csv(out / "power_tissues.csv", index=False)
    say("power\n" + P.round(4).to_string(index=False))


def score(args, out, quick):
    n_pseudo = 5 if quick else 50
    cohorts = ML.load_cohorts(args.subsample, args.seed)
    res = readouts(args, out, quick, cohorts, "", args.placebo, True)
    plans, obs, nulls, sel, T, P, Q = (res[k] for k in ("plans", "obs", "nulls", "sel", "T", "P", "Q"))
    # companion for d*: the 3 lowest-coupling tissues per cohort (by the control = young arm), same permutations and weights
    _, tw, cw, _ = B.pooled_null(ML.sub(nulls, sel), 0)
    low = {c: np.argsort(T[T.cohort == c].c_control.to_numpy())[:3] for c in cw}
    lw = {c: np.isin(np.arange(len(plans[c])), low[c]) * tw[c] for c in cw}
    lw = {c: v / v.sum() for c, v in lw.items()}
    nl = sum(cw[c] * B.wmean(nulls[c][:, :, 0], lw[c]) for c in cw)
    Dl, _ = B.pool_estimate({c: obs[c][:, 0] for c in cw}, lw, cw)
    comp = dict(companion="d_star_3_lowest_young_coupling", estimate=Dl, null_sd=float(np.nanstd(nl, ddof=1)), p_perm=A.p2(Dl, nl),
                tissues=";".join(f"{c}:{t}" for c in cw for t in np.array(list(plans[c]))[low[c]]))
    ps = ML.pseudo_clock(plans, res_anchors(res), n_pseudo, args.jobs, args.seed)
    spec, off_v, q_v = ML.pseudo_summary(ps, obs, nulls, sel, sel)
    pd.DataFrame(dict(offset_h=off_v, Q=q_v)).to_csv(out / "score_pseudo.csv", index=False)
    P.to_csv(out / "score_pooled.csv", index=False)
    T.to_csv(out / "score_tissues.csv", index=False)
    pd.DataFrame([comp]).to_csv(out / "score_companion_d_star.csv", index=False)
    pp = P[P.set == "pooled"].set_index("readout")
    p3 = [pp.loc["d_star", "p_perm"], pp.loc["offset_h", "p_perm"], Q["p_perm"]]
    adj = ML.holm(p3)
    tests = pd.DataFrame([
        dict(hypothesis="H231.1", readout="d_star", estimate=pp.loc["d_star", "estimate"], lo90=pp.loc["d_star", "lo90"],
             hi90=pp.loc["d_star", "hi90"], mde_80=pp.loc["d_star", "mde_80"], p_perm=p3[0], p_holm=adj[0],
             verdict=(("pass (equivalent within +-0.15)" if pp.loc["d_star", "equivalent"] else "not shown")
                      if pp.loc["d_star", "mde_80"] <= ML.D_BOUND else "not powered (descriptive only)"),
             companion_3_lowest_d_star=Dl, companion_p=comp["p_perm"]),
        dict(hypothesis="H231.2", readout="offset_h", estimate=pp.loc["offset_h", "estimate"], lo90=pp.loc["offset_h", "lo90"],
             hi90=pp.loc["offset_h", "hi90"], mde_80=pp.loc["offset_h", "mde_80"], p_perm=p3[1], p_holm=adj[1],
             equivalent_within_1h=bool(pp.loc["offset_h", "equivalent"]), pseudo_abs_q95=spec["offset_pseudo_abs_q95"],
             clock_specific=spec["offset_clock_specific"]),
        dict(hypothesis="H231.3", readout="Q", estimate=Q["Q"], p_perm=p3[2], p_holm=adj[2], null_q95=Q["null_q95"],
             null_mean=Q["null_mean"], pseudo_q95=spec["Q_pseudo_q95"], clock_specific=spec["Q_clock_specific"])])
    tests.to_csv(out / "score_tests.csv", index=False)
    # D231.4 / D231.6 on the full-day samples; D231.5 on fast donors who died 08-18 h
    tt = [timed_table(res, "fast_full")]
    sp = [spread_table(res, "full")]
    wh = ML.load_cohorts(args.subsample, args.seed, hours=(8, 18))
    wh = {"fast": wh["fast"]}
    r5 = readouts(args, out, quick, wh, "_wh", args.placebo, False)
    tt.append(timed_table(r5, "fast_08_18h"))
    sp.append(spread_table(r5, "fast_08_18h"))
    pd.concat(tt, ignore_index=True).to_csv(out / "score_d231_4_timed.csv", index=False)
    pd.concat(sp, ignore_index=True).to_csv(out / "score_d231_6_spread.csv", index=False)
    rows = []
    for _, r in r5["P"].iterrows():
        f = P[(P.readout == r.readout) & (P.set == "fast")].iloc[0]
        rows.append(dict(readout=r.readout, set=r.set, estimate_wh=r.estimate, se_wh=r.se_used, estimate_full=f.estimate,
                         se_full=f.se_used, difference=r.estimate - f.estimate, mde_80_wh=r.mde_80))
    qw = r5["Q"]
    rows.append(dict(readout="Q", set="fast", estimate_wh=qw["Q"], p_wh=qw["p_perm"], null_q95_wh=qw["null_q95"]))
    # per-tissue agreement of the clocksync offsets between 08-18 h and full-day fast donors
    ta, tb = T[T.cohort == "fast"].set_index("tissue").offset_h, r5["T"].set_index("tissue").offset_h
    ix = ta.index.intersection(tb.index)
    if len(ix) > 2:
        rows.append(dict(readout="per_tissue_offset_agreement", set="fast", n_tissues=len(ix),
                         pearson_r=float(np.corrcoef(ta[ix], tb[ix])[0, 1]),
                         circ_corr_js=float(CS.js_corr(ta[ix].to_numpy() * W, tb[ix].to_numpy() * W))))
    pd.DataFrame(rows).to_csv(out / "score_d231_5_workinghours.csv", index=False)
    r5["P"].to_csv(out / "score_d231_5_pooled.csv", index=False)
    say("tests\n" + tests.round(4).to_string(index=False))


def res_anchors(res):
    return {c: [a for a in S.TISSUES if a in set(pb.compartment)] for c, pb in res["pbs"].items()}


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--tissues", nargs="*", default=None, help="restrict to these tissues (check)")
    ap.add_argument("--placebo", action="store_true", help="score stage code test on random within-set labels")
    args = ap.parse_args(argv)
    quick = cli.setup(args)["quick"] or bool(args.tissues)
    out = OUT / ("subsample" if quick else ".")
    if args.placebo:
        out = OUT / "placebo"
    (out / ML.DRAW_SUB).mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", message="below the paper settings")
    t0 = time.time()
    (power if args.stage == "power" else score)(args, out, quick)
    say(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
