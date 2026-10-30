"""A164: do T2D and obesity leave the GTEx metabolic organs coupled to the shared clock? Matched design, power, equivalence.

Registered in paper/docs/registered_experiments.md (A164). Choices not in the registration: ../README.md.
Order: (1) matching, balance and power (written before any test statistic is computed or printed); (2) H164.1-3 and
sensitivities. Controlled access: only aggregates are printed or written; matched sets stay in memory.

    python compute_a164.py [--subsample N] [--seed S] [--jobs 4] [--n-perm 2000] [--n-boot 500]
Outputs: ../output/a164_{matching,balance,power,tissues,tests,sensitivity,a47_side_by_side}.csv (+ raw null / bootstrap
arrays of module statistics, no donor information, in ../output/draws/)
"""

from __future__ import annotations

import sys
import time
import warnings
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_3/scripts/gtex_cues"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _covariates as K  # noqa: E402
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

MET = ["Liver", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Pancreas", "Muscle - Skeletal"]
COMP = ["Nerve - Tibial", "Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]
TISS = MET + COMP
COH = {"vent": (0,), "fast": (1, 2)}
BOUND, Z90, MDE_K = 0.15, 1.645, 2.49
RATIO = {"t2d": 3, "obese": 1}
OUT = HERE.parent / "output"
PRIMARY = ("t2d", "obese")
SENS = ("t2d_bmi", "t2d_nosteroid", "obese_nosteroid")


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def seed_of(*parts) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode()) % (2 ** 31)


# ---- matching (blind to clock data) ----------------------------------------------------------------------------------
def pools(q: pd.DataFrame, design: str):
    """(cases, controls, ratio, bmi_caliper) donor index of a design among the donors `q` (age, bmi, t2d, steroid ...)."""
    base = design.split("_")[0]
    ok = pd.Series(True, index=q.index)
    if design.endswith("nosteroid"):
        ok &= q.steroid != 1
    if base == "t2d":
        cas, con = ok & (q.t2d == 1), ok & (q.t2d == 0)
        if design == "t2d_bmi":
            cas, con = cas & q.bmi.notna(), con & q.bmi.notna()
    else:
        ok &= q.t2d != 1
        cas, con = ok & (q.bmi >= 30), ok & (q.bmi >= 18.5) & (q.bmi < 25)
    return cas, con, RATIO[base], (3.0 if design == "t2d_bmi" else None)


def match(q: pd.DataFrame, design: str):
    """Sets [case, controls...]: exact on sex and Hardy class, age caliper 5 y, nearest age then ischaemic time, bmi
    caliper for t2d_bmi; up to `ratio` controls per case, without replacement, assigned in rounds (round r gives every
    case, fewest candidates first, its r-th control). Returns (sets, n_cases_in_pool, n_controls_in_pool)."""
    cas, con, ratio, bcal = pools(q, design)
    cases, ctrls = sorted(q.index[cas]), q.index[con]
    cand = {}
    for c in cases:
        r = q.loc[c]
        m = ctrls[(q.loc[ctrls, "male"].to_numpy() == r.male) & (q.loc[ctrls, "hardy"].to_numpy() == r.hardy)
                  & (np.abs(q.loc[ctrls, "age"].to_numpy() - r.age) <= 5)]
        if bcal is not None:
            m = m[np.abs(q.loc[m, "bmi"].to_numpy() - r.bmi) <= bcal]
        da, di = np.abs(q.loc[m, "age"].to_numpy() - r.age), np.abs(q.loc[m, "ischaemic_h"].to_numpy() - r.ischaemic_h)
        o = np.lexsort((np.array(m, dtype=str), di, da)) if len(m) else []
        cand[c] = list(np.asarray(m)[o])
    order = sorted(cases, key=lambda c: (len(cand[c]), c))
    used, got = set(), {c: [] for c in cases}
    for _ in range(ratio):
        for c in order:
            nxt = next((x for x in cand[c] if x not in used), None)
            if nxt is not None:
                used.add(nxt)
                got[c].append(nxt)
    sets = [[c] + got[c] for c in cases if got[c]]
    return sets, len(cases), len(ctrls)


def smd(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / s) if s > 0 else 0.0


def plan_of(q: pd.DataFrame, tissue_donors: dict, design: str, min_cases: int):
    """Per tissue matched sets (donors with data in the tissue), the matching and balance rows."""
    plan, mrows, brows = {}, [], []
    for t in TISS:
        qt = q.loc[q.index.isin(tissue_donors[t])]
        sets, ncase, nctrl = match(qt, design)
        cas = [s[0] for s in sets]
        ctr = [x for s in sets for x in s[1:]]
        elig = len(cas) >= min_cases
        mrows.append(dict(design=design, tissue=t, role="module" if t in MET else "companion", n_cases_pool=ncase,
                          n_controls_pool=nctrl, n_matched_cases=len(cas), n_unmatched_cases=ncase - len(cas),
                          n_controls_used=len(ctr), mean_controls_per_set=len(ctr) / len(cas) if cas else np.nan,
                          eligible=elig))
        pc, pk, *_ = pools(qt, design)
        for v in ("age", "ischaemic_h", "bmi"):
            brows.append(dict(design=design, tissue=t, variable=v,
                              smd_before=smd(qt.loc[pc, v], qt.loc[pk, v]),
                              smd_after=smd(qt.loc[cas, v], qt.loc[ctr, v]) if cas else np.nan))
        if elig:
            plan[t] = sets
    return plan, mrows, brows


# ---- estimation --------------------------------------------------------------------------------------------------------
def labels(sets, key=None, copy=None):
    """donor -> case / control; `key` (donor -> random number): in every set the member with the largest key is the
    case (label swap within sets, shared across tissues through the donor keys); `copy`: donor id suffix."""
    lab = {}
    for s in sets:
        k = 0 if key is None else int(np.argmax([key[d] for d in s]))
        for i, d in enumerate(s):
            lab[d if copy is None else f"{d}@{copy}"] = "case" if i == k else "control"
    return pd.Series(lab)


def dstar(ds, lab, t) -> float:
    try:
        tab = ds.contrast(lab, targets=[t], n_perm=0, n_boot=0, blocks="none").table
        return float(tab.loc[t, "d_star"]) if t in tab.index and "d_star" in tab.columns else np.nan
    except Exception:
        return np.nan


def stats_of(vals: dict) -> tuple:
    """(module mean, module minus companion) of finite per-tissue values."""
    m = [v for t, v in vals.items() if t in MET and np.isfinite(v)]
    c = [v for t, v in vals.items() if t in COMP and np.isfinite(v)]
    return (float(np.mean(m)) if m else np.nan, float(np.mean(m) - np.mean(c)) if m and c else np.nan)


def perm_chunk(ds, plan, seeds):
    limit_threads(1)
    donors = sorted({d for s in plan.values() for x in s for d in x})
    out = []
    for sd in seeds:
        u = dict(zip(donors, np.random.default_rng(sd).random(len(donors))))
        out.append([dstar(ds, labels(plan[t], key=u), t) for t in plan])
    return np.array(out)


def boot_chunk(pb, anchors, plan, seeds):
    """Matched-set bootstrap: cases resampled with replacement (jointly over tissues), each draw brings its set(s) in
    every tissue; donors are copied (id@draw), the Dataset and its host are refitted on the replicate."""
    limit_threads(1)
    cases = sorted({s[0] for v in plan.values() for s in v})
    out = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        draw = rng.choice(len(cases), len(cases))
        by_case = {t: {s[0]: s for s in plan[t]} for t in plan}
        sets = {t: [] for t in plan}
        mp = set()
        for j, i in enumerate(draw):
            for t in plan:
                s = by_case[t].get(cases[i])
                if s is not None:
                    sets[t].append((s, j))
                    mp.update((d, f"{d}@{j}") for d in s)
        m = pd.DataFrame(sorted(mp), columns=["orig", "copy"])
        rep = pb.merge(m, left_on="donor", right_on="orig")
        rep["donor"] = rep["copy"]
        rep = rep.drop(columns=["orig", "copy"])
        ds = cs.Dataset(rep.reset_index(drop=True), anchors=[a for a in anchors if a in set(rep.compartment)])
        row = []
        for t in plan:
            lab = pd.concat([labels([s], copy=j) for s, j in sets[t]])
            row.append(dstar(ds, lab, t))
        out.append(row)
    return np.array(out)


def run_chunks(fn, n, jobs, seed0, label, *args):
    """Parallel chunks of replicate seeds; logs after every chunk. Returns array (n, n_tissues)."""
    seeds = [seed0 + i for i in range(n)]
    size = max(1, int(np.ceil(n / (jobs * 6))))
    ch = [seeds[i:i + size] for i in range(0, n, size)]
    res, t0 = [], time.time()
    gen = Parallel(n_jobs=jobs, return_as="generator")(delayed(fn)(*args, c) for c in ch)
    for k, r in enumerate(gen, 1):
        res.append(r)
        say(f"{label}: chunk {k}/{len(ch)} done, {time.time() - t0:.0f} s")
    return np.vstack(res)


def cached(path, quick, fn, n, jobs, seed0, label, *args):
    """Replicate array of `run_chunks`, kept on disk for resuming a full run."""
    if not quick and path.exists():
        return np.load(path)
    a = run_chunks(fn, n, jobs, seed0, label, *args)
    if not quick:
        np.save(path, a)
    return a


def modules(arr, tissues):
    """(n, 2) array of module statistics per replicate from an (n, n_tissues) array."""
    return np.array([stats_of(dict(zip(tissues, r))) for r in arr])


def pool(vals, ses):
    """Inverse-variance pooling over cohorts (finite values and SEs): (estimate, SE, weights)."""
    ok = [c for c in vals if np.isfinite(vals[c]) and np.isfinite(ses[c]) and ses[c] > 0]
    if not ok:
        return np.nan, np.nan, {}
    w = {c: 1 / ses[c] ** 2 for c in ok}
    tot = sum(w.values())
    return sum(w[c] * vals[c] for c in ok) / tot, float(np.sqrt(1 / tot)), {c: w[c] / tot for c in ok}


def p2(obs, null):
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    return float((1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (1 + len(null)))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_perm = args.n_perm or (20 if quick else 2000)
    n_boot = args.n_boot or (10 if quick else 500)
    min_cases = 3 if quick else 10
    out = OUT / "subsample" if quick else OUT
    (out / "draws").mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", message="below the paper settings")
    ph = K.protected()

    # ---- data and matching -------------------------------------------------------------------------------------------
    pbs, dss, plans, mrows, brows = {}, {}, {}, [], []
    for c, h in COH.items():
        pb, don = K.cohort(h, subsample=args.subsample, seed=args.seed)
        keep = set(S.TISSUES) | set(TISS)
        pb = pb[pb.compartment.isin(keep)].reset_index(drop=True)
        q = don[["hardy", "male", "ischaemic_h"]].join(ph[["age", "bmi", "t2d", "steroid"]]).dropna(
            subset=["age", "ischaemic_h"])
        td = {t: set(pb.donor[pb.compartment == t]) for t in TISS}
        pbs[c] = pb
        dss[c] = cs.Dataset(pb, anchors=[a for a in S.TISSUES if a in set(pb.compartment)], name=f"gtex_{c}")
        for dsg in PRIMARY + SENS:
            plans[c, dsg], m, b = plan_of(q, td, dsg, min_cases)
            mrows += [dict(cohort=c, **r) for r in m]
            brows += [dict(cohort=c, **r) for r in b]
    M, B = pd.DataFrame(mrows), pd.DataFrame(brows)
    M.to_csv(out / "a164_matching.csv", index=False)
    B.to_csv(out / "a164_balance.csv", index=False)
    pd.set_option("display.width", 220)
    say("matching written (counts)\n", M[M.design.isin(PRIMARY)][["cohort", "design", "tissue", "n_cases_pool",
         "n_matched_cases", "n_unmatched_cases", "n_controls_used", "eligible"]].to_string(index=False))
    bp = B[B.design.isin(PRIMARY)].groupby(["design", "variable"]).smd_after.agg(lambda x: np.nanmax(np.abs(x)))
    say("max |SMD| after matching over cohorts and tissues\n", bp.round(3).to_string())

    # ---- null and bootstrap of the primary designs (no observed statistic yet), then power ---------------------------------
    raw_n, raw_b, nul, bt = {}, {}, {}, {}
    for c in COH:
        for dsg in PRIMARY:
            plan = plans[c, dsg]
            if not any(t in MET for t in plan):
                continue
            anchors = [a for a in S.TISSUES if a in set(pbs[c].compartment)]
            sd0 = seed_of("A164", c, dsg, args.seed)
            raw_n[c, dsg] = cached(out / "draws" / f"null_{c}_{dsg}.npy", quick, perm_chunk, n_perm, args.jobs, sd0,
                                   f"null {c} {dsg}", dss[c], plan)
            raw_b[c, dsg] = cached(out / "draws" / f"boot_{c}_{dsg}.npy", quick, boot_chunk, n_boot, args.jobs,
                                   sd0 + 10 ** 6, f"boot {c} {dsg}", pbs[c], anchors, plan)
            nul[c, dsg] = modules(raw_n[c, dsg], list(plan))
            bt[c, dsg] = modules(raw_b[c, dsg], list(plan))
    weights, prow = {}, []
    for dsg in PRIMARY:
        for j, stat in enumerate(("module", "module_minus_companion")):
            se = {c: float(np.nanstd(bt[c, dsg][:, j], ddof=1)) if (c, dsg) in bt else np.nan for c in COH}
            sd = {c: float(np.nanstd(nul[c, dsg][:, j], ddof=1)) if (c, dsg) in nul else np.nan for c in COH}
            _, se_p, w = pool({c: 0.0 for c in COH}, se)
            weights[dsg, stat] = w
            nl = sum(w[c] * nul[c, dsg][:, j] for c in w) if w else np.array([np.nan])
            for name, s_, n_ in [(c, se[c], sd[c]) for c in COH] + [("pooled", se_p, float(np.nanstd(nl, ddof=1)))]:
                prow.append(dict(design=dsg, statistic=stat, set=name, boot_se=s_, null_sd=n_, mde_80=MDE_K * n_,
                                 expected_halfwidth_90=Z90 * s_, bound=BOUND,
                                 powered=bool(Z90 * s_ <= BOUND) if np.isfinite(s_) else False,
                                 weights=";".join(f"{c}:{x:.3f}" for c, x in w.items()) if name == "pooled" else ""))
    Pw = pd.DataFrame(prow)
    Pw.to_csv(out / "a164_power.csv", index=False)
    say("power written\n", Pw[Pw.statistic == "module"].round(4).drop(columns="weights").to_string(index=False))

    # ---- observed estimates and tests ------------------------------------------------------------------------------------
    def observed(c, dsg):
        plan = plans[c, dsg]
        res = {}
        for t in plan:
            tab = dss[c].contrast(labels(plan[t]), targets=[t], n_perm=0, n_boot=0, blocks="none").table
            res[t] = tab.loc[t] if t in tab.index else None
        return plan, res

    obs, trows = {}, []
    for c in COH:
        for dsg in PRIMARY:
            if (c, dsg) not in raw_n:
                continue
            plan, res = observed(c, dsg)
            obs[c, dsg] = {t: float(r["d_star"]) if r is not None and "d_star" in r and np.isfinite(r["d_star"]) else np.nan
                           for t, r in res.items()}
            for t, r in res.items():
                if r is None or "d_star" not in r:
                    continue
                k = list(plan).index(t)
                nl, bb = raw_n[c, dsg][:, k], raw_b[c, dsg][:, k]
                trows.append(dict(cohort=c, design=dsg, tissue=t, role="module" if t in MET else "companion",
                                  n_sets=len(plan[t]), n_case=r.get("donors_case"), n_control=r.get("donors_control"),
                                  c_control=r.get("c_control"), c_case=r.get("c_case"), d=r.get("d"), d_star=r["d_star"],
                                  boot_se=float(np.nanstd(bb, ddof=1)), null_sd=float(np.nanstd(nl, ddof=1)),
                                  p_perm=p2(r["d_star"], nl)))
    pd.DataFrame(trows).to_csv(out / "a164_tissues.csv", index=False)

    tests, cohort_mod = [], {}
    for dsg in PRIMARY:
        for j, stat in enumerate(("module", "module_minus_companion")):
            est = {c: stats_of(obs[c, dsg])[j] if (c, dsg) in nul else np.nan for c in COH}
            se = {c: float(np.nanstd(bt[c, dsg][:, j], ddof=1)) if (c, dsg) in bt else np.nan for c in COH}
            D, se_p, w = pool(est, se)
            nl = sum(w[c] * nul[c, dsg][:, j] for c in w)
            tests.append(dict(design=dsg, statistic=stat, set="pooled", D=D, se=se_p, lo90=D - Z90 * se_p,
                              hi90=D + Z90 * se_p, p_perm=p2(D, nl), null_sd=float(np.nanstd(nl, ddof=1)),
                              mde_80=MDE_K * float(np.nanstd(nl, ddof=1)), n_perm=int(np.isfinite(nl).sum()),
                              n_boot=n_boot, weights=";".join(f"{c}:{x:.3f}" for c, x in w.items())))
            for c in COH:
                if np.isfinite(est[c]):
                    nlc = nul[c, dsg][:, j]
                    tests.append(dict(design=dsg, statistic=stat, set=c, D=est[c], se=se[c],
                                      lo90=est[c] - Z90 * se[c], hi90=est[c] + Z90 * se[c], p_perm=p2(est[c], nlc),
                                      null_sd=float(np.nanstd(nlc, ddof=1)), mde_80=MDE_K * float(np.nanstd(nlc, ddof=1)),
                                      n_perm=int(np.isfinite(nlc).sum()), n_boot=n_boot))
    T = pd.DataFrame(tests)
    pw = Pw.set_index(["design", "statistic", "set"]).powered
    pooled = T[T.set == "pooled"].set_index(["design", "statistic"])
    hol = pooled.loc[[(d, "module") for d in PRIMARY], "p_perm"].to_numpy()
    o = np.argsort(hol)
    adj = np.empty(2)
    adj[o] = np.minimum(1, np.maximum.accumulate((2 - np.arange(2)) * hol[o]))
    hyp = []
    for i, (dsg, name) in enumerate((("t2d", "H164.1"), ("obese", "H164.2"))):
        r = pooled.loc[(dsg, "module")]
        powered = bool(pw[(dsg, "module", "pooled")])
        inside = bool(r.lo90 > -BOUND and r.hi90 < BOUND)
        hyp.append(dict(hypothesis=name, design=dsg, D=r.D, lo90=r.lo90, hi90=r.hi90, bound=BOUND, powered=powered,
                        equivalence=("equivalent" if inside else "not shown") if powered else "not powered (no verdict)"))
    hyp = pd.DataFrame(hyp)
    h3 = []
    for i, dsg in enumerate(PRIMARY):
        pm = pooled.loc[(dsg, "module_minus_companion"), "p_perm"]
        h3.append(dict(hypothesis="H164.3", design=dsg, D=pooled.loc[(dsg, "module"), "D"], p_perm=hol[i], p_holm=adj[i],
                       effect=bool(adj[i] < .05), p_module_minus_companion=pm,
                       organ_specific=bool(adj[i] < .05 and pm < .05)))
    h3 = pd.DataFrame(h3)
    pd.concat([hyp, h3], ignore_index=True).to_csv(out / "a164_tests_hypotheses.csv", index=False)
    T.to_csv(out / "a164_tests.csv", index=False)

    # ---- sensitivities (descriptive, point estimates, weights of the primary design of the same family) -----------------------
    srows = []
    fam = {"t2d_bmi": "t2d", "t2d_nosteroid": "t2d", "obese_nosteroid": "obese"}
    for dsg in SENS:
        est = {}
        for c in COH:
            if not any(t in MET for t in plans[c, dsg]):
                continue
            plan, res = observed(c, dsg)
            vals = {t: float(r["d_star"]) if r is not None and "d_star" in r else np.nan for t, r in res.items()}
            est[c] = vals
            for t, v in vals.items():
                srows.append(dict(design=dsg, set=c, tissue=t, statistic="tissue", D=v, n_sets=len(plan[t])))
            for j, stat in enumerate(("module", "module_minus_companion")):
                srows.append(dict(design=dsg, set=c, tissue="", statistic=stat, D=stats_of(vals)[j],
                                  n_sets=sum(len(v) for t, v in plan.items() if t in MET)))
        for j, stat in enumerate(("module", "module_minus_companion")):
            w = weights.get((fam[dsg], stat), {})
            v = {c: stats_of(est[c])[j] for c in est}
            ok = [c for c in v if np.isfinite(v[c]) and c in w]
            if ok:
                srows.append(dict(design=dsg, set="pooled", tissue="", statistic=stat,
                                  D=sum(w[c] * v[c] for c in ok) / sum(w[c] for c in ok), n_sets=np.nan))
    pd.DataFrame(srows).to_csv(out / "a164_sensitivity.csv", index=False)

    # ---- the A47 unmatched estimates side by side (archive-scored tables of figure_3; loss = control minus case) ----------------
    a47 = PAPER / "results/figure_3/output/gtex_cues"
    side = []
    try:
        t47 = pd.read_csv(a47 / "a47_t2d_tissues.csv")
        side += [dict(source="A47 unmatched", design="t2d", set=r.cohort, tissue=r.tissue, D=-r.loss_star)
                 for r in t47.itertuples()]
        p47 = pd.read_csv(a47 / "a47_t2d_pooled.csv")
        side += [dict(source="A47 unmatched", design="t2d", set="pooled", tissue=f"{r.quantity} | {r.module}", D=-r.D)
                 for r in p47.itertuples() if r.quantity == "dstar"]
        b47 = pd.read_csv(a47 / "a47_primary.csv")
        side += [dict(source="A47 BMI slope per SD (not comparable in units)", design="bmi", set="pooled",
                      tissue="module minus companion", D=r.statistic) for r in b47.itertuples()]
    except FileNotFoundError:
        say("A47 tables not found")
    for r in T.itertuples():
        side.append(dict(source="A164 matched", design=r.design, set=r.set, tissue=r.statistic, D=r.D))
    pd.DataFrame(side).to_csv(out / "a164_a47_side_by_side.csv", index=False)

    say("tests\n", pd.concat([hyp, h3], ignore_index=True).round(4).to_string(index=False))
    say("pooled\n", T[T.set == "pooled"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
