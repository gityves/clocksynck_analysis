"""Matched-contrast machinery shared by A231 (age) and A233 (chronic conditions): the A208 / A164 procedure with any
binary donor label, any matching rule, a tissue set, plus the reordering statistic Q.

  label + matching   `Design`: `pools(q)` -> (case mask, control mask) among the donors of a tissue; `match(q, cas, con)` ->
                     sets [case, control, ...] (blind to clock data). `match_pairs` (A231: 1:1, exact sex and Hardy class,
                     nearest ischaemic time) and `match_ratio` (A233, the A164 T2D rule: up to 3 controls, exact sex and
                     Hardy class, age caliper 5 y, nearest age then ischaemic time).
  estimation         A208 functions by import (`stat2`, `perm_chunk`, `boot_chunk`, pooling, `pseudo_sets`); one
                     `Dataset.contrast` per tissue; permutations swap the case label within matched sets, shared across tissues.
  Q                  Q = sum_t w_t * offset_diff_t^2, w_t = 1 / null variance of the tissue offset (all tissues of all
                     cohorts in one sum). Null: the same permutations. MDE_Q: the SD s (h) of independent N(0, s^2) shifts
                     added to every tissue's null draw such that 80 % of the shifted draws exceed the null q95.
Controlled access: ages, conditions, times of death live in memory only; every file written is an aggregate.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER / "exploratory/A208_sex_coupling/scripts"))
import compute_a208 as B  # noqa: E402  (also puts paper/, A164 and gtex_cues on the path)
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

A, S, K = B.A, B.S, B.A.K
COH, say, MDE_K, Z90 = B.COH, B.say, B.MDE_K, B.Z90
COMP = A.COMP  # companion tissues of A164 (nerve, the two skins)
D_BOUND, OFF_BOUND = B.D_BOUND, B.OFF_BOUND
DRAW_SUB = "draws"


# ---- protected labels (in memory only) ---------------------------------------------------------------------------------
def protected() -> pd.DataFrame:
    """age (exact), bmi, t2d and the five A233 conditions (1 case, 0 control, NaN unknown) by donor."""
    ph = K.protected()[["age", "bmi", "t2d"]]
    p = gtex.phenotypes(["MHHTN", "MHHRTDIS", "MHCOPD", "MHDPRSSN", "MHSMKSTS"])
    ok = lambda x: pd.to_numeric(x, errors="coerce").where(lambda v: v.isin([0, 1]))
    c = pd.DataFrame({"htn": ok(p.MHHTN), "ihd": ok(p.MHHRTDIS), "copd": ok(p.MHCOPD), "depression": ok(p.MHDPRSSN),
                      "smoking": p.MHSMKSTS.map({"Yes": 1.0, "No": 0.0})})
    return ph.join(c)


# ---- matching ----------------------------------------------------------------------------------------------------------
def match_pairs(q: pd.DataFrame, cas: pd.Series, con: pd.Series) -> list:
    """1:1 sets [case, control]: exact on sex and Hardy class, nearest ischaemic time (missing = median). The smaller arm is
    processed first, fewest candidates first (ties by id); each takes the closest unused member of the other arm."""
    q = q.assign(ischaemic_h=q.ischaemic_h.fillna(q.ischaemic_h.median()))
    a, b = sorted(q.index[cas]), sorted(q.index[con])
    case_first = len(a) <= len(b)
    first, other = (a, b) if case_first else (b, a)
    cand = {}
    for x in first:
        r = q.loc[x]
        m = [y for y in other if q.male[y] == r.male and q.hardy[y] == r.hardy]
        cand[x] = sorted(m, key=lambda y: (abs(q.ischaemic_h[y] - r.ischaemic_h), y))
    used, sets = set(), []
    for x in sorted(first, key=lambda x: (len(cand[x]), x)):
        y = next((z for z in cand[x] if z not in used), None)
        if y is not None:
            used.add(y)
            sets.append([x, y] if case_first else [y, x])
    return sets


def match_ratio(q: pd.DataFrame, cas: pd.Series, con: pd.Series, ratio: int = 3, caliper: float = 5.0) -> list:
    """A164 rule (`compute_a164.match`): up to `ratio` controls per case without replacement, exact on sex and Hardy
    class, age caliper, nearest age then ischaemic time, controls assigned in rounds (fewest candidates first)."""
    q = q.dropna(subset=["age", "ischaemic_h"])
    cases, ctrls = sorted(q.index[cas.reindex(q.index).fillna(False)]), q.index[con.reindex(q.index).fillna(False)]
    cand = {}
    for c in cases:
        r = q.loc[c]
        m = ctrls[(q.loc[ctrls, "male"].to_numpy() == r.male) & (q.loc[ctrls, "hardy"].to_numpy() == r.hardy)
                  & (np.abs(q.loc[ctrls, "age"].to_numpy() - r.age) <= caliper)]
        da, di = np.abs(q.loc[m, "age"].to_numpy() - r.age), np.abs(q.loc[m, "ischaemic_h"].to_numpy() - r.ischaemic_h)
        o = np.lexsort((np.array(m, dtype=str), di, da)) if len(m) else []
        cand[c] = list(np.asarray(m)[o])
    used, got = set(), {c: [] for c in cases}
    for _ in range(ratio):
        for c in sorted(cases, key=lambda c: (len(cand[c]), c)):
            nxt = next((x for x in cand[c] if x not in used), None)
            if nxt is not None:
                used.add(nxt)
                got[c].append(nxt)
    return [[c] + got[c] for c in cases if got[c]]


@dataclass
class Design:
    name: str
    pools: callable            # q -> (case mask, control mask)
    match: callable            # (q, cas, con) -> sets
    min_arm: int = 0           # donors per arm with data in the tissue
    min_cases: int = 10        # matched cases (sets) per tissue
    balance: tuple = ("ischaemic_h", "bmi", "t2d", "male")
    modules: dict = field(default_factory=dict)  # name -> tissue names or callable(tissue) -> bool


def smd_rows(design, q, qt, cas, con, sets, tissue, cohort) -> list:
    case_d = [s[0] for s in sets]
    ctl_d = [x for s in sets for x in s[1:]]
    return [dict(cohort=cohort, design=design.name, tissue=t_, variable=v,
                 smd_before=A.smd(qt.loc[cas, v], qt.loc[con, v]),
                 smd_after=A.smd(q.loc[case_d, v], q.loc[ctl_d, v]) if sets else np.nan)
            for t_ in [tissue] for v in design.balance]


def make_plan(design: Design, q: pd.DataFrame, pb: pd.DataFrame, cohort: str, tissues=None, quick=False):
    """({tissue: sets}, matching rows, balance rows) of one cohort and design. Eligible: `min_arm` donors per arm with data
    in the tissue and `min_cases` matched cases (quick checks: 4 and 3)."""
    min_arm, min_cases = (4, 3) if quick else (design.min_arm, design.min_cases)
    plan, rows, brows = {}, [], []
    for t in sorted(set(pb.compartment)):
        if tissues and t not in tissues:
            continue
        qt = q.loc[q.index.isin(set(pb.donor[pb.compartment == t]))]
        cas, con = design.pools(qt)
        n1, n0 = int(cas.sum()), int(con.sum())
        if min(n1, n0) < max(min_arm, 1):
            if n1 + n0:
                rows.append(dict(cohort=cohort, design=design.name, tissue=t, n_case_pool=n1, n_control_pool=n0,
                                 n_sets=0, n_controls_used=0, eligible=False))
            continue
        sets = design.match(qt, cas, con)
        ok = len(sets) >= min_cases
        rows.append(dict(cohort=cohort, design=design.name, tissue=t, n_case_pool=n1, n_control_pool=n0, n_sets=len(sets),
                         n_controls_used=sum(len(s) - 1 for s in sets), eligible=ok))
        brows += smd_rows(design, q, qt, cas, con, sets, t, cohort)
        if ok:
            plan[t] = sets
    return plan, rows, brows


def eligibility_summary(M: pd.DataFrame) -> pd.DataFrame:
    """Per design and cohort: tissues with data, eligible tissues, matched cases (sets) over eligible tissues, controls per set."""
    rows = []
    for (d, c), g in M.groupby(["design", "cohort"]):
        e = g[g.eligible]
        rows.append(dict(design=d, cohort=c, n_tissues_with_both_arms=len(g), n_eligible_tissues=len(e),
                         n_sets_total=int(e.n_sets.sum()), min_sets=e.n_sets.min() if len(e) else np.nan,
                         median_sets=e.n_sets.median() if len(e) else np.nan, max_sets=e.n_sets.max() if len(e) else np.nan,
                         max_case_pool=int(g.n_case_pool.max()), max_control_pool=int(g.n_control_pool.max()),
                         mean_controls_per_set=e.n_controls_used.sum() / e.n_sets.sum() if len(e) else np.nan))
    return pd.DataFrame(rows)


def balance_summary(Bt: pd.DataFrame, M: pd.DataFrame) -> pd.DataFrame:
    """Per cohort, design and variable over eligible tissues: median and max |SMD| before and after matching."""
    e = Bt.merge(M[["cohort", "design", "tissue", "eligible"]], on=["cohort", "design", "tissue"])
    e = e[e.eligible]
    g = e.assign(a=e.smd_after.abs(), b=e.smd_before.abs()).groupby(["design", "cohort", "variable"])
    return g.agg(n_tissues=("a", "size"), median_abs_smd_before=("b", "median"), median_abs_smd_after=("a", "median"),
                 max_abs_smd_after=("a", "max")).reset_index()


# ---- data --------------------------------------------------------------------------------------------------------------
def load_cohorts(subsample=None, seed=0, hours=None) -> dict:
    """{cohort: (pseudobulk, donor table q with covariates, donor table with tod (controlled))}. `hours`: (lo, hi) keeps fast
    donors with lo <= time of death < hi."""
    ph = protected()
    out = {}
    for c, h in COH.items():
        pb, don = gtex.load(tissues=None, hardy=h, genes="clock", depth="library_size", require_tod=(c == "fast"))
        if hours is not None and c == "fast":
            keep = don.index[(don.tod >= hours[0]) & (don.tod < hours[1])]
            pb, don = pb[pb.donor.isin(keep)].reset_index(drop=True), don.loc[keep]
        if subsample:
            keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
            pb, don = pb[pb.donor.isin(keep)].reset_index(drop=True), don.loc[don.index.isin(keep)]
        q = don[["hardy", "male", "ischaemic_h"]].join(ph).dropna(subset=["age"])
        out[c] = (pb, q, don)
    return out


def build(design: Design, cohorts: dict, tissues=None, quick=False, tag="gtex", anchor_set=None):
    """plans, datasets, anchors, pseudobulks (plan tissues + anchors), matching and balance tables of one design.
    anchor_set: host anchor tissues (default: the A164 anchors `S.TISSUES`)."""
    anchor_set = S.TISSUES if anchor_set is None else tuple(anchor_set)
    plans, dss, anchors, pbs, mrows, brows = {}, {}, {}, {}, [], []
    for c, (pb, q, _) in cohorts.items():
        plan, m, b = make_plan(design, q, pb, c, tissues, quick)
        mrows += m
        brows += b
        if not plan:
            continue
        pbc = pb[pb.compartment.isin(set(plan) | set(anchor_set))].reset_index(drop=True)
        anchors[c] = [a for a in anchor_set if a in set(pbc.compartment)]
        plans[c], pbs[c] = plan, pbc
        dss[c] = cs.Dataset(pbc, anchors=anchors[c], name=f"{tag}_{design.name}_{c}")
    return plans, dss, anchors, pbs, pd.DataFrame(mrows), pd.DataFrame(brows)


def nulls_of(design_name, plans, dss, out, quick, n_perm, jobs, seed, suffix=""):
    """(replicates, tissues, 2) permutation null per cohort: d*, offset_diff_h. Cached in out/draws."""
    n_perm = n_perm or (20 if quick else 1000)
    return {c: A.cached(out / DRAW_SUB / f"null_{design_name}{suffix}_{c}.npy", quick, B.perm_chunk, n_perm, jobs,
                        A.seed_of("A231A233", design_name, suffix, c, seed), f"null {design_name}{suffix} {c}", dss[c],
                        plans[c]) for c in plans}


# ---- selections of tissues ---------------------------------------------------------------------------------------------
def selection(plans: dict, names) -> dict:
    """{cohort: column indices of the plan tissues that are in `names` (list or predicate)}, empty cohorts dropped."""
    f = names if callable(names) else (lambda t: t in set(names))
    sel = {c: [i for i, t in enumerate(p) if f(t)] for c, p in plans.items()}
    return {c: v for c, v in sel.items() if v}


def sub(arrs: dict, sel: dict) -> dict:
    return {c: arrs[c][:, sel[c], :] for c in sel}


# ---- Q -----------------------------------------------------------------------------------------------------------------
def q_weights(X: np.ndarray) -> np.ndarray:
    """1 / null variance of each column of an offset null (replicates x tissues); 0 where not estimable."""
    v = np.nanvar(X, axis=0, ddof=1)
    return np.where(np.isfinite(v) & (v > 0), 1 / v, 0.0)


def q_of(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Q = sum_t w_t x_t^2 per row; non-finite tissues dropped and the sum rescaled to the full weight."""
    ok = np.isfinite(X)
    num = (np.where(ok, X, 0.0) ** 2 * w).sum(1)
    den = (ok * w).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num * w.sum() / den


def q_mde(X: np.ndarray, w: np.ndarray, seed: int = 0, level: float = 0.8, alpha: float = 0.05):
    """(s, q95, null mean of Q): the SD s (h) of independent per-tissue N(0, s^2) shifts added to the null draws X for which
    a share `level` of the shifted draws exceeds the null q95 of Q (bisection on a fixed set of standard normals)."""
    q0 = q_of(X, w)
    q95 = float(np.nanquantile(q0, 1 - alpha))
    eps = np.random.default_rng(seed).standard_normal(X.shape)
    pw = lambda s: float(np.nanmean(q_of(X + s * eps, w) > q95))
    lo, hi = 0.0, 1.0
    while pw(hi) < level and hi < 1e3:
        hi *= 2
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pw(mid) < level else (lo, mid)
    return (lo + hi) / 2, q95, float(np.nanmean(q0))


def q_block(nulls: dict, seed: int = 0):
    """Offset-null columns of all cohorts side by side, their weights, and a cohort index per column."""
    cols = [nulls[c][:, :, 1] for c in nulls]
    X = np.hstack(cols)
    return X, q_weights(X), np.concatenate([[c] * a.shape[1] for c, a in zip(nulls, cols)])


# ---- power tables ------------------------------------------------------------------------------------------------------
def power_rows(nulls: dict, plans: dict, M: pd.DataFrame, sel: dict, label: str, design: str, seed: int = 0) -> list:
    """Rows (readouts d*, offset, Q) for the tissues `sel` (cohort -> column indices): per cohort and pooled."""
    sn = sub(nulls, sel)
    rows = []
    n_case = {c: int(M[(M.cohort == c) & M.eligible & M.tissue.isin([list(plans[c])[i] for i in sel[c]])].n_sets.sum())
              for c in sel}
    for j, stat in enumerate(B.STATS):
        null, tw, cw, coh = B.pooled_null(sn, j)
        bound = D_BOUND if j == 0 else OFF_BOUND
        for c in sn:
            sd = float(np.nanstd(coh[c], ddof=1))
            rows.append(dict(design=design, selection=label, readout=stat, set=c, n_tissues=len(sel[c]), n_sets=n_case[c],
                             null_sd=sd, mde_80=MDE_K * sd, weight=cw.get(c, np.nan), bound=bound,
                             powered_vs_bound=MDE_K * sd <= bound))
        sd = float(np.nanstd(null, ddof=1))
        rows.append(dict(design=design, selection=label, readout=stat, set="pooled", n_tissues=sum(map(len, sel.values())),
                         n_sets=sum(n_case.values()), null_sd=sd, mde_80=MDE_K * sd, weight=1.0, bound=bound,
                         powered_vs_bound=MDE_K * sd <= bound))
    Xs = {c: sn[c][:, :, 1] for c in sn}
    for name, parts in [(c, [c]) for c in sn] + [("pooled", list(sn))]:
        X = np.hstack([Xs[c] for c in parts])
        w = q_weights(X)
        s, q95, m = q_mde(X, w, seed)
        q0 = q_of(X, w)
        rows.append(dict(design=design, selection=label, readout="Q", set=name, n_tissues=X.shape[1],
                         n_sets=sum(n_case[c] for c in parts), null_sd=float(np.nanstd(q0, ddof=1)), mde_80=s,
                         q_null_mean=m, q_null_q95=q95, bound=np.nan, powered_vs_bound=np.nan))
    return rows


def tissue_rows(nulls: dict, plans: dict, M: pd.DataFrame, design: str) -> list:
    rows = []
    for j, stat in enumerate(B.STATS):
        for c in nulls:
            tw = B.tissue_weights(nulls[c][:, :, j])
            for k, t in enumerate(plans[c]):
                rows.append(dict(design=design, cohort=c, tissue=t, readout=stat,
                                 n_sets=int(M[(M.cohort == c) & (M.tissue == t)].n_sets.iloc[0]),
                                 null_sd=float(np.nanstd(nulls[c][:, k, j], ddof=1)), weight=float(tw[k])))
    return rows


# ---- score helpers (not run before the power note is dated) -------------------------------------------------------------
def placebo_plans(plans: dict, seed: int) -> dict:
    """Plans whose case label is a random within-set draw (shared across tissues): a test of the score code that never
    touches the real labels."""
    out = {}
    for c, plan in plans.items():
        donors = sorted({d for s in plan.values() for x in s for d in x})
        u = dict(zip(donors, np.random.default_rng(seed).random(len(donors))))
        out[c] = {t: [sorted(s, key=lambda d: -u[d]) for s in sets] for t, sets in plan.items()}
    return out


def observed(dss: dict, plans: dict) -> tuple:
    """Observed per-tissue contrast tables per cohort: ({cohort: array (tissues, 2) of d*, offset}, long table)."""
    obs, rows = {}, []
    for c in plans:
        a = np.full((len(plans[c]), 2), np.nan)
        for k, t in enumerate(plans[c]):
            try:
                r = dss[c].contrast(A.labels(plans[c][t]), targets=[t], n_perm=0, n_boot=0, blocks="none").table.loc[t]
            except Exception:
                continue
            a[k] = [r.get("d_star", np.nan), r.get("offset_diff_h", np.nan)]
            rows.append(dict(cohort=c, tissue=t, n_sets=len(plans[c][t]), c_control=r.get("c_control"),
                             c_case=r.get("c_case"), d=r.get("d"), d_star=a[k, 0], offset_h=a[k, 1],
                             V_control=r.get("V_control"), V_case=r.get("V_case"), z2_control=r.get("z2_control"),
                             z2_case=r.get("z2_case"), V0=r.get("V0")))
        obs[c] = a
    return obs, pd.DataFrame(rows)


def boot_chunk(pb, anchors, plan, seeds):
    """`compute_a208.boot_chunk` (cases resampled jointly over tissues, each brings its set(s) in every tissue, host refitted on
    the replicate) with a guard: a tissue none of the resampled cases has a set in gives NaN for that replicate."""
    limit_threads(1)
    cases = sorted({s[0] for v in plan.values() for s in v})
    by_case = {t: {s[0]: s for s in plan[t]} for t in plan}
    out = []
    for sd in seeds:
        draw = np.random.default_rng(sd).choice(len(cases), len(cases))
        sets, mp = {t: [] for t in plan}, set()
        for j, i in enumerate(draw):
            for t in plan:
                s = by_case[t].get(cases[i])
                if s is not None:
                    sets[t].append((s, j))
                    mp.update((d, f"{d}@{j}") for d in s)
        m = pd.DataFrame(sorted(mp), columns=["orig", "copy"])
        rep = pb.merge(m, left_on="donor", right_on="orig")
        rep["donor"] = rep["copy"]
        rep = rep.drop(columns=["orig", "copy"]).reset_index(drop=True)
        ds = cs.Dataset(rep, anchors=[a for a in anchors if a in set(rep.compartment)])
        out.append([B.stat2(ds, pd.concat([A.labels([s], copy=j) for s, j in sets[t]]), t) if sets[t] else (np.nan, np.nan)
                    for t in plan])
    return np.array(out)


def boots_of(name, pbs, anchors, plans, out, quick, n_boot, jobs, seed, suffix=""):
    n_boot = n_boot or (10 if quick else 200)
    return {c: A.cached(out / DRAW_SUB / f"boot_{name}{suffix}_{c}.npy", quick, boot_chunk, n_boot, jobs,
                        A.seed_of("A231A233", name, suffix, c, seed) + 10 ** 6, f"boot {name}{suffix} {c}", pbs[c],
                        anchors[c], plans[c]) for c in plans}


def pooled_tests(obs, nulls, boots, sel, label, design, case_label="case") -> pd.DataFrame:
    """Pooled estimates for the tissues `sel` as A208 `score`: estimate, bootstrap SE, null SD, A202 interval, MDE, p."""
    sn, sb = sub(nulls, sel), sub(boots, sel)
    so = {c: obs[c][sel[c]] for c in sel}
    rows = []
    for j, stat in enumerate(B.STATS):
        null, tw, cw, _ = B.pooled_null(sn, j)
        bound = D_BOUND if j == 0 else OFF_BOUND
        D, per = B.pool_estimate({c: so[c][:, j] for c in cw}, tw, cw)
        bt = sum(cw[c] * B.wmean(sb[c][:, :, j], tw[c]) for c in cw)
        for name, est, nl, bb in [(c, per[c], B.wmean(sn[c][:, :, j], tw[c]), B.wmean(sb[c][:, :, j], tw[c])) for c in cw] \
                + [("pooled", D, null, bt)]:
            se_b, sd_n = float(np.nanstd(bb, ddof=1)), float(np.nanstd(nl, ddof=1))
            se = max(se_b, sd_n)
            rows.append(dict(design=design, selection=label, readout=stat, set=name, estimate=est, boot_se=se_b,
                             null_sd=sd_n, se_used=se, lo90=est - Z90 * se, hi90=est + Z90 * se, bound=bound,
                             mde_80=MDE_K * sd_n, p_perm=A.p2(est, nl),
                             equivalent=bool(est - Z90 * se > -bound and est + Z90 * se < bound)))
    return pd.DataFrame(rows)


def q_test(obs: dict, nulls: dict, sel: dict) -> dict:
    """Observed Q over the tissues `sel` (all cohorts in one sum), one-sided add-one permutation p, weights from the null."""
    X = np.hstack([nulls[c][:, sel[c], 1] for c in sel])
    w = q_weights(X)
    xo = np.concatenate([obs[c][sel[c], 1] for c in sel])[None]
    q_obs = float(q_of(xo, w)[0])
    q0 = q_of(X, w)
    q0 = q0[np.isfinite(q0)]
    return dict(Q=q_obs, p_perm=float((1 + np.sum(q0 >= q_obs)) / (1 + len(q0))), null_q95=float(np.quantile(q0, .95)),
                null_mean=float(q0.mean()), w=w)


def pseudo_tissue_offsets(i, pseudo, plans, anchors) -> dict:
    """Per-tissue offset differences of pseudo-clock set i (whole pipeline on the pseudo genes), per cohort."""
    limit_threads(1)
    out = {}
    for c in plans:
        w = pseudo[c][i]
        pb, _ = gtex.load(tissues=list(plans[c]) + anchors[c], hardy=COH[c], genes=list(w.genes), depth="library_size",
                          require_tod=(c == "fast"))
        ds = cs.Dataset(pb, anchors=anchors[c], waveform=w)
        out[c] = np.array([B.stat2(ds, A.labels(plans[c][t]), t)[1] for t in plans[c]], float)
    return out


def pseudo_clock(plans, anchors, n_pseudo, jobs, seed):
    """n_pseudo sets of per-tissue offset differences per cohort: list of {cohort: array}. Observed values only."""
    pseudo = {c: B.pseudo_sets(COH[c], n_pseudo, A.seed_of("A231A233pseudo", c, seed)) for c in plans}
    gen = Parallel(n_jobs=jobs, return_as="generator")(delayed(pseudo_tissue_offsets)(i, pseudo, plans, anchors)
                                                          for i in range(n_pseudo))
    res = []
    for k, v in enumerate(gen, 1):
        say(f"pseudo set {k}/{n_pseudo}")
        res.append(v)
    return res


def pseudo_summary(pseudo, obs, nulls, sel_off, sel_q) -> dict:
    """Pooled offset and Q of every pseudo-clock set (same weights as the clock null) and the clock-specificity flags.
    sel_off: selection of the offset readout (pooled by the clock null weights); sel_q: selection of Q."""
    sn = sub(nulls, sel_off)
    _, tw, cw, _ = B.pooled_null(sn, 1)
    off = np.array([B.pool_estimate({c: p[c][sel_off[c]] for c in cw}, tw, cw)[0] for p in pseudo], float)
    qt = q_test(obs, nulls, sel_q)
    qv = np.array([float(q_of(np.concatenate([p[c][sel_q[c]] for c in sel_q])[None], qt["w"])[0]) for p in pseudo])
    Do = B.pool_estimate({c: obs[c][sel_off[c], 1] for c in cw}, tw, cw)[0]
    q95o, q95q = float(np.nanquantile(np.abs(off), .95)), float(np.nanquantile(qv, .95))
    return dict(offset_pseudo_abs_q95=q95o, offset_clock_specific=bool(abs(Do) > q95o), Q_pseudo_q95=q95q,
                Q_clock_specific=bool(qt["Q"] > q95q), n_sets=len(pseudo)), off, qv


def holm(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    adj[o] = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    return adj
