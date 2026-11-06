"""A208: do women and men differ in organ coupling or offsets? GTEx, A164 procedure with sex as the label (female = case).

Registered in paper/docs/registered_experiments.md (A208). Two stages:
  power  permutation null only (no observed statistic is computed or printed): null SD of the pooled d* and of the pooled
         offset difference per cohort and overall, MDE_80 = 2.49 x null SD -> output/power.csv (+ matching, per-tissue null SD).
  score  observed d* and offset per tissue, matched-set bootstrap, pooled estimates, equivalence (+-0.15 on d*, +-1 h on the
         offset), permutation p, specificity companions (3 lowest-coupling tissues; pseudo-clock offsets) -> output/score_*.csv.
         Run only after the power note is dated in the registration.

    python compute_a208.py --stage power [--subsample N] [--tissues T1 T2 ..] [--jobs 8] [--n-perm 1000]
Controlled access: sex, Hardy class, age bracket and ischaemic time are the public GTEx annotations (`common.gtex.load`);
no protected phenotype is read. Only aggregates are written; matched sets stay in memory.

Choices not fixed by the registration
- Cohorts as registered: vent = Hardy 0; fast = Hardy 1-2 with time of death (not used, only the filter).
- Tissues: every non-cultured GTEx tissue with >= 20 female and >= 20 male donors (donors with data in the tissue) in the
  cohort; sex-specific organs excluded (Ovary, Uterus, Vagina, Cervix x2, Fallopian Tube, Prostate, Testis, Breast).
  A tissue is also dropped when fewer than 10 matched pairs remain (min-pairs, as A164 "eligible"). A tissue may enter one
  cohort and not the other.
- Matching per tissue, among donors with data in that tissue and a known age bracket: exact on Hardy class and age bracket
  (the bracket midpoint `age` of `common.gtex.load`, so exact equality = same bracket), 1:1 without replacement; women are
  processed fewest-candidates first (ties by donor id), each takes the man with the closest ischaemic time (ties by id;
  a missing ischaemic time counts as the cohort median). Unmatched donors are dropped.
- d* and offset_diff_h come from one `Dataset.contrast` call per tissue (levels control = male, case = female; targets=[t],
  blocks="none", n_perm=0, n_boot=0) on the matched donors of that tissue; host = the 15 GTEx anchors, as A164.
- Permutations (1,000): one uniform key per donor and replicate, the set member with the largest key is the case (A164
  `labels`), shared across tissues; the same draws feed d*, offset and the companions.
- Pooling: within a cohort over tissues, then over cohorts, inverse-variance weights 1 / null_sd^2 (null SD of each
  tissue's d* or offset; of each cohort's pooled value), fixed by the permutation null (not by observed data) and reused in
  `score`, so that the MDE is that of the statistic scored. Offsets are pooled as a weighted arithmetic mean of the
  per-tissue circular differences (hours, |value| << 12). A non-finite tissue in a replicate is dropped from that
  replicate and weights renormalised.
- score: bootstrap as A164 (women resampled jointly over tissues, each brings its pair in every tissue, host refitted);
  pooled 90 % interval = estimate +- 1.645 x max(pooled bootstrap SE, pooled null SD) (A202 rule, as
  `Contrast.equivalence`). Two-sided add-one permutation p. Companion 1: the 3 tissues of lowest mean (c_male, c_female)
  per cohort. Companion 2: 50 expression-matched pseudo-clock sets (`clocksync.paired.pseudo_waveforms`, match on the mean
  log CPM of 500 random cohort samples, genes with mean count >= 1), whole pipeline on the pseudo genes; clock-specific if
  |pooled clock offset| exceeds the 95th percentile of |pooled pseudo offsets|.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
A164 = PAPER / "exploratory/A164_gtex_metabolic_matched/scripts"
sys.path.insert(0, str(A164))  # compute_a164 puts paper/ and gtex_cues on the path
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import compute_a164 as A  # noqa: E402  (labels, run_chunks, p2, seed_of, say, Z90, MDE_K)
import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402

OUT = HERE.parent / "output"
COH = {"vent": (0,), "fast": (1, 2)}
SEX_SPECIFIC = ("Ovary", "Uterus", "Vagina", "Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", "Prostate",
                "Testis", "Breast - Mammary Tissue")
MIN_PER_SEX, MIN_PAIRS = 20, 10
D_BOUND, OFF_BOUND, A99_OFFSET = 0.15, 1.0, 0.35
STATS = ("d_star", "offset_h")
say, Z90, MDE_K = A.say, A.Z90, A.MDE_K


# ---- matching (public annotations only) --------------------------------------------------------------------------------
def match_sex(q: pd.DataFrame) -> list:
    """1:1 sets [woman, man] within Hardy class x age bracket (donors of `q`: male, hardy, age, ischaemic_h)."""
    q = q.dropna(subset=["age"]).copy()
    q["ischaemic_h"] = q.ischaemic_h.fillna(q.ischaemic_h.median())
    sets = []
    for _, g in q.groupby(["hardy", "age"]):
        women, men = sorted(g.index[g.male == 0]), list(g.index[g.male == 1])
        if not women or not men:
            continue
        cand = {w: sorted(men, key=lambda m: (abs(g.ischaemic_h[m] - g.ischaemic_h[w]), m)) for w in women}
        used = set()
        for w in sorted(women, key=lambda w: (len(cand[w]), w)):
            m = next((x for x in cand[w] if x not in used), None)
            if m is not None:
                used.add(m)
                sets.append([w, m])
    return sets


def make_plan(q, pb, tissues=None):
    """({tissue: sets}, matching rows) of one cohort; tissues with < MIN_PER_SEX donors of a sex or < MIN_PAIRS pairs out."""
    plan, rows = {}, []
    for t in sorted(set(pb.compartment)):
        if t in SEX_SPECIFIC or (tissues and t not in tissues):
            continue
        qt = q.loc[q.index.isin(set(pb.donor[pb.compartment == t]))]
        nf, nm = int((qt.male == 0).sum()), int((qt.male == 1).sum())
        if min(nf, nm) < MIN_PER_SEX:
            continue
        sets = match_sex(qt)
        rows.append(dict(tissue=t, n_female=nf, n_male=nm, n_pairs=len(sets), eligible=len(sets) >= MIN_PAIRS))
        if len(sets) >= MIN_PAIRS:
            plan[t] = sets
    return plan, rows


# ---- estimation (adapts A164 dstar / perm_chunk / boot_chunk: two readouts) -------------------------------------------
def stat2(ds, lab, t) -> tuple:
    try:
        tab = ds.contrast(lab, targets=[t], n_perm=0, n_boot=0, blocks="none").table
        if t in tab.index and "d_star" in tab.columns:
            return float(tab.loc[t, "d_star"]), float(tab.loc[t, "offset_diff_h"])
    except Exception:
        pass
    return np.nan, np.nan


def perm_chunk(ds, plan, seeds):
    limit_threads(1)
    donors = sorted({d for s in plan.values() for x in s for d in x})
    out = []
    for sd in seeds:
        u = dict(zip(donors, np.random.default_rng(sd).random(len(donors))))
        out.append([stat2(ds, A.labels(plan[t], key=u), t) for t in plan])
    return np.array(out)  # (replicates, tissues, 2)


def boot_chunk(pb, anchors, plan, seeds):
    """Matched-set bootstrap (A164 `boot_chunk`): women resampled jointly over tissues, host refitted on the replicate."""
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
        out.append([stat2(ds, pd.concat([A.labels([s], copy=j) for s, j in sets[t]]), t) for t in plan])
    return np.array(out)


# ---- pooling ----------------------------------------------------------------------------------------------------------
def wmean(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Row-wise weighted mean over the last axis, dropping non-finite entries (weights renormalised)."""
    ok = np.isfinite(x)
    num = np.where(ok, x, 0) @ w if x.ndim == 2 else np.where(ok, x, 0) * w
    den = ok.astype(float) @ w if x.ndim == 2 else ok * w
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def tissue_weights(null_c: np.ndarray) -> np.ndarray:
    """Inverse-variance weights over tissues from the permutation null (replicates x tissues)."""
    sd = np.nanstd(null_c, axis=0, ddof=1)
    w = np.where(np.isfinite(sd) & (sd > 0), 1 / sd ** 2, 0.0)
    return w / w.sum()


def pooled_null(nulls: dict, j: int):
    """(pooled null (replicates,), tissue weights per cohort, cohort weights, cohort nulls) for readout j."""
    tw = {c: tissue_weights(a[:, :, j]) for c, a in nulls.items()}
    coh = {c: wmean(a[:, :, j], tw[c]) for c, a in nulls.items()}
    sd = {c: float(np.nanstd(v, ddof=1)) for c, v in coh.items()}
    cw = {c: 1 / s ** 2 for c, s in sd.items() if np.isfinite(s) and s > 0}
    tot = sum(cw.values())
    cw = {c: v / tot for c, v in cw.items()}
    return sum(cw[c] * coh[c] for c in cw), tw, cw, coh


def pool_estimate(est: dict, tw: dict, cw: dict) -> tuple:
    """Pooled value and per-cohort values of per-tissue estimates (cohort -> array over tissues)."""
    per = {c: float(wmean(np.asarray(est[c], float), tw[c])) for c in cw}
    return sum(cw[c] * per[c] for c in cw if np.isfinite(per[c])) / sum(cw[c] for c in cw if np.isfinite(per[c])), per


# ---- stages -----------------------------------------------------------------------------------------------------------
def setup_data(args, out):
    warnings.filterwarnings("ignore", message="below the paper settings")
    pbs, dss, plans, mrows, anchors = {}, {}, {}, [], {}
    for c, h in COH.items():
        pb, don = gtex.load(tissues=None, hardy=h, genes="clock", depth="library_size", require_tod=(c == "fast"))
        if args.subsample:
            keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
            pb, don = pb[pb.donor.isin(keep)].reset_index(drop=True), don.loc[don.index.isin(keep)]
        plan, rows = make_plan(don[["hardy", "male", "age", "ischaemic_h"]], pb, args.tissues)
        keep = set(plan) | set(S.TISSUES)
        pb = pb[pb.compartment.isin(keep)].reset_index(drop=True)
        anchors[c] = [a for a in S.TISSUES if a in set(pb.compartment)]
        pbs[c], plans[c] = pb, plan
        dss[c] = cs.Dataset(pb, anchors=anchors[c], name=f"gtex_sex_{c}")
        mrows += [dict(cohort=c, **r) for r in rows]
        say(f"{c}: {don.shape[0]} donors ({int((don.male == 0).sum())} women, {int((don.male == 1).sum())} men), "
            f"{len(plan)} eligible tissues, {sum(len(s) for s in plan.values())} tissue pairs")
    M = pd.DataFrame(mrows)
    M.to_csv(out / "matching.csv", index=False)
    return pbs, dss, plans, anchors, M


def null_draws(args, out, quick, dss, plans):
    n_perm = args.n_perm or (20 if quick else 1000)
    nulls = {}
    for c in COH:
        if plans[c]:
            nulls[c] = A.cached(out / "draws" / f"null_{c}.npy", quick, perm_chunk, n_perm, args.jobs,
                                A.seed_of("A208", c, args.seed), f"null {c}", dss[c], plans[c])
    return nulls


def power(nulls, plans, M, out):
    rows, trows = [], []
    for j, stat in enumerate(STATS):
        null, tw, cw, coh = pooled_null(nulls, j)
        for c in nulls:
            sd = float(np.nanstd(coh[c], ddof=1))
            rows.append(dict(statistic=stat, set=c, n_tissues=len(plans[c]),
                             n_pairs=int(M[(M.cohort == c) & M.eligible].n_pairs.sum()), null_sd=sd, mde_80=MDE_K * sd,
                             weight=cw.get(c, np.nan)))
            for t, w, s in zip(plans[c], tw[c], np.nanstd(nulls[c][:, :, j], axis=0, ddof=1)):
                trows.append(dict(statistic=stat, cohort=c, tissue=t, null_sd=float(s), weight=float(w)))
        sd = float(np.nanstd(null, ddof=1))
        rows.append(dict(statistic=stat, set="pooled", n_tissues=sum(len(p) for p in plans.values()),
                         n_pairs=int(M[M.eligible].n_pairs.sum()), null_sd=sd, mde_80=MDE_K * sd, weight=1.0))
    P = pd.DataFrame(rows)
    P["bound"] = P.statistic.map({"d_star": D_BOUND, "offset_h": OFF_BOUND})
    P["powered_vs_bound"] = P.mde_80 <= P.bound  # H208.1 reads this for d* (bound 0.15); the offset bound is 1 h
    P.to_csv(out / "power.csv", index=False)
    pd.DataFrame(trows).to_csv(out / "power_tissues.csv", index=False)
    say("power\n" + P.round(4).to_string(index=False))


def pseudo_sets(hardy, n_sets, seed):
    """n_sets expression-matched pseudo-clock waveforms (genes outside the clock set; match on mean log fraction)."""
    import anndata as ad
    from clocksync.paired import pseudo_waveforms
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        idx = np.flatnonzero(o.DTHHRDY.isin(list(hardy)).to_numpy() & ~o.tissue.astype(str).str.startswith("Cells -").to_numpy())
        idx = np.sort(np.random.default_rng(seed).choice(idx, 500, replace=False))
        X = a.X[idx]
        lib = o.library_size.to_numpy(float)[idx]
        frac = np.asarray((X.multiply(1 / lib[:, None]) if hasattr(X, "multiply") else X / lib[:, None]).mean(0)).ravel()
        cnt = np.asarray(X.mean(0)).ravel()
        lm = pd.Series(np.log(frac + 1e-12), index=a.var_names)[cnt >= 1]
    finally:
        a.file.close()
    return pseudo_waveforms(lm, n_sets, np.random.default_rng(seed), cs.gtex())


def score(args, out, quick, pbs, dss, plans, anchors, M, nulls):
    n_boot = args.n_boot or (10 if quick else 200)
    n_pseudo = 5 if quick else 50
    boots, obs, trows = {}, {}, []
    for c in nulls:
        boots[c] = A.cached(out / "draws" / f"boot_{c}.npy", quick, boot_chunk, n_boot, args.jobs,
                            A.seed_of("A208", c, args.seed) + 10 ** 6, f"boot {c}", pbs[c], anchors[c], plans[c])
        res = {t: dss[c].contrast(A.labels(plans[c][t]), targets=[t], n_perm=0, n_boot=0, blocks="none").table.loc[t]
               for t in plans[c]}
        obs[c] = np.array([[r.get("d_star", np.nan), r.get("offset_diff_h", np.nan)] for r in res.values()], float)
        base = {t: np.nanmean([r.get("c_control", np.nan), r.get("c_case", np.nan)]) for t, r in res.items()}
        for k, (t, r) in enumerate(res.items()):
            trows.append(dict(cohort=c, tissue=t, n_pairs=len(plans[c][t]), c_male=r.get("c_control"),
                              c_female=r.get("c_case"), d=r.get("d"), d_star=obs[c][k, 0], offset_h=obs[c][k, 1],
                              null_sd_d_star=float(np.nanstd(nulls[c][:, k, 0], ddof=1)),
                              null_sd_offset=float(np.nanstd(nulls[c][:, k, 1], ddof=1)), baseline_c=base[t]))
    T = pd.DataFrame(trows)
    T.to_csv(out / "score_tissues.csv", index=False)
    rows = []
    for j, stat in enumerate(STATS):
        null, tw, cw, _ = pooled_null(nulls, j)
        bound = D_BOUND if j == 0 else OFF_BOUND
        D, per = pool_estimate({c: obs[c][:, j] for c in cw}, tw, cw)
        bt = sum(cw[c] * wmean(boots[c][:, :, j], tw[c]) for c in cw)
        se_b, sd_n = float(np.nanstd(bt, ddof=1)), float(np.nanstd(null, ddof=1))
        hw = Z90 * max(se_b, sd_n)
        rows.append(dict(statistic=stat, set="pooled", estimate=D, boot_se=se_b, null_sd=sd_n, lo90=D - hw, hi90=D + hw,
                         bound=bound, mde_80=MDE_K * sd_n, p_perm=A.p2(D, null), equivalent=bool(D - hw > -bound and D + hw < bound)))
        for c in cw:
            nc = wmean(nulls[c][:, :, j], tw[c])
            bc = wmean(boots[c][:, :, j], tw[c])
            rows.append(dict(statistic=stat, set=c, estimate=per[c], boot_se=float(np.nanstd(bc, ddof=1)),
                             null_sd=float(np.nanstd(nc, ddof=1)), p_perm=A.p2(per[c], nc), bound=bound))
        if j == 0:  # companion 1: the 3 lowest-coupling tissues per cohort, same permutations
            low = {c: np.argsort(T[T.cohort == c].baseline_c.to_numpy())[:3] for c in cw}
            lw = {c: np.isin(np.arange(len(plans[c])), low[c]).astype(float) * tw[c] for c in cw}
            lw = {c: v / v.sum() for c, v in lw.items()}
            nl = sum(cw[c] * wmean(nulls[c][:, :, 0], lw[c]) for c in cw)
            Dl, _ = pool_estimate({c: obs[c][:, 0] for c in cw}, lw, cw)
            rows.append(dict(statistic="d_star_3_lowest_coupling", set="pooled", estimate=Dl, null_sd=float(np.nanstd(nl, ddof=1)),
                             p_perm=A.p2(Dl, nl), bound=D_BOUND))
    S_ = pd.DataFrame(rows)
    # companion 2: pseudo-clock offsets, whole pipeline on expression-matched gene sets (observed values only)
    off = S_[(S_.statistic == "offset_h") & (S_.set == "pooled")].iloc[0]
    tw, cw = pooled_null(nulls, 1)[1:3]
    pseudo = {c: pseudo_sets(COH[c], n_pseudo, A.seed_of("A208pseudo", c, args.seed)) for c in cw}
    vals = []
    for i in range(n_pseudo):
        est = {}
        for c in cw:
            w = pseudo[c][i]
            pb, _ = gtex.load(tissues=list(plans[c]) + anchors[c], hardy=COH[c], genes=list(w.genes), depth="library_size",
                              require_tod=False)
            ds = cs.Dataset(pb, anchors=anchors[c], waveform=w)
            est[c] = np.array([stat2(ds, A.labels(plans[c][t]), t)[1] for t in plans[c]], float)
        vals.append(pool_estimate(est, tw, cw)[0])
    vals = np.abs(np.array(vals, float))
    q95 = float(np.nanquantile(vals, 0.95))
    S_ = pd.concat([S_, pd.DataFrame([dict(statistic="offset_h_pseudo_abs_q95", set="pooled", estimate=q95,
                                           bound=OFF_BOUND, equivalent=bool(abs(off.estimate) > q95))])], ignore_index=True)
    S_["vs_A99_same_sign"] = np.where(S_.statistic == "offset_h", bool(np.sign(off.estimate) == np.sign(A99_OFFSET)), np.nan)
    S_.to_csv(out / "score_pooled.csv", index=False)
    pd.DataFrame(dict(abs_pseudo_offset=vals)).to_csv(out / "score_pseudo.csv", index=False)
    say("score\n" + S_.round(4).to_string(index=False))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--tissues", nargs="*", default=None, help="restrict to these tissues (check)")
    args = ap.parse_args(argv)
    quick = cli.setup(args)["quick"] or bool(args.tissues)
    out = OUT / "subsample" if quick else OUT
    (out / "draws").mkdir(parents=True, exist_ok=True)
    pbs, dss, plans, anchors, M = setup_data(args, out)
    say("matching (counts)\n" + M.to_string(index=False))
    nulls = null_draws(args, out, quick, dss, plans)
    if args.stage == "power":
        power(nulls, plans, M, out)
    else:
        score(args, out, quick, pbs, dss, plans, anchors, M, nulls)


if __name__ == "__main__":
    main()
