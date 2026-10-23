"""A132 (secondary, own alpha 0.05): brain death or ventilation? Collapse of organ phase relations in Hardy 0.

Registered in paper/docs/registered_experiments.md, section A132 (H132.1, H132.2, D132.3). A109 states and tissues
(GTEx v10, 12 tissues, donors with >= 2 tissues and time of death), fast deaths (Hardy 1-2) as reference.
  H132.1  Hardy 0 organ donors (COHORT "Organ Donor (OPO)") vs ventilated post-mortem cases (COHORT "Postmortem"):
          Fisher-Lee rho of the pair offsets vs fast (10,000 joint tissue-label permutations) and folded distance to
          0 / 12 h with p_closer (1,000 donor-label permutations vs fast), per group; contrast = folded distance of
          organ donors minus ventilated post-mortem cases at equal n (one draw of n_eq donors from the larger group),
          donor-label permutation p (1,000, one-sided, prediction < 0). Pairs: those with >= 20 shared donors in the
          groups compared (fixed mask), as the feasibility rule counts them.
  H132.2  Organ donors, ventilation duration (DTHVNTD, hours) top vs bottom tertile: difference in folded distance
          (top minus bottom), donor-label permutation p (one-sided, prediction < 0), equal n.
  D132.3  Long-term steroids (MHSTRDLT 1 vs 0), Hardy 0 and Hardy 4: folded distance and share of the cross-covariance
          along the glucocorticoid direction (A109 non-temporal output), arms >= 15 donors.
MDE (80 % power, one-sided 0.05) = (z.95 + z.80) x SD of the donor-label permutation null of the difference; no
observed difference is computed in --scout.
Controlled access (dbGaP phs000424.v8): phenotypes and time of death in memory only; aggregates in the outputs.

    python compute_brain_death.py [--scout] [--subsample N] [--seed S]
Outputs (paper/exploratory/A132_brain_death/output/): a132_groups.csv, a132_contrast.csv, a132_steroids.csv,
scout_variables.csv, scout_cohort_hardy.csv, scout_groups.csv, scout_mde.csv, scout_tertiles.csv, scout_steroids.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

A109 = Path(__file__).resolve().parents[2] / "A109_coupling_nulls" / "scripts"
sys.path.insert(0, str(A109))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

FOLDER = "A132_brain_death"
BODY = H.BODY
OPO, PM = "Organ Donor (OPO)", "Postmortem"
Z_MDE = 1.6449 + 0.8416  # one-sided 0.05, power 0.80
MIN_ARM = 15
PH_COLS = ["COHORT", "DTHHRDY", "DTHVNT", "DTHVNTD", "DTHVNTDU", "MHSTRDLT"]


def sizes(quick: bool) -> dict:
    return dict(n_fl=200 if quick else 10_000, n_perm=100 if quick else 1000)


# ---- statistics --------------------------------------------------------------------------------------------------------

def pair_counts(M, idx) -> np.ndarray:
    """Shared donors of the 66 pairs a < b among donors idx."""
    m = M[idx]
    n = m.T @ m
    return np.array([n[a, b] for a, b in H.PAIRS])


def pair_mask(M, *groups) -> np.ndarray:
    return np.all([pair_counts(M, g) >= H.MIN_N for g in groups], axis=0)


def phase(Zraw, M, idx, mask) -> np.ndarray:
    return H.pair_phase(N._H_of(Zraw, M, idx))[mask]


def fold(theta) -> float:
    return H.folded_distance(theta) if len(theta) else np.nan


def diff_null(Zraw, M, a, b, mask, n_perm, rng) -> np.ndarray:
    """Donor-label permutation null of d(a) - d(b) over the pooled donors (fixed pair mask)."""
    pool, na = np.concatenate([a, b]), len(a)
    out = []
    for _ in range(n_perm):
        p = rng.permutation(pool)
        out.append(fold(phase(Zraw, M, p[:na], mask)) - fold(phase(Zraw, M, p[na:], mask)))
    return np.array(out)


def mde_of(null) -> dict:
    sd = float(np.std(null, ddof=1))
    return dict(null_sd=sd, null_q05=float(np.quantile(null, 0.05)), mde_80pct=Z_MDE * sd)


def vs_fast(Zraw, M, fast, gd, group, rng, sz) -> dict:
    """Per-group readouts against fast deaths (A109 H109.2 form) over the pairs with >= 20 shared donors in both."""
    mask = pair_mask(M, fast, gd)
    th_f, th_g = phase(Zraw, M, fast, mask), phase(Zraw, M, gd, mask)
    Hg = N._H_of(Zraw, M, gd)
    iu = np.flatnonzero(mask)
    r = ST.fisher_lee(th_f, th_g)
    null = np.array([ST.fisher_lee(th_f, H.permute_labels(Hg, rng.permutation(len(BODY)))[mask])
                     for _ in range(sz["n_fl"])])
    d_f, d_g = fold(th_f), fold(th_g)
    pool, nf = np.concatenate([fast, gd]), len(fast)
    nl = []
    for _ in range(sz["n_perm"]):
        p = rng.permutation(pool)
        nl.append(fold(phase(Zraw, M, p[:nf], mask)) - fold(phase(Zraw, M, p[nf:], mask)))
    return dict(group=group, n=len(gd), rho=r, p=N.pval(r, null), folded_h=d_g, folded_fast_h=d_f,
                p_closer=N.pval(d_f - d_g, nl), n_pairs=int(len(iu)), n_fast=nf, n_perm_rho=len(null))


def contrast(Zraw, M, a, b, rng, sz, label, rng_draw) -> dict:
    """Difference of folded distance a minus b at equal n (one draw), donor-label permutation, one-sided < 0."""
    n = min(len(a), len(b))
    a_, b_ = (rng_draw.choice(g, n, replace=False) if len(g) > n else g for g in (a, b))
    mask = pair_mask(M, a_, b_)
    obs = fold(phase(Zraw, M, a_, mask)) - fold(phase(Zraw, M, b_, mask))
    null = diff_null(Zraw, M, a_, b_, mask, sz["n_perm"], rng)
    return dict(contrast=label, difference_h=obs, p=N.pval(obs, null, upper=False), n_each=n, n_pairs=int(mask.sum()),
                n_perm=len(null), **mde_of(null))


def gc_share(Zraw, M, idx) -> float:
    Zc, m = N._centre(Zraw, M, idx)
    return H.share_along(H.cross_cov_matrix(Zc, m), H.GC_H % 12)


# ---- scout -------------------------------------------------------------------------------------------------------------

def tertiles(x: np.ndarray):
    lo, hi = np.quantile(x, [1 / 3, 2 / 3])
    return lo, hi, np.flatnonzero(x <= lo), np.flatnonzero(x >= hi)


def scout(T, ph, Zraw, M, rng, sz, out):
    full = gtex.phenotypes(PH_COLS)
    var = [dict(variable=c, value=str(k), n=int(v)) for c in ("COHORT", "DTHVNT", "DTHVNTDU", "MHSTRDLT")
           for k, v in full[c].value_counts(dropna=False).items()]
    pd.DataFrame(var).to_csv(out / "scout_variables.csv", index=False)
    full = full.assign(DTHHRDY=full.DTHHRDY.fillna(-9).astype(int), COHORT=full.COHORT.fillna("NA"))
    rows = []
    for name, d in (("all_donors", full), ("analysis_universe", ph.assign(COHORT=ph.COHORT.fillna("NA")))):
        ct = pd.crosstab(d.COHORT, d.DTHHRDY.fillna(-9).astype(int))
        rows += [dict(universe=name, cohort=c, hardy=int(h), n=int(ct.loc[c, h])) for c in ct.index for h in ct.columns]
    pd.DataFrame(rows).to_csv(out / "scout_cohort_hardy.csv", index=False)

    h0 = T.group_donors("vent")
    od, pm = h0[ph.COHORT.to_numpy()[h0] == OPO], h0[ph.COHORT.to_numpy()[h0] == PM]
    fast = T.group_donors("fast")
    rows = []
    for g, d in (("organ_donor", od), ("ventilated_postmortem", pm), ("fast", fast)):
        pc = pair_counts(M, d)
        rows.append(dict(group=g, n_donors=len(d), pairs_ge20=int((pc >= H.MIN_N).sum()),
                         median_shared=float(np.median(pc)), dthvnt1=int((ph.DTHVNT.to_numpy()[d] == 1).sum())))
    ns = min(len(od), len(pm))
    pc_eq = pair_counts(M, rng.choice(od, ns, replace=False) if len(od) > ns else od)
    pc_eq2 = pair_counts(M, rng.choice(pm, ns, replace=False) if len(pm) > ns else pm)
    g = pd.DataFrame(rows)
    small = g.iloc[:2].sort_values("n_donors").iloc[0]
    g["smaller_group"] = small.group
    g["feasible_ge40_donors_ge40_pairs"] = bool(small.n_donors >= 40 and small.pairs_ge20 >= 40)
    g["pairs_ge20_at_equal_n"] = [int((pc_eq >= H.MIN_N).sum()), int((pc_eq2 >= H.MIN_N).sum()), np.nan]
    g.to_csv(out / "scout_groups.csv", index=False)
    print(g.to_string(index=False), flush=True)

    # MDE of the OD - VPM contrast: spread of the donor-label permutation null at equal n (observed not computed)
    rows = []
    if min(len(od), len(pm)) >= 4:
        a_ = rng.choice(od, ns, replace=False) if len(od) > ns else od
        b_ = rng.choice(pm, ns, replace=False) if len(pm) > ns else pm
        mk = pair_mask(M, a_, b_)
        rows.append(dict(contrast="H132.1", n_each=ns, n_pairs=int(mk.sum()),
                         **mde_of(diff_null(Zraw, M, a_, b_, mk, sz["n_perm"], rng))))
    # tertiles
    vd = ph.DTHVNTD.to_numpy(float)
    ok = od[np.isfinite(vd[od])]
    trow = []
    if len(ok) >= 9:
        lo, hi, il, ih = tertiles(vd[ok])
        bot, top = ok[il], ok[ih]
        mk = pair_mask(M, bot, top)
        trow.append(dict(n_with_duration=len(ok), n_missing=len(od) - len(ok), cut_low_h=lo, cut_high_h=hi,
                         median_h=float(np.median(vd[ok])), n_bottom=len(bot), n_top=len(top),
                         pairs_ge20_bottom=int((pair_counts(M, bot) >= H.MIN_N).sum()),
                         pairs_ge20_top=int((pair_counts(M, top) >= H.MIN_N).sum()), pairs_ge20_both=int(mk.sum())))
        n_ = min(len(bot), len(top))
        a_, b_ = (rng.choice(x, n_, replace=False) if len(x) > n_ else x for x in (top, bot))
        mk = pair_mask(M, a_, b_)
        rows.append(dict(contrast="H132.2", n_each=n_, n_pairs=int(mk.sum()),
                         **mde_of(diff_null(Zraw, M, a_, b_, mk, sz["n_perm"], rng))))
    pd.DataFrame(trow).to_csv(out / "scout_tertiles.csv", index=False)
    pd.DataFrame(rows).to_csv(out / "scout_mde.csv", index=False)
    print(pd.DataFrame(trow).round(2).to_string(index=False), "\n", pd.DataFrame(rows).round(3).to_string(index=False),
          flush=True)
    srow = []
    for g_, d in (("H0", T.group_donors("vent")), ("H4", T.group_donors("H4"))):
        s = ph.MHSTRDLT.to_numpy()[d]
        srow.append(dict(group=g_, n_yes=int((s == 1).sum()), n_no=int((s == 0).sum()),
                         n_other_or_missing=int(((s != 1) & (s != 0)).sum()),
                         pairs_ge20_yes=int((pair_counts(M, d[s == 1]) >= H.MIN_N).sum())))
    pd.DataFrame(srow).to_csv(out / "scout_steroids.csv", index=False)
    print(pd.DataFrame(srow).to_string(index=False), flush=True)


# ---- scoring -----------------------------------------------------------------------------------------------------------

def score(T, ph, Zraw, M, rng, sz, out, seed):
    rd = np.random.default_rng(132 + seed)
    h0, fast = T.group_donors("vent"), T.group_donors("fast")
    coh = ph.COHORT.to_numpy()
    od, pm = h0[coh[h0] == OPO], h0[coh[h0] == PM]
    vd = ph.DTHVNTD.to_numpy(float)
    ok = od[np.isfinite(vd[od])]
    _, _, il, ih = tertiles(vd[ok])
    bot, top = ok[il], ok[ih]
    ns = min(len(od), len(pm))
    eq = lambda d, n: rd.choice(d, n, replace=False) if len(d) > n else d  # noqa: E731
    gr = [vs_fast(Zraw, M, fast, eq(d, ns), g, rng, sz)
          for g, d in (("organ_donor", od), ("ventilated_postmortem", pm))]
    gr += [vs_fast(Zraw, M, fast, d, g + "_full", rng, sz) for g, d in (("organ_donor", od), ("ventilated_postmortem", pm))]
    nt = min(len(top), len(bot))
    gr += [vs_fast(Zraw, M, fast, eq(d, nt), g, rng, sz) for g, d in (("duration_top_tertile", top), ("duration_bottom_tertile", bot))]
    g = pd.DataFrame(gr)[["group", "n", "rho", "p", "folded_h", "p_closer", "n_pairs", "folded_fast_h", "n_fast", "n_perm_rho"]]
    g.to_csv(out / "a132_groups.csv", index=False)
    c = pd.DataFrame([contrast(Zraw, M, od, pm, rng, sz, "H132.1_organ_donor_minus_ventilated_postmortem", rd),
                      contrast(Zraw, M, top, bot, rng, sz, "H132.2_top_minus_bottom_tertile", rd)])
    c.to_csv(out / "a132_contrast.csv", index=False)
    rows = []
    for gname, grp in (("H0", T.group_donors("vent")), ("H4", T.group_donors("H4"))):
        s = ph.MHSTRDLT.to_numpy()[grp]
        for arm, v in (("yes", 1), ("no", 0)):
            d = grp[s == v]
            fl = len(d) >= MIN_ARM
            mk = pair_mask(M, d) if fl else None
            rows.append(dict(group=gname, steroids=arm, n=len(d), reported=fl,
                             folded_h=fold(phase(Zraw, M, d, mk)) if fl else np.nan,
                             n_pairs=int(mk.sum()) if fl else np.nan, gc_share=gc_share(Zraw, M, d) if fl else np.nan))
    pd.DataFrame(rows).to_csv(out / "a132_steroids.csv", index=False)
    for t in (g, c, pd.DataFrame(rows)):
        print(t.round(4).to_string(index=False), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--scout", action="store_true", help="group sizes, pairs, tertiles, steroids, MDE; no scoring")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    q = run["quick"]
    sz = sizes(q)
    out = exploratory_dir(FOLDER, q)
    rng = np.random.default_rng(1320 + args.seed)
    T = H.Tissues(1, args.seed, args.subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    ph = gtex.phenotypes(PH_COLS).reindex(T.donors)
    Zraw, M, _ = N._raw_matrix(T)
    scout(T, ph, Zraw, M, rng, sz, out)
    if not args.scout:
        score(T, ph, Zraw, M, rng, sz, out, args.seed)


if __name__ == "__main__":
    main()
