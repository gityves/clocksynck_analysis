"""A110 (what keeps organs coupled when daily cues are disrupted by illness?): scout and scoring.

Registered in paper/docs/registered_experiments.md, section A110. Data, tissues, states and C as A109
(GTEx v10, `A109_coupling_nulls/scripts/_helpers.py`); groups: Hardy 4 (main), Hardy 3, fast deaths (Hardy 1-2),
Hardy 0 organ donors (COHORT "Organ Donor (OPO)").
  H110.1  Hardy 4 donors split by the interval from onset of the immediate cause of death to death (DTHCODD, unit
          DTHCODDU; DTHFUCODD if < 60 % are complete). C(top tertile) - C(bottom tertile), equal n, median split if a
          tertile has < 40 pairs with >= 20 shared donors; one-sided (top < bottom) donor-label permutation p (1,000).
          MDE = 2.49 x SD of the permutation null; memory-model drop C0 (1 - exp(-(2 pi sigma t / 24)^2 / 2)).
  H110.2  (a) Spearman(own amplitude, coupling) over peripheral tissues in Hardy 4 (A104 statistic, one-sided > 0,
          tissue-label permutation); (b) retained coupling Hardy 4 / fast at the common spread d* vs the Eraslan
          stromal-vascular nucleus fraction (one-sided Spearman > 0); Holm over (a, b); (c) descriptive hub removal.
  H110.3  Temperature scores (CIRBP, RBM3; HSPA1A, HSPA1B, HSPH1, DNAJB1; residualised on ischaemic time and its
          square per tissue) regressed out of the Hardy 4 states: drop in C, donor bootstrap, 200 expression-matched
          random score sets.
  H110.4  Brain - Hypothalamus minus the median other brain region, coupling to the peripheral host at the common
          spread d*, Hardy 3-4 vs Hardy 0 organ donors (one-sided, donor-label permutation).
  H110.5  Within-digestive minus digestive-to-other coherence, Hardy 4 minus fast (group-label permutation).
--scout computes counts, completeness, pair structure and permutation-null spreads only: no observed statistic.

    python compute_mechanism.py [--scout] [--subsample N] [--seed S] [--only h1 h2 h3 h4 h5]
Outputs (paper/exploratory/A110_mechanism/output/): scout_*.csv, a110_*.csv. Controlled access (dbGaP phs000424):
phenotypes, time of death and the interval ill in memory only; aggregates in the outputs.
"""

from __future__ import annotations

import sys
from pathlib import Path

A109 = Path(__file__).resolve().parents[2] / "A109_coupling_nulls" / "scripts"
sys.path.insert(0, str(A109))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import _helpers as H  # noqa: E402
import clocksync as cs  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import DATA, exploratory_dir  # noqa: E402

FOLDER = "A110_mechanism"
BODY = H.BODY
OPO = "Organ Donor (OPO)"
Z_MDE = 1.6449 + 0.8416  # one-sided 0.05, power 0.80
C0_FILE = A109.parent / "output" / "h109_1_null.csv"  # published (scored) A109 C of Hardy 4
SIGMA_T = 0.3  # h/day, memory-model period spread (registered)
MIN_PAIRS_TERTILE = 40
HOURS_PER_UNIT = {"hour": 1.0, "day": 24.0, "week": 168.0, "month": 730.0, "year": 8766.0}
SUBSTRUCT = {"Colon - Transverse - Mixed Cell", "Colon - Transverse - Mucosa", "Colon - Transverse - Muscularis",
             "Liver - Hepatocyte", "Liver - Mixed Cell", "Liver - Portal Tract", "Pancreas - Acini", "Pancreas - Islets",
             "Pancreas - Mixed Cell", "Stomach - Mixed Cell", "Stomach - Mucosa", "Stomach - Muscularis"}
DIGESTIVE = ["Colon - Transverse", "Colon - Sigmoid", "Stomach", "Esophagus - Mucosa", "Esophagus - Muscularis", "Liver",
             "Pancreas"]
NOT_BRAIN_PERIPHERY = {"Pituitary"}
SPINAL = "Brain - Spinal cord (cervical c-1)"
HYPO = "Brain - Hypothalamus"
COLD, HEAT = ["CIRBP", "RBM3"], ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1"]
GROUP_HARDY = {"fast": (1, 2), "H3": (3,), "H4": (4,), "H0": (0,)}
ERASLAN = [DATA.parent / "CoPhaser_data" / "data" / "gtex" / "snrna" / "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad",
           DATA.parent / "data" / "gtex" / "snrna" / "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"]
MIN_DONORS_H2 = 20  # donors per tissue and group for the H110.2 fits (A33 floor; A104 used 40: see the scout)


def sizes(quick: bool) -> dict:
    return dict(n_perm=100 if quick else 1000, n_perm_tissue=1000 if quick else 10_000, n_boot=20 if quick else 200,
                n_rand=20 if quick else 200, n_hub=20 if quick else 100)


# ---- shared data -------------------------------------------------------------------------------------------------------

def obs_table(require_tod: bool = True) -> pd.DataFrame:
    """One row per (donor, tissue): Hardy class, ischaemic time (min), time-of-death present flag. Counts only (no
    expression). Cultured cells and the sub-structure labels of other studies excluded; first sample kept."""
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["subject", "tissue", "DTHHRDY", "ischemic_time", "tod"]].copy()
    finally:
        a.file.close()
    o["subject"], o["tissue"] = o.subject.astype(str), o.tissue.astype(str)
    o = o[o.DTHHRDY.isin([0, 1, 2, 3, 4]) & ~o.tissue.str.startswith("Cells -") & ~o.tissue.isin(SUBSTRUCT)]
    if require_tod:
        o = o[o.tod.notna()]
    o = o.drop_duplicates(["subject", "tissue"])
    o["has_tod"] = o.tod.notna()
    return o.drop(columns="tod").reset_index(drop=True)


def periphery(o: pd.DataFrame) -> list:
    return sorted(t for t in o.tissue.unique() if not t.startswith("Brain") and t not in NOT_BRAIN_PERIPHERY)


def brain_regions(o: pd.DataFrame) -> list:
    return sorted(t for t in o.tissue.unique() if t.startswith("Brain") and t != SPINAL)


def donor_sets(o: pd.DataFrame, hardy) -> dict:
    """tissue -> set of donors of the Hardy classes `hardy`."""
    s = o[o.DTHHRDY.isin(hardy)]
    return {t: set(g.subject) for t, g in s.groupby("tissue")}


def shared(sets: dict, a: str, b: str) -> int:
    return len(sets.get(a, set()) & sets.get(b, set()))


def eraslan_sv() -> pd.DataFrame:
    """Stromal-vascular nucleus fraction per tissue (fibroblast, myofibroblast, endothelial, pericyte; also with
    smooth-muscle myocytes), public Eraslan snRNA atlas obs (as compute_step7_confirmations.gtex)."""
    path = next(p for p in ERASLAN if p.exists())
    o = ad.read_h5ad(path, backed="r").obs[["Tissue Site Detail", "Broad cell type"]].astype(str)
    sv = o["Broad cell type"].str.match(r"Fibroblast|Myofibroblast|Endothelial|Pericyte")
    sm = o["Broad cell type"].str.contains("smooth muscle")
    tis = o["Tissue Site Detail"]
    return pd.DataFrame(dict(sv=sv.groupby(tis).mean(), sv_smooth=(sv | sm).groupby(tis).mean(),
                             nuclei=tis.value_counts()))


def interval_hours(ph: pd.DataFrame, var: str) -> np.ndarray:
    """Interval (hours) of DTHCODD / DTHFUCODD with its unit column; NaN when missing or the unit is unknown."""
    val = pd.to_numeric(ph[var], errors="coerce").to_numpy(float)
    unit = ph[var + "U"].astype(str).str.lower().str.strip().map(HOURS_PER_UNIT).to_numpy(float)
    return val * unit


def memory_drop(c0: float, t_days: float, sigma_h_per_day: float = SIGMA_T) -> float:
    """C0 - C(t), C(t) = C0 exp(-(2 pi sigma t / 24)^2 / 2)."""
    return float(c0 * (1 - np.exp(-((2 * np.pi * sigma_h_per_day * t_days / 24) ** 2) / 2)))


def sigma_for_drop(c0: float, t_days: float, drop: float) -> float:
    """Period spread (h/day) at which the memory model predicts `drop` at t_days."""
    if t_days <= 0 or drop >= c0:
        return np.nan
    return float(24 / (2 * np.pi * t_days) * np.sqrt(-2 * np.log(1 - drop / c0)))


# ---- C on donor sets (A109 helpers) ------------------------------------------------------------------------------------

def pair_counts(M, idx) -> np.ndarray:
    m = M[idx]
    n = m.T @ m
    return np.array([n[a, b] for a, b in H.PAIRS])


def pair_idx(M, *groups) -> np.ndarray:
    """Indices of the pairs with >= 20 shared donors in every group given (fixed mask)."""
    return np.flatnonzero(np.all([pair_counts(M, g) >= H.MIN_N for g in groups], axis=0))


def C_of(T, D, idx, kind="clock") -> float:
    """Mean over the pairs `idx` of the Poisson-corrected |complex correlation| of the donors D (one draw)."""
    return float(H.coherence(T, kind, [D], pair_idx=list(idx))[0])


def tertiles(x: np.ndarray):
    lo, hi = np.quantile(x, [1 / 3, 2 / 3])
    return lo, hi, np.flatnonzero(x <= lo), np.flatnonzero(x >= hi)


def split_donors(d: np.ndarray, x: np.ndarray, M) -> dict:
    """Tertile split of donors d with interval x; median split (x <= median vs > median) when either tertile has
    < 40 pairs with >= 20 shared donors."""
    lo, hi, il, ih = tertiles(x)
    bot, top = d[il], d[ih]
    pb, pt = int((pair_counts(M, bot) >= H.MIN_N).sum()), int((pair_counts(M, top) >= H.MIN_N).sum())
    out = dict(design="tertile", lo=lo, hi=hi, bot=bot, top=top, pairs_bottom=pb, pairs_top=pt)
    if min(pb, pt) < MIN_PAIRS_TERTILE:
        med = np.median(x)
        out.update(design="median", lo=med, hi=med, bot=d[x <= med], top=d[x > med])
        out.update(pairs_bottom=int((pair_counts(M, out["bot"]) >= H.MIN_N).sum()),
                   pairs_top=int((pair_counts(M, out["top"]) >= H.MIN_N).sum()))
    return out


def equal_n(a, b, rng):
    n = min(len(a), len(b))
    return tuple(rng.choice(g, n, replace=False) if len(g) > n else g for g in (a, b))


def c_diff_null(T, top, bot, idx, n_perm, rng) -> np.ndarray:
    """Donor-label permutation null of C(top) - C(bottom) (equal n, fixed pair mask `idx`)."""
    pool, n = np.concatenate([top, bot]), len(top)
    out = []
    for _ in range(n_perm):
        p = rng.permutation(pool)
        out.append(C_of(T, p[:n], idx) - C_of(T, p[n:], idx))
    return np.array(out)


def c0_hardy4() -> float:
    d = pd.read_csv(C0_FILE)
    return float(d[(d.group == "H4") & (d.null_type == "pseudo_clock")].C_clock.iloc[0])


# ---- scout -------------------------------------------------------------------------------------------------------------

def scout_h1(T, ph, M, rng, sz, out):
    h4 = T.group_donors("H4")
    # reference: null spread and MDE if every Hardy 4 donor had the interval (random split labels, pure donor-label null)
    ref = []
    for design, k in (("tertile_n_each", len(h4) // 3), ("median_split_n_each", len(h4) // 2)):
        sel = rng.permutation(h4)
        a_, b_ = sel[:k], sel[k:2 * k]
        idx = pair_idx(M, a_, b_)
        nl = c_diff_null(T, a_, b_, idx, sz["n_perm"], rng)
        ref.append(dict(design=design, n_each=k, n_pairs=len(idx), null_sd=float(np.std(nl, ddof=1)),
                        mde_80pct=Z_MDE * float(np.std(nl, ddof=1)), n_perm=len(nl)))
    pd.DataFrame(ref).to_csv(out / "scout_h1_reference_complete.csv", index=False)
    print("scout_h1_reference_complete\n", pd.DataFrame(ref).round(4).to_string(index=False), flush=True)
    c0_ = c0_hardy4()
    curve = [dict(t_diff_days=t, sigma_T_h_per_day=sg, predicted_drop=memory_drop(c0_, t, sg))
             for sg in (0.1, SIGMA_T, 1.0) for t in (1, 3, 7, 14, 30, 100, 257)]
    pd.DataFrame(curve).to_csv(out / "scout_h1_memory_curve.csv", index=False)
    n_h4 = len(h4)
    c0 = c0_hardy4()
    comp, tert, nul, mem, quant = [], [], [], [], []
    units = []
    for var in ("DTHCODD", "DTHFUCODD", "DTHLUCODD"):
        x = interval_hours(ph, var)[h4]
        raw_u = ph[var + "U"].astype(str).to_numpy()[h4]
        for k, v in pd.Series(raw_u).value_counts(dropna=False).items():
            units.append(dict(variable=var, unit=str(k), n=int(v)))
        ok = np.isfinite(x)
        comp.append(dict(variable=var, n_h4=n_h4, n_complete=int(ok.sum()), frac_complete=ok.mean(),
                         n_zero=int((x[ok] == 0).sum()), registered_primary=(var == "DTHCODD"),
                         fallback_needed=(var == "DTHCODD" and ok.mean() < 0.6)))
        d, xv = h4[ok], x[ok]
        q = np.quantile(xv, [0, 0.1, 0.25, 1 / 3, 0.5, 2 / 3, 0.75, 0.9, 1])
        quant.append(dict(variable=var, n=len(xv), **{f"q{int(round(p * 100))}_days": v / 24 for p, v in
                                                       zip([0, .1, .25, 1 / 3, .5, 2 / 3, .75, .9, 1], q)}))
        s = split_donors(d, xv, M)
        bot, top = s["bot"], s["top"]
        xb, xt = xv[np.isin(d, bot)], xv[np.isin(d, top)]
        tert.append(dict(variable=var, design=s["design"], cut_low_days=s["lo"] / 24, cut_high_days=s["hi"] / 24,
                         n_bottom=len(bot), n_top=len(top), pairs_ge20_bottom=s["pairs_bottom"],
                         pairs_ge20_top=s["pairs_top"], median_bottom_days=float(np.median(xb)) / 24,
                         median_top_days=float(np.median(xt)) / 24,
                         tertile_pairs_rule_met=bool(min(pair_counts_ge(M, tertiles(xv), d)) >= MIN_PAIRS_TERTILE)))
        a_, b_ = equal_n(top, bot, rng)
        idx = pair_idx(M, a_, b_)
        null = c_diff_null(T, a_, b_, idx, sz["n_perm"], rng) if len(idx) else np.full(2, np.nan)
        sd = float(np.std(null, ddof=1)) if len(idx) else np.nan
        t_days = (float(np.median(xt)) - float(np.median(xb))) / 24
        nul.append(dict(variable=var, design=s["design"], n_each=len(a_), n_pairs=len(idx), n_perm=len(null),
                        null_sd=sd, null_q05=float(np.quantile(null, 0.05)) if len(idx) else np.nan, mde_80pct=Z_MDE * sd))
        mem.append(dict(variable=var, design=s["design"], C0_hardy4=c0, sigma_T_h_per_day=SIGMA_T,
                        t_top_median_days=float(np.median(xt)) / 24, t_bottom_median_days=float(np.median(xb)) / 24,
                        t_diff_days=t_days, predicted_drop=memory_drop(c0, t_days),
                        predicted_drop_from_zero_to_top=memory_drop(c0, float(np.median(xt)) / 24),
                        mde_80pct=Z_MDE * sd,
                        sigma_T_needed_for_mde=sigma_for_drop(c0, t_days, Z_MDE * sd) if t_days > 0 else np.nan))
    cg = []
    hardy_ = T.hardy
    vent = pd.to_numeric(ph.reindex(T.donors).DTHVNTD, errors="coerce").to_numpy(float) if "DTHVNTD" in ph else None
    for g, hs in (("H0", (0,)), ("H3", (3,)), ("H4", (4,)), ("fast", (1, 2))):
        d = np.flatnonzero(np.isin(hardy_, hs))
        r = dict(group=g, n=len(d))
        for var in ("DTHCODD", "DTHFUCODD", "DTHLUCODD"):
            r[f"frac_{var}"] = float(np.isfinite(interval_hours(ph, var)[d]).mean())
        r["frac_any_interval"] = float(np.any([np.isfinite(interval_hours(ph, v)[d]) for v in
                                              ("DTHCODD", "DTHFUCODD", "DTHLUCODD")], axis=0).mean())
        if vent is not None:
            r["frac_DTHVNTD"] = float(np.isfinite(vent[d]).mean())
        cg.append(r)
    for name, rows in (("scout_h1_completeness", comp), ("scout_h1_completeness_by_group", cg), ("scout_h1_units", units), ("scout_h1_quantiles", quant),
                       ("scout_h1_tertiles", tert), ("scout_h1_null", nul), ("scout_h1_memory", mem)):
        df = pd.DataFrame(rows)
        df.to_csv(out / f"{name}.csv", index=False)
        print(name, "\n", df.round(4).to_string(index=False), flush=True)


def pair_counts_ge(M, tert, d):
    """Pairs with >= 20 shared donors in the bottom and top tertile (tertile definition, before any fallback)."""
    lo, hi, il, ih = tert
    return [int((pair_counts(M, d[i]) >= H.MIN_N).sum()) for i in (il, ih)]


def scout_h2(o4, o_fast, out):
    s4, sf = donor_sets(o4, (4,)), donor_sets(o_fast, (1, 2))
    tis = periphery(o4)
    e = eraslan_sv()
    rows = []
    for t in tis:
        n4, nf = len(s4.get(t, ())), len(sf.get(t, ()))
        rows.append(dict(tissue=t, n_h4=n4, n_fast=nf, ge40_both=min(n4, nf) >= 40, ge20_both=min(n4, nf) >= MIN_DONORS_H2,
                         eraslan_sv=e.sv.get(t, np.nan), eraslan_sv_smooth=e.sv_smooth.get(t, np.nan),
                         eraslan_nuclei=e.nuclei.get(t, np.nan)))
    df = pd.DataFrame(rows)
    df.to_csv(out / "scout_h2_tissues.csv", index=False)
    both = df[df.ge20_both]
    summ = dict(peripheral_tissues=len(df), tissues_ge40_h4_and_fast=int(df.ge40_both.sum()),
                tissues_ge20_h4_and_fast=int(df.ge20_both.sum()), tissues_ge15_h4_and_fast=int((df[["n_h4", "n_fast"]].min(axis=1) >= 15).sum()),
                eraslan_tissues=len(e), eraslan_matched_ge20=int(both.eraslan_sv.notna().sum()),
                eraslan_matched_ge40=int(df[df.ge40_both].eraslan_sv.notna().sum()),
                eraslan_unmatched=";".join(sorted(set(e.index) - set(df.tissue))))
    pd.DataFrame([summ]).to_csv(out / "scout_h2_summary.csv", index=False)
    print(df.to_string(index=False), "\n", pd.DataFrame([summ]).T.to_string(), flush=True)


def scout_h3(o4, out):
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        names = set(a.var_names)
        lg = a.obs[["subject", "tissue", "DTHHRDY", "ischemic_time"]].copy()
    finally:
        a.file.close()
    rows = [dict(item=f"gene:{g}", value=int(g in names)) for g in COLD + HEAT]
    lg = lg[lg.DTHHRDY == 4]
    lg = lg[~lg.tissue.astype(str).str.startswith("Cells -") & ~lg.tissue.astype(str).isin(SUBSTRUCT)]
    ok = lg.ischemic_time.notna()
    per = ok.groupby(lg.tissue.astype(str)).mean()
    rows += [dict(item="ischaemic_time_complete_frac_h4_samples", value=float(ok.mean())),
             dict(item="ischaemic_time_complete_frac_min_tissue", value=float(per.min())),
             dict(item="h4_samples", value=int(len(lg))),
             dict(item="ischaemic_time_min_h4_median", value=float(lg.ischemic_time.median())),
             dict(item="ischaemic_time_min_h4_q25", value=float(lg.ischemic_time.quantile(0.25))),
             dict(item="ischaemic_time_min_h4_q75", value=float(lg.ischemic_time.quantile(0.75)))]
    df = pd.DataFrame(rows)
    df.to_csv(out / "scout_h3.csv", index=False)
    print(df.to_string(index=False), flush=True)


def scout_h4(o_all, o_notod, ph_all, out):
    """Donors with hypothalamus and >= 3 other brain regions per group, with and without the time-of-death filter."""
    rows = []
    br = brain_regions(o_notod)
    others = [b for b in br if b != HYPO]
    cohort = ph_all.COHORT.fillna("NA")
    anyb = set(o_notod[o_notod.tissue.isin(br)].subject)
    nbrain = []
    for label, o in (("tod_required_A109", o_all), ("no_tod_filter", o_notod)):
        for g, hardy, opo in (("H34", (3, 4), None), ("H0_organ_donor", (0,), True), ("H0_all", (0,), None),
                              ("fast", (1, 2), None), ("H4", (4,), None), ("H3", (3,), None)):
            s = o[o.DTHHRDY.isin(hardy)]
            if opo:
                s = s[s.subject.map(cohort).eq(OPO)]
            hyp = set(s[s.tissue == HYPO].subject)
            n_grp, n_anyb = s.subject.nunique(), len(set(s.subject) & anyb)
            nb = s[s.tissue.isin(others)].groupby("subject").size()
            ok = {d for d in hyp if nb.get(d, 0) >= 3}
            nbody = s[s.tissue.isin(BODY)].groupby("subject").size()
            ok2 = {d for d in ok if nbody.get(d, 0) >= 2}
            rows.append(dict(universe=label, group=g, n_donors_in_group=n_grp, n_with_any_brain_region=n_anyb,
                             n_donors_with_hypothalamus=len(hyp),
                             n_hyp_and_ge3_other_brain=len(ok), n_hyp_ge3_brain_and_ge2_peripheral=len(ok2),
                             feasible_ge20=bool(len(ok2) >= 20)))
    df = pd.DataFrame(rows)
    df.to_csv(out / "scout_h4_groups.csv", index=False)
    pd.DataFrame([dict(brain_regions_other=len(others), regions=";".join(others), excluded=SPINAL)]).to_csv(
        out / "scout_h4_regions.csv", index=False)
    print(df.to_string(index=False), f"\nother brain regions: {len(others)}", flush=True)


def scout_h5(o, out):
    s4, sf = donor_sets(o, (4,)), donor_sets(o, (1, 2))
    other = [t for t in BODY if t not in DIGESTIVE]
    present = [t for t in DIGESTIVE if t in s4 or t in sf]
    rows = [dict(tissue=t, digestive=True, n_h4=len(s4.get(t, ())), n_fast=len(sf.get(t, ()))) for t in DIGESTIVE]
    rows += [dict(tissue=t, digestive=False, n_h4=len(s4.get(t, ())), n_fast=len(sf.get(t, ()))) for t in other]
    pd.DataFrame(rows).to_csv(out / "scout_h5_tissues.csv", index=False)
    pr = []
    for kind, pairs in (("within_digestive", [(a, b) for i, a in enumerate(DIGESTIVE) for b in DIGESTIVE[i + 1:]]),
                        ("digestive_to_other", [(a, b) for a in DIGESTIVE for b in other])):
        n4 = np.array([shared(s4, a, b) for a, b in pairs])
        nf = np.array([shared(sf, a, b) for a, b in pairs])
        pr.append(dict(kind=kind, pairs=len(pairs), pairs_ge20_h4=int((n4 >= 20).sum()), pairs_ge20_fast=int((nf >= 20).sum()),
                       pairs_ge20_both=int(((n4 >= 20) & (nf >= 20)).sum()), median_shared_h4=float(np.median(n4)),
                       median_shared_fast=float(np.median(nf))))
    df = pd.DataFrame(pr)
    df.to_csv(out / "scout_h5_pairs.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False), "\n", df.to_string(index=False), flush=True)
    print("digestive tissues in the A109 set:", [t for t in BODY if t in DIGESTIVE], flush=True)


def scout_h5_equal_n(sz, out, rng, sub, seed):
    """Pairs and permutation-null spread at equal n (Hardy 4 vs an equal draw of fast deaths); no observed statistic."""
    P = Panel(DIGESTIVE + [t for t in BODY if t not in DIGESTIVE], sub=sub, seed=seed)
    h4, fast, wi, cr, stat = h5_setup(P, rng)
    row = dict(n_each=len(h4), n_h4_all=int((P.hardy == 4).sum()), n_fast_all=int(np.isin(P.hardy, (1, 2)).sum()),
               pairs_within_ge20_both=len(wi), pairs_cross_ge20_both=len(cr))
    if wi and cr:
        nul = h5_null(h4, fast, stat, sz["n_perm"], rng)
        sd = float(np.std(nul, ddof=1))
        row.update(null_sd=sd, null_q05=float(np.quantile(nul, 0.05)), mde_80pct=Z_MDE * sd, n_perm=len(nul))
    pd.DataFrame([row]).to_csv(out / "scout_h5_equal_n.csv", index=False)
    print(pd.DataFrame([row]).round(4).T.to_string(), flush=True)


# ---- scoring -----------------------------------------------------------------------------------------------------------

def score_h1(T, ph, M, rng, sz, out, seed, group="H4"):
    rd = np.random.default_rng(110 + seed)
    h4 = T.group_donors(group)
    c0 = c0_hardy4()
    x = interval_hours(ph, "DTHCODD")[h4]
    var = "DTHCODD"
    if np.isfinite(x).mean() < 0.6:
        var, x = "DTHFUCODD", interval_hours(ph, "DTHFUCODD")[h4]
    ok = np.isfinite(x)
    d, xv = h4[ok], x[ok]
    s = split_donors(d, xv, M)
    top, bot = equal_n(s["top"], s["bot"], rd)
    idx = pair_idx(M, top, bot)
    if len(idx) == 0:
        row = dict(variable=var, n_complete=int(ok.sum()), n_h4=len(h4), pairs_ge20_both=0, reading="infeasible")
        pd.DataFrame([row]).to_csv(out / "a110_h1.csv", index=False)
        print("H110.1 infeasible:", row, flush=True)
        return
    obs = C_of(T, top, idx) - C_of(T, bot, idx)
    null = c_diff_null(T, top, bot, idx, sz["n_perm"], rng)
    xb, xt = xv[np.isin(d, s["bot"])], xv[np.isin(d, s["top"])]
    t_days = (np.median(xt) - np.median(xb)) / 24
    sd = float(np.std(null, ddof=1))
    row = dict(variable=var, design=s["design"], n_each=len(top), n_pairs=len(idx), diff_C_top_minus_bottom=obs,
               p_less=N.pval(obs, null, upper=False), null_sd=sd, mde_80pct=Z_MDE * sd, t_diff_days=float(t_days),
               predicted_memory_drop=memory_drop(c0, float(t_days)), n_perm=len(null),
               reading=("decay" if N.pval(obs, null, upper=False) < 0.05 else
                        "maintained" if Z_MDE * sd < memory_drop(c0, float(t_days)) else "uninformative"))
    pd.DataFrame([row]).to_csv(out / "a110_h1.csv", index=False)
    print(pd.DataFrame([row]).round(4).T.to_string(), flush=True)


def h2_data(sub: int | None, seed: int):
    pb, donors = gtex.load(hardy=(1, 2, 4), depth="library_size")
    pb = pb[~pb.compartment.isin(SUBSTRUCT) & ~pb.compartment.str.startswith("Brain") & (pb.compartment != "Pituitary")]
    if sub:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:sub]
        pb = pb[pb.donor.isin(keep)]
    grp = donors.hardy.map(lambda h: "H4" if h == 4 else "fast")
    n = pb.groupby(["compartment", pb.donor.map(grp)]).donor.nunique().unstack()
    min_d = 5 if sub else MIN_DONORS_H2
    tis = sorted(n.index[(n.min(axis=1) >= min_d)])
    return pb[pb.compartment.isin(tis)].reset_index(drop=True), grp, tis


def _rho(t, x="own_amplitude", y="coupling"):
    t = t[[x, y]].dropna()
    return spearmanr(t[x], t[y]).statistic if len(t) >= 5 else np.nan


def perm_rho(x, y, n, rng):
    return np.array([spearmanr(rng.permutation(x), y).statistic for _ in range(n)])


def score_h2(sz, out, seed, sub):
    rng = np.random.default_rng(1102 + seed)
    pb, grp, tis = h2_data(sub, seed)
    print("H110.2 tissues:", len(tis), flush=True)
    dss, tabs = {}, {}
    for g in ("H4", "fast"):
        d = pb[pb.donor.map(grp) == g]
        d = d[d.compartment.isin(d.groupby("compartment").donor.nunique().loc[lambda n: n >= 2].index)]
        ds = cs.Dataset(d, anchors=sorted(set(d.compartment)), min_donors=2, waveform="human")
        dss[g] = ds
        f = ds.fit(n_boot=0, gate=False).table
        tabs[g] = f[["coupling", "own_amplitude", "z2", "V", "psi", "boundary", "n"]]
    t4 = tabs["H4"]
    # (a)
    r_a = _rho(t4)
    null = perm_rho(t4.own_amplitude.to_numpy(), t4.coupling.to_numpy(), sz["n_perm_tissue"], rng)
    p_a = ST.perm_p(null, r_a, "greater", tie=0.0)
    # (b)
    donors = pd.Series(grp.reindex(pb.donor.unique()).to_numpy(), index=pb.donor.unique())
    both = cs.Dataset(pb, anchors=sorted(set(tabs["H4"].index) & set(tabs["fast"].index)), min_donors=2, waveform="human")
    ct = both.contrast(donors, levels=("fast", "H4"), targets=list(both.anchors), n_perm=0, n_boot=0).table
    V0 = ct.V0
    cs_ = {g: ct[f"z2_{g}"] * V0 / (ct[f"z2_{g}"] * V0 + 2 * ct[f"psi_{g}"]) for g in ("fast", "H4")}
    retained = (cs_["H4"] / cs_["fast"]).rename("retained")
    e = eraslan_sv()
    rows_b = []
    for col in ("sv", "sv_smooth"):
        j = pd.concat([retained, e[col]], axis=1).dropna()
        j.columns = ["retained", "frac"]
        if len(j) >= 5:
            r = spearmanr(j.frac, j.retained).statistic
            nl = perm_rho(j.frac.to_numpy(), j.retained.to_numpy(), sz["n_perm_tissue"], rng)
            rows_b.append(dict(fraction=col, tissues=len(j), rho=r, p_greater=ST.perm_p(nl, r, "greater", tie=0.0)))
        else:
            rows_b.append(dict(fraction=col, tissues=len(j), rho=np.nan, p_greater=np.nan))
    b = pd.DataFrame(rows_b)
    prim = pd.DataFrame([dict(test="H110.2a", n=int(t4[["own_amplitude", "coupling"]].notna().all(axis=1).sum()), rho=r_a,
                              p=p_a, rho_fast=_rho(tabs["fast"])),
                         dict(test="H110.2b", n=int(b.tissues.iloc[0]), rho=b.rho.iloc[0], p=b.p_greater.iloc[0],
                              rho_fast=np.nan)])
    prim["p_holm"] = ST.holm(prim.p.to_numpy())
    prim.to_csv(out / "a110_h2.csv", index=False)
    b.to_csv(out / "a110_h2b_fractions.csv", index=False)
    pd.DataFrame({"retained": retained, "c_h4_dstar": cs_["H4"], "c_fast_dstar": cs_["fast"]}).join(e).to_csv(
        out / "a110_h2b_tissues.csv")
    print(prim.round(4).to_string(index=False), "\n", b.round(4).to_string(index=False), flush=True)
    # (c) hubs, descriptive
    tis = sorted(set(tabs["H4"].index) & set(tabs["fast"].index))
    rows = []
    for g in ("H4", "fast"):
        dn = np.array(sorted(donors.index[donors == g]))
        dn = rng.permutation(dn)
        h1_, h2_ = dn[: len(dn) // 2], dn[len(dn) // 2:]
        f1 = dss[g].subset(donors=h1_).fit(n_boot=0, gate=False).table
        hubs = list(f1.own_amplitude.nlargest(4).index)
        rows.append(dict(group=g, h1=h1_, h2=h2_, hubs=hubs))

    half2 = {r["group"]: dss[r["group"]].subset(donors=r["h2"]) for r in rows}
    c_all = {g: ds2.fit(n_boot=0, gate=False).table.coupling for g, ds2 in half2.items()}  # host: all other tissues

    def drop(g, removed, hubs):
        """Mean over the targets (not removed, not the Hardy 4 hubs) of c(all-other host) - c(host without `removed`)."""
        have = set(half2[g].compartments)
        targets = [t for t in tis if t in have and t not in removed and t not in hubs]
        keep = [t for t in tis if t not in removed]
        c_wo = half2[g].subset(compartments=keep).fit(n_boot=0, gate=False, targets=targets).table.coupling
        return float(np.nanmean((c_all[g].reindex(targets) - c_wo.reindex(targets))))

    # the four hubs of Hardy 4 half 1 removed in both groups (drop defined on the tissues that are not hubs)
    hub4 = rows[0]["hubs"]
    obs = {r["group"]: drop(r["group"], hub4, hub4) for r in rows}
    rd = [dict(group=g, hubs=";".join(hub4), mean_drop=v) for g, v in obs.items()]
    nul = []
    for _ in range(sz["n_hub"]):
        rh = list(rng.choice([t for t in tis if t not in hub4], 4, replace=False))
        nul.append(drop("H4", rh, hub4) - drop("fast", rh, hub4))
    rd.append(dict(group="H4_minus_fast", hubs=";".join(hub4), mean_drop=obs["H4"] - obs["fast"],
                   null_median=float(np.median(nul)), p_greater=ST.perm_p(nul, obs["H4"] - obs["fast"], "greater"),
                   n_null=len(nul)))
    pd.DataFrame(rd).to_csv(out / "a110_h2c_hubs.csv", index=False)
    print(pd.DataFrame(rd).round(4).to_string(index=False), flush=True)


def row_expression(T, donors_sel, genes_all: bool = True):
    """Per tissue: log CPM (rows of T.rows[t] belonging to donors_sel) of all genes (float32), gene names, ischaemic
    time (min) of the rows. Row order as `H.Tissues` (sorted obs position); asserted against `T.row_donor`."""
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        keep = (o.DTHHRDY.isin([0, 1, 2, 3, 4]) & o.tissue.astype(str).isin(BODY) & o.tod.notna()).to_numpy()
        sub = o.iloc[np.flatnonzero(keep)]
        subj = sub.subject.astype(str).to_numpy()
        ok = ~pd.Series(list(zip(subj, sub.tissue.astype(str)))).duplicated().to_numpy()
        sub, subj = sub[ok], subj[ok]
        names = np.asarray(a.var_names)
        res = {}
        selset = set(T.donors[donors_sel])
        for t in BODY:
            m = (sub.tissue.astype(str) == t).to_numpy() & np.isin(subj, T.donors)
            ids = sub.index[m]
            pos = a.obs_names.get_indexer(ids)
            o_ = np.argsort(pos)
            pos, sd = pos[o_], subj[m][o_]
            assert (sd == T.row_donor[t]).all(), "row order differs from H.Tissues"
            rows = np.flatnonzero(np.isin(sd, list(selset)))
            X = a.X[pos[rows]]
            X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            lib = o.library_size.to_numpy(float)[pos[rows]]
            lg = np.log(X / lib[:, None] * 1e6 + 0.5).astype(np.float32)
            isch = o.ischemic_time.to_numpy(float)[pos[rows]]
            res[t] = (rows, lg, isch)
    finally:
        a.file.close()
    return res, names


def residual_scores(lg, isch, cols_sets):
    """Scores (mean z of log CPM within the tissue over the rows given, per gene set) residualised on ischaemic time and
    its square."""
    z = (lg - lg.mean(0)) / np.maximum(lg.std(0), 1e-6)
    s = np.column_stack([z[:, c].mean(1) for c in cols_sets])
    t = np.where(np.isfinite(isch), isch, np.nanmean(isch))
    t = (t - t.mean()) / max(t.std(), 1e-9)
    Xc = np.column_stack([np.ones(len(t)), t, t ** 2])
    return s - Xc @ np.linalg.lstsq(Xc, s, rcond=None)[0]


def regress_out(Zt, S):
    """Clock-plane states (complex) with the linear part on the scores S removed (means kept)."""
    S = S - S.mean(0)
    Y = np.column_stack([Zt.real, Zt.imag])
    B = np.linalg.lstsq(np.column_stack([np.ones(len(S)), S]), Y, rcond=None)[0][1:]
    R = Y - S @ B
    return R[:, 0] + 1j * R[:, 1]


def score_h3(T, M, rng, sz, out, seed):
    h4 = T.group_donors("H4")
    res, names = row_expression(T, h4)
    col = {g: i for i, g in enumerate(names)}
    clock = set(H.WF.genes) | H.CLOCK_OUTPUT
    mean_all = np.mean(np.vstack([lg.mean(0, keepdims=True) for _, lg, _ in res.values()]), axis=0)
    real = [[col[g] for g in COLD], [col[g] for g in HEAT]]
    pool = np.array([i for i, g in enumerate(names) if g not in clock and g not in COLD + HEAT])
    T.Z["resid"], T.V["resid"] = {}, T.V["clock"]
    Zc = T.Z["clock"]
    idx = list(range(len(H.PAIRS)))

    def set_resid(cols_sets):
        for t in BODY:
            rows, lg, isch = res[t]
            Zr = Zc[t][0].copy()
            Zr[rows] = regress_out(Zc[t][0][rows], residual_scores(lg, isch, cols_sets))
            T.Z["resid"][t] = Zr[None, :]

    c_raw = C_of(T, h4, idx)
    set_resid(real)
    c_res = C_of(T, h4, idx)
    obs = c_raw - c_res
    boot = []
    for _ in range(sz["n_boot"]):
        D = rng.choice(h4, len(h4))
        boot.append(C_of(T, D, idx) - C_of(T, D, idx, "resid"))
    drops = []
    for _ in range(sz["n_rand"]):
        sel = []
        for g in COLD + HEAT:
            near = pool[np.abs(mean_all[pool] - mean_all[col[g]]) < 0.25]
            sel.append(int(rng.choice(near if len(near) else pool)))
        set_resid([sel[:2], sel[2:]])
        drops.append(c_raw - C_of(T, h4, idx, "resid"))
    drops = np.array(drops)
    q95 = float(np.quantile(drops, 0.95))
    row = dict(C_raw=c_raw, C_resid=c_res, drop=obs, boot_lo=float(np.quantile(boot, 0.025)),
               boot_hi=float(np.quantile(boot, 0.975)), n_boot=len(boot), random_median=float(np.median(drops)),
               random_q95=q95, p_random=float((1 + np.sum(drops >= obs)) / (1 + len(drops))), n_random=len(drops),
               excluded=bool(obs < 0.05 and np.quantile(boot, 0.975) < 0.05 and obs <= q95))
    pd.DataFrame([row]).to_csv(out / "a110_h3.csv", index=False)
    print(pd.DataFrame([row]).round(4).T.to_string(), flush=True)


def score_h4(sz, out, seed, sub):
    rng_seed = 1104 + seed
    ph = gtex.phenotypes(["COHORT"])
    o = obs_table()
    br = brain_regions(o)
    others = [b for b in br if b != HYPO]
    pb, donors = gtex.load(tissues=BODY + br, hardy=(0, 1, 2, 3, 4), depth="library_size")
    coh = ph.COHORT.reindex(donors.index).fillna("NA")
    nb = pb[pb.compartment.isin(others)].groupby("donor").size()
    hyp = set(pb[pb.compartment == HYPO].donor)
    ok = {d for d in hyp if nb.get(d, 0) >= 3}
    lab = {d: ("h34" if donors.hardy[d] in (3, 4) else "opo" if donors.hardy[d] == 0 and coh[d] == OPO else
               "fast" if donors.hardy[d] in (1, 2) else None) for d in ok}
    rows = []
    for ref, label in (("opo", "H110.4_h34_vs_organ_donor"), ("fast", "H110.4_fallback_h34_vs_fast_descriptive")):
        use = {d: g for d, g in lab.items() if g in ("h34", ref)}
        if sub:
            keep = np.random.default_rng(seed).permutation(sorted(use))[:sub]
            use = {d: use[d] for d in keep}
        n = pd.Series(use).value_counts()
        if len(n) < 2 or n.min() < (3 if sub else 20):
            rows.append(dict(contrast=label, n_h34=int(n.get("h34", 0)), n_ref=int(n.get(ref, 0)), status="infeasible"))
            continue
        d = pb[pb.donor.isin(use) & pb.compartment.isin(BODY + [HYPO] + others)].reset_index(drop=True)
        ds = cs.Dataset(d, anchors=sorted(set(BODY) & set(d.compartment)), min_donors=2)
        c = ds.contrast(pd.Series(use), levels=(ref, "h34"), targets=[HYPO] + others, n_perm=sz["n_perm"], n_boot=0,
                        seed=rng_seed)
        for sc in ("d_star", "d"):
            obs = c.observed[sc]
            nl = c.null[sc]
            stat = obs[HYPO] - np.nanmedian(obs[others])
            ns = nl[HYPO] - nl[others].median(axis=1)
            rows.append(dict(contrast=label, scale=sc, n_h34=int(n["h34"]), n_ref=int(n[ref]), statistic=float(stat),
                             p_greater=ST.perm_p(ns.to_numpy(), float(stat), "greater"), n_perm=int(ns.notna().sum()),
                             null_sd=float(ns.std(ddof=1)), mde_80pct=Z_MDE * float(ns.std(ddof=1)),
                             hyp_d=float(obs[HYPO]), median_other_d=float(np.nanmedian(obs[others])), status="ok"))
        c.table.to_csv(out / f"a110_h4_regions_{ref}.csv")
    df = pd.DataFrame(rows)
    df.to_csv(out / "a110_h4.csv", index=False)
    print(df.round(4).to_string(index=False), flush=True)


class Panel:
    """Clock-plane states of an arbitrary tissue list (digestive + A109 tissues) for the donors of the given Hardy classes
    (local: `H.Tissues` is fixed to the 12 A109 tissues)."""

    def __init__(self, tissues, hardy=(1, 2, 4), sub=None, seed=0):
        pb, donors = gtex.load(tissues=tissues, hardy=hardy, depth="library_size")
        n = pb.groupby("donor").size()
        keep = sorted(n.index[n >= 2])
        if sub:
            keep = sorted(np.random.default_rng(seed).permutation(keep)[:sub])
        pb = pb[pb.donor.isin(keep)]
        self.tissues = [t for t in tissues if t in set(pb.compartment)]
        self.donors = np.array(keep)
        gi = {d: i for i, d in enumerate(self.donors)}
        self.hardy = donors.hardy.reindex(self.donors).to_numpy(int)
        self.R, self.Z, self.V = {}, {}, {}
        for t in self.tissues:
            x = pb[pb.compartment == t].sort_values("donor").reset_index(drop=True)
            z, v = H._state_set(x[list(H.WF.genes)].assign(depth=x.depth), [H.WF], x.depth.to_numpy(float))
            self.Z[t], self.V[t] = z[0], v[0]
            r = np.full(len(self.donors), -1)
            r[[gi[d] for d in x.donor]] = np.arange(len(x))
            self.R[t] = r

    def pair_C(self, a, b, D):
        """Poisson-corrected |correlation| of tissues a, b over donors D (global indices); NaN if < 20 shared."""
        ra, rb = self.R[a][D], self.R[b][D]
        m = (ra >= 0) & (rb >= 0)
        n = int(m.sum())
        if n < H.MIN_N:
            return np.nan
        ia, ib = ra[m], rb[m]
        A, B = self.Z[a][ia], self.Z[b][ib]
        A, B = A - A.mean(), B - B.mean()
        k = n - 1
        ea = (abs(A) ** 2).sum() / k - self.V[a][ia].mean()
        eb = (abs(B) ** 2).sum() / k - self.V[b][ib].mean()
        return float(abs((A * B.conj()).sum() / k / np.sqrt(max(ea, 1e-12) * max(eb, 1e-12))))

    def shared(self, a, b, D):
        return int(((self.R[a][D] >= 0) & (self.R[b][D] >= 0)).sum())


def h5_setup(P: Panel, rng):
    """Equal-n groups (Hardy 4, a random draw of as many fast donors), the within-digestive and digestive-to-other pairs
    with >= 20 shared donors in both groups (fixed mask) and the statistic of one donor group."""
    other = [t for t in BODY if t not in DIGESTIVE]
    dig = [t for t in DIGESTIVE if t in P.tissues]
    oth = [t for t in other if t in P.tissues]
    wi = [(a, b) for i, a in enumerate(dig) for b in dig[i + 1:]]
    cr = [(a, b) for a in dig for b in oth]
    h4, fast = np.flatnonzero(P.hardy == 4), np.flatnonzero(np.isin(P.hardy, (1, 2)))
    fast = rng.choice(fast, min(len(h4), len(fast)), replace=False)
    ok = lambda pairs: [p for p in pairs if min(P.shared(*p, h4), P.shared(*p, fast)) >= H.MIN_N]  # noqa: E731
    wi, cr = ok(wi), ok(cr)

    def stat(D):
        return (np.nanmean([P.pair_C(a, b, D) for a, b in wi]) - np.nanmean([P.pair_C(a, b, D) for a, b in cr]))

    return h4, fast, wi, cr, stat


def h5_null(h4, fast, stat, n_perm, rng):
    pool, n = np.concatenate([h4, fast]), len(h4)
    nul = []
    for _ in range(n_perm):
        p = rng.permutation(pool)
        nul.append(stat(p[:n]) - stat(p[n:]))
    return np.array(nul)


def score_h5(sz, out, seed, sub):
    rng = np.random.default_rng(1105 + seed)
    P = Panel(DIGESTIVE + [t for t in BODY if t not in DIGESTIVE], sub=sub, seed=seed)
    h4, fast, wi, cr, stat = h5_setup(P, rng)
    if not wi or not cr:
        print("H110.5 infeasible: no within-digestive or digestive-to-other pairs with >= 20 shared donors", flush=True)
        pd.DataFrame([dict(status="infeasible", pairs_within=len(wi), pairs_cross=len(cr))]).to_csv(out / "a110_h5.csv", index=False)
        return
    obs = stat(h4) - stat(fast)
    nul = h5_null(h4, fast, stat, sz["n_perm"], rng)
    mean = lambda pairs, D: float(np.nanmean([P.pair_C(a, b, D) for a, b in pairs]))  # noqa: E731
    p = ST.perm_p(nul, obs, "greater")
    row = dict(n_each=len(h4), pairs_within=len(wi), pairs_cross=len(cr), diff_h4_minus_fast=float(obs), p_greater=p,
               null_sd=float(np.std(nul, ddof=1)), mde_80pct=Z_MDE * float(np.std(nul, ddof=1)), n_perm=len(nul),
               within_h4=mean(wi, h4), cross_h4=mean(cr, h4), within_fast=mean(wi, fast), cross_fast=mean(cr, fast),
               excluded=bool(not (p < 0.05 and obs > 0)))
    pd.DataFrame([row]).to_csv(out / "a110_h5.csv", index=False)
    print(pd.DataFrame([row]).round(4).T.to_string(), flush=True)


# ---- main ----------------------------------------------------------------------------------------------------------------

def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--scout", action="store_true", help="counts, completeness, pair structure, null spreads; no observed statistic")
    ap.add_argument("--h1-group", default="H4", choices=["H4", "H3", "fast", "vent"],
                    help="group split by the interval (smoke tests only; the registration is H4)")
    ap.add_argument("--only", nargs="*", default=None, choices=["h1", "h2", "h3", "h4", "h5"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    q = run["quick"]
    sz = sizes(q)
    out = exploratory_dir(FOLDER, q)
    rng = np.random.default_rng(1100 + args.seed)
    todo = set(args.only or ["h1", "h2", "h3", "h4", "h5"])
    if args.scout:
        o = obs_table(True)
        o_all = obs_table(False)
        ph_all = gtex.phenotypes(["COHORT"])
        scout_h2(o, o, out)
        scout_h3(o, out)
        scout_h4(o, o_all, ph_all, out)
        scout_h5(o, out)
        scout_h5_equal_n(sz, out, rng, args.subsample, args.seed)
    if args.scout or {"h1", "h3"} & todo:
        T = H.Tissues(1, args.seed, args.subsample)
        T.build_states(np.random.default_rng(1), 0, real=True)
        ph = gtex.phenotypes(["DTHCODD", "DTHCODDU", "DTHFUCODD", "DTHFUCODDU", "DTHLUCODD", "DTHLUCODDU", "DTHVNTD"]).reindex(T.donors)
        Zraw, M, _ = N._raw_matrix(T)
        if args.scout:
            scout_h1(T, ph, M, rng, sz, out)
            return
        if "h1" in todo:
            score_h1(T, ph, M, rng, sz, out, args.seed, args.h1_group)
        if "h3" in todo:
            score_h3(T, M, rng, sz, out, args.seed)
    if "h2" in todo:
        score_h2(sz, out, args.seed, args.subsample)
    if "h4" in todo:
        score_h4(sz, out, args.seed, args.subsample)
    if "h5" in todo:
        score_h5(sz, out, args.seed, args.subsample)


if __name__ == "__main__":
    main()
