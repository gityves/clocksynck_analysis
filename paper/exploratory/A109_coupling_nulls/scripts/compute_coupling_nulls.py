"""A109: is the coupling of organs without cues a clock? (GTEx v10, 12 A92 tissues)

Pre-registered in paper/docs/registered_experiments.md, section A109 (the data-dependent completion, group sizes, equal n
and MDE, is added after the scout, before scoring).

Hypotheses
  H109.1 (confirmatory, Holm over the family of five): in ventilated donors (Hardy 0) the inter-organ coherence
         C_clock exceeds the pseudo-clock null. Estimand: C_clock minus the median C of 200 expression-matched 14-gene
         sets that carry the GTEx waveform loadings (`clocksync.paired.pseudo_waveforms`, matched per tissue on mean log
         CPM, clock, clock-like and D6 clock-output genes excluded); one-sided empirical p (add-one).
  H109.2 (confirmatory): offset preservation. Fisher-Lee circular correlation, over the 66 tissue pairs, of the phase of
         the complex cross-covariance arg E[z_a conj z_b] in ventilated donors vs fast deaths; null: tissue labels
         permuted jointly in the ventilated group (10,000), one-sided. Second condition: the median absolute pair
         phase folded to [0, 12 h] is not closer to 0 / 12 h than in fast deaths (donor labels permuted between the
         groups, one-sided p for "closer").
  Secondary: (a) H109.1 in Hardy 3, Hardy 4, fast deaths (positive reference) and against the permuted-acrophase null;
  (b) excess coherence per group with donor bootstrap and contrasts; (c) non-temporal shared direction (fast
  deaths: states regressed on cos, sin of time of death; residual cross-tissue covariance) and the share of each
  group's cross-covariance along it; (d) pairs collected with a different kit (whole blood, SMCAT) vs the other
  pairs; (e) H109.1 and H109.2 in Hardy 4; (f) circularity ratio |E[z_a z_b]| / |E[z_a conj z_b]| against simulations of
  rotating clocks, displacement along the glucocorticoid direction (17.7 h) and 25/50/75 % mixtures at the observed
  design; descriptive host statistics.
Null types: pseudo_clock (primary), acrophase_perm (secondary), shuffled_pairing (sanity, expected C about 0).
Clock states: GTEx waveform projection per tissue, Poisson variance, centred within tissue and group (A92/A93);
equal n per group by subsampling as A93 (n = number of Hardy 4 donors; median C over the draws; same draws for the data
and every null set). GTEx time of death is controlled access: aggregates only.

    python compute_coupling_nulls.py [--subsample N] [--seed S] [--jobs 1] [--scout]
Outputs (paper/exploratory/A109_coupling_nulls/output/, .../subsample/ for checks):
  h109_1_null.csv, clock_plane_variance.csv, h109_2_offsets.csv, pair_relations.csv, excess_bootstrap.csv,
  nontemporal_share.csv, circularity_sims.csv, procedure_check.csv, h109_2_power.csv (and scout_*.csv with --scout:
  group sizes, pair counts, MDE; the scout never computes C on the real clock genes of ventilated donors).
Base code: figure_3/scripts/gtex_cues/compute_coupling_sources.py (A92/A93), compute_glucocorticoid.py (A27), A26.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

FOLDER = "A109_coupling_nulls"
BODY, GROUPS = H.BODY, H.GROUPS


def sizes(quick: bool) -> dict:
    return dict(n_sets=20 if quick else 200, n_sub=5 if quick else 20, n_shuf=20 if quick else 200,
                n_boot=20 if quick else 200, n_fl=200 if quick else 10_000, n_closer=100 if quick else 1000,
                n_sim=20 if quick else 200, n_power=30 if quick else 200)


def pval(obs, null, upper=True):
    null = np.asarray(null)
    null = null[np.isfinite(null)]
    return float((1 + np.sum(null >= obs if upper else null <= obs)) / (1 + len(null)))


# ---- H109.1 ----------------------------------------------------------------------------------------------------------

def h109_1(T, draws, n_eq, rng, sz, pair_idx=None, label=""):
    rows = []
    for g in GROUPS:
        D = draws[g]
        c = H.coherence(T, "clock", D, pair_idx=pair_idx)[0]
        npair = len(H.PAIRS) if pair_idx is None else len(pair_idx)
        for kind, name in (("pseudo", "pseudo_clock"), ("acro", "acrophase_perm")):
            nul = H.coherence(T, kind, D, pair_idx=pair_idx)
            rows.append(dict(group=g, null_type=name, C_clock=c, null_median=np.median(nul), excess=c - np.median(nul),
                             p_one_sided=pval(c, nul), null_sd=nul.std(ddof=1), null_q95=np.quantile(nul, 0.95),
                             n_donors=len(T.group_donors(g)), n_equal=min(n_eq, len(T.group_donors(g))),
                             n_sets=len(nul), n_pairs=npair))
        if pair_idx is None:
            gd = T.group_donors(g)
            sh = [H.coherence(T, "clock", D, R=H.shuffled_R(T, gd, rng))[0] for _ in range(sz["n_shuf"])]
            rows.append(dict(group=g, null_type="shuffled_pairing", C_clock=c, null_median=np.median(sh),
                             excess=c - np.median(sh), p_one_sided=pval(c, sh), null_sd=np.std(sh, ddof=1),
                             null_q95=np.quantile(sh, 0.95), n_donors=len(gd), n_equal=min(n_eq, len(gd)),
                             n_sets=len(sh), n_pairs=npair))
    return pd.DataFrame(rows).assign(subset=label)


# ---- clock-plane variance (A104: z2 V + 2 psi) -------------------------------------------------------------------------

def clock_plane_variance(T, draws):
    rows = []
    for g in GROUPS:
        gd = T.group_donors(g)
        for tag, D in (("full", gd), ("equal_n_draw", draws[g][0])):
            parts = []
            for t in BODY:
                r = T.R[t][D]
                r = r[r >= 0]
                x = T.rows[t].iloc[r].copy()
                x.insert(0, "donor", T.row_donor[t][r])
                x.insert(1, "compartment", t)
                parts.append(x)
            pb = pd.concat(parts, ignore_index=True)
            tab = cs.Dataset(pb, anchors=BODY, name=g).fit(n_boot=0, gate=False).table
            tab = tab[["coupling", "own_amplitude", "z2", "V", "psi", "boundary", "n"]].copy()
            tab["clock_plane_variance"] = tab.z2 * tab.V + 2 * tab.psi
            rows.append(tab.reset_index(names="tissue").assign(group=g, donors=tag))
    return pd.concat(rows, ignore_index=True)


# ---- H109.2 ----------------------------------------------------------------------------------------------------------

def _H_of(Zraw, M, idx):
    m = M[idx]
    mu = (Zraw[idx] * m).sum(0) / np.maximum(m.sum(0), 1)
    Zc = (Zraw[idx] - mu) * m
    return H.cross_moments(Zc, m)[0]


def h109_2(T, rng, sz, n_eq):
    """Full groups (primary) and an equal-n sensitivity row (one draw of n_eq donors per group, fast included)."""
    Zraw, M, _ = _raw_matrix(T)
    fast_all = T.group_donors("fast")
    rows = []
    for design in ("full", "equal_n"):
        fast = fast_all if design == "full" else rng.choice(fast_all, min(n_eq, len(fast_all)), replace=False)
        th_f = H.pair_phase(_H_of(Zraw, M, fast))
        for g in ("vent", "H4", "H3"):
            gd = T.group_donors(g)
            if design == "equal_n":
                gd = rng.choice(gd, min(n_eq, len(gd)), replace=False)
            Hg = _H_of(Zraw, M, gd)
            th_g = H.pair_phase(Hg)
            r = ST.fisher_lee(th_f, th_g)
            null = np.array([ST.fisher_lee(th_f, H.permute_labels(Hg, rng.permutation(len(BODY))))
                             for _ in range(sz["n_fl"])])
            d_f, d_g = H.folded_distance(th_f), H.folded_distance(th_g)
            pool = np.concatenate([fast, gd])
            nf = len(fast)
            obs = d_f - d_g  # > 0: the group is closer to 0 / 12 h than fast deaths
            nl = []
            for _ in range(sz["n_closer"]):
                p = rng.permutation(pool)
                nl.append(H.folded_distance(H.pair_phase(_H_of(Zraw, M, p[:nf]))) -
                          H.folded_distance(H.pair_phase(_H_of(Zraw, M, p[nf:]))))
            rows.append(dict(group=g, design=design, rho_fisher_lee=r, p=pval(r, null), null_median=np.median(null),
                             median_folded_phase_h=H.folded_median(th_g), median_folded_phase_fast_h=H.folded_median(th_f),
                             folded_distance_h=d_g, folded_distance_fast_h=d_f, p_closer=pval(obs, nl),
                             n_pairs=len(th_g), n_donors=len(gd), n_fast=nf))
    return pd.DataFrame(rows)


def _raw_matrix(T, D=None):
    """Donor x tissue states (not centred), mask, variance for all donors (or D)."""
    D = np.arange(len(T.donors)) if D is None else D
    Z = np.zeros((len(D), len(BODY)), complex)
    M = np.zeros((len(D), len(BODY)))
    V = np.zeros((len(D), len(BODY)))
    for j, t in enumerate(BODY):
        r = T.R[t][D]
        m = r >= 0
        Z[m, j] = T.Z["clock"][t][0][r[m]]
        V[m, j] = T.V["clock"][t][0][r[m]]
        M[m, j] = 1
    return Z, M, V


def h109_2_power(T, rng, sz):
    """MDE of the Fisher-Lee test at 66 pairs: the fast-death pair phases (positive reference, no ventilated data) plus
    wrapped-normal noise of sd sigma per pair stand for the ventilated phases; labels permuted jointly (500)."""
    Zraw, M, _ = _raw_matrix(T)
    Hf = _H_of(Zraw, M, T.group_donors("fast"))
    th_f = H.pair_phase(Hf)
    iu = np.triu_indices(len(BODY), 1)
    rows = []
    for sigma in (0.3, 0.6, 0.9, 1.2, 1.6, 2.0, 2.5):
        rs, ps = [], []
        for _ in range(sz["n_power"]):
            e = np.zeros((len(BODY),) * 2)
            e[iu] = rng.normal(0, sigma, len(iu[0]))
            e = e - e.T
            Hv = Hf * np.exp(1j * e)
            r = ST.fisher_lee(th_f, H.pair_phase(Hv))
            null = [ST.fisher_lee(th_f, H.permute_labels(Hv, rng.permutation(len(BODY)))) for _ in range(500)]
            rs.append(r)
            ps.append(pval(r, null))
        ps = np.array(ps)
        rows.append(dict(sigma_rad=sigma, mean_rho=float(np.mean(rs)), power_p05=float(np.mean(ps < 0.05)),
                         power_p01=float(np.mean(ps < 0.01)), n_sim=sz["n_power"]))
    return pd.DataFrame(rows)


# ---- pair relations, procedure, bootstrap ------------------------------------------------------------------------------

def pair_relations(T, kits):
    Zraw, M, _ = _raw_matrix(T)
    rows = []
    for g in GROUPS:
        gd = T.group_donors(g)
        Hm, Pc, N = H.cross_moments(*_centre(Zraw, M, gd))
        cc = H.pair_corr(T, "clock", gd)[0]
        for j, (a, b) in enumerate(H.PAIRS):
            rows.append(dict(group=g, tissue_a=BODY[a], tissue_b=BODY[b], n_shared=int(N[a, b]),
                             offset_h=float(np.angle(Hm[a, b]) * H.HOURS), cov_modulus=float(abs(Hm[a, b])),
                             corr=float(abs(cc[j])), pseudo_cov_ratio=float(abs(Pc[a, b]) / abs(Hm[a, b])),
                             kit_a=kits.get(BODY[a], ""), kit_b=kits.get(BODY[b], ""),
                             kit_differs=kits.get(BODY[a]) != kits.get(BODY[b])))
    return pd.DataFrame(rows)


def _centre(Zraw, M, idx):
    m = M[idx]
    mu = (Zraw[idx] * m).sum(0) / np.maximum(m.sum(0), 1)
    return (Zraw[idx] - mu) * m, m


def procedure_check(T, draws, n_eq, rng, sz, rel):
    differ = [j for j, (a, b) in enumerate(H.PAIRS)
              if rel.kit_differs.iloc[j]]  # `rel` of the first group: same pair order for every group
    same = [j for j in range(len(H.PAIRS)) if j not in differ]
    out = []
    for label, idx in (("kit_differs", differ), ("kit_same", same)):
        d = h109_1(T, draws, n_eq, rng, sz, pair_idx=idx, label=label)
        out.append(d[d.null_type == "pseudo_clock"])
    return pd.concat(out, ignore_index=True)


def excess_bootstrap(T, n_eq, rng, sz):
    reps = {}
    obs = {}
    for g in GROUPS:
        gd = T.group_donors(g)
        n = min(n_eq, len(gd))
        ex = []
        for _ in range(sz["n_boot"]):
            D = rng.choice(gd, len(gd))
            if len(gd) > n_eq:
                D = rng.choice(D, n, replace=False)
            c = np.nanmean(abs(H.pair_corr(T, "pseudo", D)), axis=1)
            c0 = np.nanmean(abs(H.pair_corr(T, "clock", D)), axis=1)[0]
            ex.append(c0 - np.nanmedian(c))
        reps[g] = np.array(ex)
    return reps


def bootstrap_table(reps, h1):
    rows = []
    ex = h1[h1.null_type == "pseudo_clock"].set_index("group").excess
    for g, e in reps.items():
        rows.append(dict(contrast=g, estimate=ex[g], boot_mean=e.mean(), boot_se=e.std(ddof=1),
                         lo=np.nanquantile(e, 0.025), hi=np.nanquantile(e, 0.975), n_boot=len(e)))
    for a, b in (("vent", "fast"), ("H4", "fast"), ("vent", "H4"), ("H3", "fast")):
        d = reps[a] - reps[b]
        rows.append(dict(contrast=f"{a}-{b}", estimate=ex[a] - ex[b], boot_mean=d.mean(), boot_se=d.std(ddof=1),
                         lo=np.nanquantile(d, 0.025), hi=np.nanquantile(d, 0.975), n_boot=len(d)))
    return pd.DataFrame(rows)


# ---- non-temporal direction, host descriptives, circularity ----------------------------------------------------------------

def nontemporal(T):
    fast = T.group_donors("fast")
    Zr, Mr = H.tod_residual_states(T, fast)
    Mx = H.cross_cov_matrix(Zr, Mr)
    w, v = np.linalg.eigh(Mx)
    nt_h = float((np.arctan2(v[1, -1], v[0, -1]) * H.HOURS) % 12)
    rows = []
    Zraw, M, V = _raw_matrix(T)
    for g in GROUPS:
        gd = T.group_donors(g)
        Zc, m = _centre(Zraw, M, gd)
        Cm = H.cross_cov_matrix(Zc, m)
        Hm, _, _ = H.hermitian(Zc, m, V[gd])
        u, lead = H.host_frame(Hm)
        w_ = H.host_states(Zc, m, u)
        phi = np.angle(w_)
        ax = H.axis_angle_h(w_)
        d_gc = abs(((ax - H.GC_H % 12) + 6) % 12 - 6)
        rows.append(dict(group=g, n_donors=len(gd), nt_axis_h=nt_h, nt_eig_share=float(w[-1] / max(w.sum(), 1e-12)),
                         share_along_nt=H.share_along(Cm, nt_h), share_orth_nt=1 - H.share_along(Cm, nt_h),
                         share_along_gc=H.share_along(Cm, H.GC_H % 12), host_axis_h=ax,
                         host_axis_to_gc_h=d_gc, nt_axis_to_gc_h=abs(((nt_h - H.GC_H % 12) + 6) % 12 - 6),
                         leading_eig_share=lead, rayleigh_R=float(abs(np.mean(np.exp(1j * phi)))),
                         axial_R=float(abs(np.mean(np.exp(2j * phi)))),
                         nt_learned_in_fast=(g == "fast")))
    return pd.DataFrame(rows)


def circularity(T, rng, sz):
    Zraw, M, V = _raw_matrix(T)
    gd = T.group_donors("vent")
    Zc, m = _centre(Zraw, M, gd)
    obs = H.circ_ratio(Zc, m)
    c = float(np.nanmean(abs(H.pair_corr(T, "clock", gd)[0])))
    rows = [dict(model="observed", fraction=np.nan, rep=-1, ratio=obs, C_target=c)]
    for name, f in (("rotating", 0.0), ("mix25", 0.25), ("mix50", 0.5), ("mix75", 0.75), ("displacement", 1.0)):
        for r in range(sz["n_sim"]):
            zs = H.simulate_group(Zc, m, V[gd], c, f, H.GC_H, rng)
            rows.append(dict(model=name, fraction=f, rep=r, ratio=H.circ_ratio(zs, m), C_target=c))
    d = pd.DataFrame(rows)
    med = d[d.model != "observed"].groupby("fraction").ratio.median()
    d["best_fraction"] = float(med.index[np.argmin(abs(med - obs))])
    d["best_fraction_interp"] = float(np.interp(obs, med.values, med.index.values)) if np.all(np.diff(med.values) > 0) else np.nan
    return d


# ---- scout ----------------------------------------------------------------------------------------------------------

def scout(T, n_eq, rng, sz, out):
    rows = []
    for g in GROUPS:
        gd = T.group_donors(g)
        for t in BODY:
            rows.append(dict(group=g, tissue=t, n_donors=int((T.R[t][gd] >= 0).sum())))
    cnt = pd.DataFrame(rows)
    cnt.to_csv(out / "scout_counts.csv", index=False)
    print(cnt.pivot(index="tissue", columns="group", values="n_donors").to_string(), flush=True)
    rows = []
    for g in GROUPS:
        gd = T.group_donors(g)
        n = min(n_eq, len(gd))
        dr = H.draws_equal_n(gd, n, sz["n_sub"], rng)
        np_full = int(np.sum([np.isfinite(x) for x in H.pair_corr(T, "pseudo", gd)[0]]))
        np_eq = [int(np.sum(np.isfinite(H.pair_corr(T, "pseudo", D)[0]))) for D in dr]
        nul = H.coherence(T, "pseudo", dr)
        rows.append(dict(group=g, donors_ge2_tissues=len(gd), n_equal=n, pairs_full=np_full,
                         pairs_equal_n_min=min(np_eq), pairs_equal_n_median=int(np.median(np_eq)),
                         null_median=np.median(nul), null_sd=nul.std(ddof=1), null_q95=np.quantile(nul, 0.95),
                         null_q99=np.quantile(nul, 0.99), mde_p05=np.quantile(nul, 0.95) - np.median(nul),
                         mde_p01=np.quantile(nul, 0.99) - np.median(nul)))
    s = pd.DataFrame(rows)
    s.to_csv(out / "scout_summary.csv", index=False)
    print(s.round(4).to_string(index=False), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--scout", action="store_true", help="group sizes, pair counts and MDE only (no real clock C of "
                                                         "ventilated donors)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    sz = sizes(run["quick"])
    out = exploratory_dir(FOLDER, run["quick"])
    rng = np.random.default_rng(109 + args.seed)
    T = H.Tissues(sz["n_sets"], args.seed, args.subsample)
    T.build_states(np.random.default_rng(2109 + args.seed), sz["n_sets"], real=True)
    n_eq = len(T.group_donors("H4"))
    print({g: len(T.group_donors(g)) for g in GROUPS}, "equal n (A93 rule, Hardy 4):", n_eq, flush=True)
    if args.scout:
        scout(T, n_eq, rng, sz, out)
        pw = h109_2_power(T, rng, sz)
        pw.to_csv(out / "scout_h109_2_power.csv", index=False)
        print(pw.round(3).to_string(index=False))
        return
    draws = {g: H.draws_equal_n(T.group_donors(g), min(n_eq, len(T.group_donors(g))), sz["n_sub"], rng)
             for g in GROUPS}
    h1 = h109_1(T, draws, n_eq, rng, sz)
    h1.to_csv(out / "h109_1_null.csv", index=False)
    print(h1.round(4).to_string(index=False), flush=True)
    clock_plane_variance(T, draws).to_csv(out / "clock_plane_variance.csv", index=False)
    h2 = h109_2(T, rng, sz, n_eq)
    h2.to_csv(out / "h109_2_offsets.csv", index=False)
    print(h2.round(4).to_string(index=False), flush=True)
    rel = pair_relations(T, H.procedure_kits())
    rel.to_csv(out / "pair_relations.csv", index=False)
    pc = procedure_check(T, draws, n_eq, rng, sz, rel[rel.group == "fast"].reset_index(drop=True))
    pc.to_csv(out / "procedure_check.csv", index=False)
    reps = excess_bootstrap(T, n_eq, rng, sz)
    bt = bootstrap_table(reps, h1)
    bt.to_csv(out / "excess_bootstrap.csv", index=False)
    print(bt.round(4).to_string(index=False), flush=True)
    nt = nontemporal(T)
    nt.to_csv(out / "nontemporal_share.csv", index=False)
    print(nt.round(3).to_string(index=False), flush=True)
    cs_ = circularity(T, rng, sz)
    cs_.to_csv(out / "circularity_sims.csv", index=False)
    print(cs_.groupby("model", dropna=False).ratio.describe().round(3).to_string(), flush=True)
    h109_2_power(T, rng, sz).to_csv(out / "h109_2_power.csv", index=False)


if __name__ == "__main__":
    main()
