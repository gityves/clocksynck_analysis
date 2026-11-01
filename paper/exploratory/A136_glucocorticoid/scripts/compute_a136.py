"""A136 (registered 2026-11-01): is the loss of organ phase relations in brain-dead organ donors a drug effect?

Registered in paper/docs/registered_experiments.md, section A136 (draft: A136_draft.md). Reuses A132 / A132b
(`compute_brain_death`, `compute_posthoc_rho`) and the A109 states. Controlled access (dbGaP phs000424.v8): phenotypes,
time of death and per-donor scores stay in memory; outputs are aggregates only.

Stages: (0) A132b check (14 genes); (1) GR and hypoxia scores; (2) power for H136.2 (tertile-label permutation null,
written to power.csv before any observed H136.2 statistic); (3) tests H136.1-H136.3, companions, Holm, MHSTRDLT
sensitivity.

    python compute_a136.py [--subsample N] [--seed S]
Output (paper/exploratory/A136_glucocorticoid/output/): reproduction.csv, power.csv, scores_summary.csv, tests.csv,
companions.csv, sensitivity.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

A132 = Path(__file__).resolve().parents[2] / "A132_brain_death" / "scripts"
sys.path.insert(0, str(A132))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import compute_brain_death as B  # noqa: E402
import compute_posthoc_rho as R  # noqa: E402
from compute_brain_death import H, N, ST, cli, exploratory_dir, gtex  # noqa: E402
import clocksync as cs  # noqa: E402

FOLDER = "A136_glucocorticoid"
GR = ["FKBP5", "TSC22D3", "ZBTB16", "KLF15", "DDIT4"]
HYPOXIA = ["VEGFA", "ADM", "BNIP3", "PGK1", "SLC2A1"]
NOPER = [g for g in cs.gtex().genes if g not in ("PER1", "PER2")]
Z_MDE = 2.49  # registered: MDE = 2.49 x null SD


# ---- exposure scores ---------------------------------------------------------------------------------------------------

def tissue_scores(gene_sets: dict) -> dict:
    """{name: Series indexed by donor}: per tissue mean log CPM of the genes, residualised on ischaemic time, z-scored
    within tissue over all donors; donor score = median over the donor's tissues (the 12 A92 tissues)."""
    import anndata as ad
    from clocksync.genes import gene_columns

    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    per = {k: [] for k in gene_sets}
    try:
        o = a.obs
        keep = o.tissue.astype(str).isin(H.BODY).to_numpy()
        sub = o.iloc[np.flatnonzero(keep)]
        subj = sub.subject.astype(str).to_numpy()
        ok = ~pd.Series(list(zip(subj, sub.tissue.astype(str)))).duplicated().to_numpy()
        sub, subj = sub[ok], subj[ok]
        for t in H.BODY:
            m = (sub.tissue.astype(str) == t).to_numpy()
            ids = sub.index[m]
            pos = a.obs_names.get_indexer(ids)
            srt = np.argsort(pos)
            pos, ids = pos[srt], ids[srt]
            sd = subj[m][srt]
            X = a.X[pos]
            X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            lib = a.obs.library_size.to_numpy(float)[pos]
            isch = a.obs.ischemic_time.to_numpy(float)[pos]
            lc = np.log(X / lib[:, None] * 1e6 + 0.5)
            for k, genes in gene_sets.items():
                s = lc[:, gene_columns(a.var_names, genes, missing="raise")].mean(1)
                f = np.isfinite(isch) & np.isfinite(s)
                A = np.c_[np.ones(f.sum()), isch[f]]
                r = s[f] - A @ np.linalg.lstsq(A, s[f], rcond=None)[0]
                per[k].append(pd.Series((r - r.mean()) / r.std(ddof=1), index=sd[f]))
    finally:
        a.file.close()
    return {k: pd.concat(v).groupby(level=0).median() for k, v in per.items()}


# ---- statistics --------------------------------------------------------------------------------------------------------

def rho(Zraw, M, fast, g, mask):
    return R.rho(Zraw, M, fast, g, mask)


def split_idx(score: np.ndarray, d: np.ndarray):
    """Bottom and top tertile of donors d by score (equal n = floor(n / 3), by rank)."""
    o = d[np.argsort(score[d], kind="stable")]
    k = len(d) // 3
    return o[:k], o[-k:]


def tertile_stat(Zraw, M, fast, bot, top, mask):
    return rho(Zraw, M, fast, top, mask) - rho(Zraw, M, fast, bot, mask)


def perm_null(Zraw, M, fast, d, k, masks: dict, n_perm, seed):
    """Tertile-label permutation null of rho(top) - rho(bottom): random disjoint k-subsets of d; one null per mask
    (same permutations)."""
    rng = np.random.default_rng(seed)
    out = {n: [] for n in masks}
    for _ in range(n_perm):
        p = rng.permutation(d)
        for n, mk in masks.items():
            out[n].append(tertile_stat(Zraw, M, fast, p[:k], p[-k:], mk))
    return {n: np.asarray(v) for n, v in out.items()}


def a132b(Zraw, M, fast, od, pm, seed, n_perm):
    """A132b statistic (rho post-mortem minus organ donor, equal-n A132 draw, two-sided donor-label null)."""
    rng = np.random.default_rng(13200 + seed)
    ns = min(len(od), len(pm))
    rd = np.random.default_rng(132 + seed)
    od_ = rd.choice(od, ns, replace=False) if len(od) > ns else od
    pm_ = rd.choice(pm, ns, replace=False) if len(pm) > ns else pm
    mask = B.pair_mask(M, fast, od_, pm_)
    r_od, r_pm = rho(Zraw, M, fast, od_, mask), rho(Zraw, M, fast, pm_, mask)
    obs = r_pm - r_od
    pool, null = np.concatenate([pm_, od_]), []
    for _ in range(n_perm):
        p = rng.permutation(pool)
        null.append(rho(Zraw, M, fast, p[:ns], mask) - rho(Zraw, M, fast, p[ns:], mask))
    null = np.asarray(null)
    return dict(rho_postmortem=r_pm, rho_organ_donor=r_od, difference=obs, p_two_sided=N.pval(abs(obs), np.abs(null)),
                p_one_sided=N.pval(obs, null), null_sd=float(np.std(null, ddof=1)), n_each=ns, n_pairs=int(mask.sum()),
                n_fast=len(fast), n_perm=len(null))


def build(subsample, seed, genes=None):
    if genes is not None:
        H.WF = cs.gtex().subset(list(genes))
    T = H.Tissues(1, seed, subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    Zraw, M, _ = N._raw_matrix(T)
    return T, Zraw, M


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(run, 1.0)
    return adj


# ---- main --------------------------------------------------------------------------------------------------------------

def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    q = run["quick"]
    n_perm_h2, n_perm_h1, n_perm_b = (50, 500, 50) if q else (1000, 10_000, 1000)
    out = exploratory_dir(FOLDER, q)
    seed = args.seed

    T, Zraw, M = build(args.subsample, seed)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh, steroid = ph.COHORT.to_numpy(), ph.MHSTRDLT.to_numpy()
    h0, fast = T.group_donors("vent"), T.group_donors("fast")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]

    # (0) reproduction of A132b (14 genes)
    rep = pd.DataFrame([dict(analysis="A132b_14gene_reproduction", **a132b(Zraw, M, fast, od, pm, seed, n_perm_b))])
    rep.to_csv(out / "reproduction.csv", index=False)
    print(rep.round(4).to_string(index=False), flush=True)

    # (1) scores (per-donor values stay in memory)
    sc = tissue_scores({"GR": GR, "hypoxia": HYPOXIA})
    gr = sc["GR"].reindex(T.donors).to_numpy()
    hy = sc["hypoxia"].reindex(T.donors).to_numpy()
    summ = []
    for lab, d in (("organ_donor", od), ("ventilated_postmortem", pm), ("fast", fast)):
        summ.append(dict(group=lab, n=len(d), n_gr=int(np.isfinite(gr[d]).sum()), gr_median=np.nanmedian(gr[d]),
                         gr_q25=np.nanquantile(gr[d], .25), gr_q75=np.nanquantile(gr[d], .75),
                         hypoxia_median=np.nanmedian(hy[d]), n_steroid_longterm=int((steroid[d] == 1).sum())))
    pd.DataFrame(summ).to_csv(out / "scores_summary.csv", index=False)
    print(pd.DataFrame(summ).round(3).to_string(index=False), flush=True)
    od = od[np.isfinite(gr[od]) & np.isfinite(hy[od])]
    pm = pm[np.isfinite(gr[pm]) & np.isfinite(hy[pm])]

    # (2) power for H136.2 and its companions: null of the tertile-label permutation (no observed statistic)
    def masks_for(d, scores):
        res = {}
        for n, s in scores.items():
            bot, top = split_idx(s, d)
            res[n] = (B.pair_mask(M, fast, bot, top), len(bot))
        return res
    mk_od = masks_for(od, {"GR": gr, "hypoxia": hy})
    mk_pm = masks_for(pm, {"GR": gr})
    k_od, k_pm = mk_od["GR"][1], mk_pm["GR"][1]
    nulls_od = perm_null(Zraw, M, fast, od, k_od, {n: v[0] for n, v in mk_od.items()}, n_perm_h2, 13600 + seed)
    nulls_pm = perm_null(Zraw, M, fast, pm, k_pm, {"GR": mk_pm["GR"][0]}, n_perm_h2, 13601 + seed)
    pw = []
    for lab, nl, mk, k in (("H136.2_organ_donor_GR", nulls_od["GR"], mk_od["GR"][0], k_od),
                           ("companion_organ_donor_hypoxia", nulls_od["hypoxia"], mk_od["hypoxia"][0], k_od),
                           ("companion_ventilated_GR", nulls_pm["GR"], mk_pm["GR"][0], k_pm)):
        sd = float(np.std(nl, ddof=1))
        pw.append(dict(test=lab, n_per_tertile=k, n_pairs=int(mk.sum()), null_sd=sd, mde=Z_MDE * sd,
                       null_q05=float(np.quantile(nl, .05)), n_perm=len(nl),
                       min_p=1 / (1 + len(nl))))
    pd.DataFrame(pw).to_csv(out / "power.csv", index=False)
    print(pd.DataFrame(pw).round(4).to_string(index=False), flush=True)

    # (3) tests
    rows, comp, sens = [], [], []
    rng = np.random.default_rng(13602 + seed)

    def h1(o_, p_, sc_):
        obs = np.median(sc_[o_]) - np.median(sc_[p_])
        pool, n_o = np.concatenate([o_, p_]), len(o_)
        null = np.array([np.median(sc_[(pp := rng.permutation(pool))[:n_o]]) - np.median(sc_[pp[n_o:]])
                         for _ in range(n_perm_h1)])
        return obs, N.pval(obs, null), float(np.std(null, ddof=1))

    obs1, p1, sd1 = h1(od, pm, gr)
    rows.append(dict(test="H136.1_GR_organ_donor_minus_ventilated", statistic=obs1, p_one_sided=p1, null_sd=sd1,
                     n_a=len(od), n_b=len(pm), n_perm=n_perm_h1))

    def h2(d, s, nl_, mk_):
        bot, top = split_idx(s, d)
        r_t, r_b = rho(Zraw, M, fast, top, mk_), rho(Zraw, M, fast, bot, mk_)
        return dict(rho_top=r_t, rho_bottom=r_b, statistic=r_t - r_b, p_one_sided=N.pval(r_t - r_b, nl_, upper=False),
                    n_per_tertile=len(bot), n_pairs=int(mk_.sum()), null_sd=float(np.std(nl_, ddof=1)), n_perm=len(nl_))

    r2 = h2(od, gr, nulls_od["GR"], mk_od["GR"][0])
    rows.append(dict(test="H136.2_GR_tertile_rho_top_minus_bottom_organ_donors", **r2))
    comp.append(dict(test="companion_ventilated_GR_tertile", **h2(pm, gr, nulls_pm["GR"], mk_pm["GR"][0])))
    comp.append(dict(test="companion_organ_donor_hypoxia_tertile", **h2(od, hy, nulls_od["hypoxia"], mk_od["hypoxia"][0])))

    # H136.3: 12-gene waveform without PER1 / PER2, A132b statistic
    T12, Z12, M12 = build(args.subsample, seed, NOPER)
    assert (T12.donors == T.donors).all()
    h0_12, fast12 = T12.group_donors("vent"), T12.group_donors("fast")
    od12, pm12 = h0_12[coh[h0_12] == B.OPO], h0_12[coh[h0_12] == B.PM]
    b12 = a132b(Z12, M12, fast12, od12, pm12, seed, n_perm_b)
    rows.append(dict(test="H136.3_A132b_difference_no_PER1_PER2", statistic=b12["difference"],
                     p_one_sided=b12["p_one_sided"], p_two_sided=b12["p_two_sided"], null_sd=b12["null_sd"],
                     n_per_tertile=b12["n_each"], n_pairs=b12["n_pairs"], n_perm=b12["n_perm"],
                     rho_postmortem=b12["rho_postmortem"], rho_organ_donor=b12["rho_organ_donor"]))
    t = pd.DataFrame(rows)
    # Holm over H136.1-H136.3; H136.3 uses its two-sided p (A132b form; "stays" is read on p < 0.05)
    praw = [rows[0]["p_one_sided"], rows[1]["p_one_sided"], b12["p_two_sided"]]
    t["p_used"] = praw
    t["p_holm"] = holm(praw)
    t.to_csv(out / "tests.csv", index=False)
    pd.DataFrame(comp).to_csv(out / "companions.csv", index=False)
    print(t.round(4).to_string(index=False), "\n", pd.DataFrame(comp).round(4).to_string(index=False), flush=True)

    # MHSTRDLT sensitivity: donors on long-term steroids removed from the score and the groups
    keep = steroid != 1
    gr_s = tissue_scores_excl(keep_donors=T.donors[keep])
    od_s = od[keep[od]]
    pm_s = pm[keep[pm]]
    gs = gr_s.reindex(T.donors).to_numpy()
    od_s, pm_s = od_s[np.isfinite(gs[od_s])], pm_s[np.isfinite(gs[pm_s])]
    obs1s, p1s, sd1s = h1(od_s, pm_s, gs)
    bot, top = split_idx(gs, od_s)
    mk = B.pair_mask(M, fast, bot, top)
    nl = perm_null(Zraw, M, fast, od_s, len(bot), {"GR": mk}, n_perm_h2, 13603 + seed)["GR"]
    r2s = h2_generic(Zraw, M, fast, od_s, gs, nl, mk)
    sens.append(dict(test="H136.1_no_MHSTRDLT", statistic=obs1s, p_one_sided=p1s, null_sd=sd1s, n_a=len(od_s),
                     n_b=len(pm_s), n_removed=int((~keep[np.concatenate([od, pm])]).sum())))
    sens.append(dict(test="H136.2_no_MHSTRDLT", **r2s))
    pd.DataFrame(sens).to_csv(out / "sensitivity.csv", index=False)
    print(pd.DataFrame(sens).round(4).to_string(index=False), flush=True)


def tissue_scores_excl(keep_donors) -> pd.Series:
    """GR score recomputed on the donors in keep_donors only (residualisation, z-scoring and median)."""
    import anndata as ad
    from clocksync.genes import gene_columns

    kd = set(map(str, keep_donors))
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    per = []
    try:
        o = a.obs
        sub = o.iloc[np.flatnonzero(o.tissue.astype(str).isin(H.BODY).to_numpy())]
        subj = sub.subject.astype(str).to_numpy()
        ok = ~pd.Series(list(zip(subj, sub.tissue.astype(str)))).duplicated().to_numpy() & np.isin(subj, list(kd))
        sub, subj = sub[ok], subj[ok]
        for t in H.BODY:
            m = (sub.tissue.astype(str) == t).to_numpy()
            pos = a.obs_names.get_indexer(sub.index[m])
            srt = np.argsort(pos)
            pos, sd = pos[srt], subj[m][srt]
            X = a.X[pos]
            X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            lib = a.obs.library_size.to_numpy(float)[pos]
            isch = a.obs.ischemic_time.to_numpy(float)[pos]
            s = np.log(X / lib[:, None] * 1e6 + 0.5)[:, gene_columns(a.var_names, GR, missing="raise")].mean(1)
            f = np.isfinite(isch) & np.isfinite(s)
            A = np.c_[np.ones(f.sum()), isch[f]]
            r = s[f] - A @ np.linalg.lstsq(A, s[f], rcond=None)[0]
            per.append(pd.Series((r - r.mean()) / r.std(ddof=1), index=sd[f]))
    finally:
        a.file.close()
    return pd.concat(per).groupby(level=0).median()


def h2_generic(Zraw, M, fast, d, s, nl_, mk_):
    bot, top = split_idx(s, d)
    r_t, r_b = rho(Zraw, M, fast, top, mk_), rho(Zraw, M, fast, bot, mk_)
    return dict(rho_top=r_t, rho_bottom=r_b, statistic=r_t - r_b, p_one_sided=N.pval(r_t - r_b, nl_, upper=False),
                n_per_tertile=len(bot), n_pairs=int(mk_.sum()), null_sd=float(np.std(nl_, ddof=1)), n_perm=len(nl_))


if __name__ == "__main__":
    main()
