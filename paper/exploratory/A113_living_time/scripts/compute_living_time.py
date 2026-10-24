"""A113 (internal time of single living biopsies) and A130 (internal minus wall time against age): scoring.

Registered in paper/docs/registered_experiments.md (sections A113-A131 and the WP4 note). Fixed choices of the note: exon layer,
13 genes for muscle, adipose scored against ZT, lean + T2D pooled, both cohorts must pass, A130 pool = V25.5 + GSE104674 with
cohort as stratum, no sex descriptive.
  H113.1 per dataset (GSE108539 muscle, GSE104674 adipose): circular correlation (Jammalamadaka-SenGupta) of the single-sample
         clock-plane phase with recorded time (muscle: clock time; adipose: ZT), permutation p (1,000), pass: rho > 0.3, p < 0.05.
         Single sample = A35 arm (a): states of the GTEx projector (A28 E1), centred on the cohort mean. Permutation: the time
         labels are shuffled within each subject (subject = block; keeps subject offsets, identical schedules make a permutation of
         whole subjects uninformative); lean / T2D strata are respected because the shuffle is within subject.
         Two estimators are reported: JS (Jammalamadaka-SenGupta, A35 convention, referenced to the circular mean of each sample) and
         Fisher-Lee (invariant to the mean direction). JS is ill-defined when the recorded times are evenly spread (clocksync
         docstring): muscle times are a 4-hourly grid over 24 h, so the JS value of muscle is arbitrary; Fisher-Lee is the valid
         statistic there (choice made after seeing the JS result: post hoc, flagged).
  Secondary: median |error| after the common offset (A35 `score`); within-subject ordering (A35 arm c: states centred per
         subject; pooled rho, per-subject rho, same permutation); muscle window subsets (08-18 h, 09-17 h) with the A25 estimator on
         the sample states (single-tissue host, frame = GTEx frame, shift 0), Spearman of the tangent coordinate with time and the
         hospital-prior direction; the common offset (muscle).
  Muscle projector without ARNTL: the GTEx waveform is REFIT on the 13 genes (A28 recipe: `fit_waveform`, 15 tissues, fast deaths
         Hardy 1-2 with time of death, y = log((count + 0.5) / library size), tissue as group; all donors, no split); a sensitivity
         drops the ARNTL row of the published 14-gene waveform (`Waveform.subset`, projector renormalised).
         The 14-gene refit is compared with the published waveform as a check (waveform_refit.csv).
  H130.1 per subject internal minus wall time (V25.5: A25 estimator with the hospital prior, `compute_breast_timed.run`;
         adipose: circular mean over the subject's biopsies of single-sample phase minus ZT); cohort offset removed (stratum);
         within-cohort OLS slope on age (h/year); one-sided permutation (slope < 0) of age within cohort (1,000).
Controlled access: GTEx time of death used in memory for the waveform refit; only waveform parameters are written.

    python compute_living_time.py [--subsample N] [--seed S]
Outputs (paper/exploratory/A113_living_time/output/): a113_summary.csv, a113_subjects.csv, a113_windows.csv, a113_phases.csv,
waveform_refit.csv, a130_subjects.csv, a130_summary.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_1" / "scripts" / "windows"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402
from common import cli, data, gtex  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

import _internal_time as IT  # noqa: E402
import compute_breast_timed as BT  # noqa: E402

FOLDER = "A113_living_time"
W = 2 * np.pi / 24
CORE = list(cs.gtex().genes)
TIS = ("Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",
       "Esophagus - Mucosa", "Nerve - Tibial", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
       "Adipose - Subcutaneous", "Lung", "Pituitary", "Brain - Cerebellum", "Brain - Cortex", "Whole Blood")  # A28 pool
WINDOWS = {"08-18 h": (8, 18), "09-17 h": (9, 17)}
N_PERM = 1000


def refit(genes, tag):
    """A28 recipe on all fast-death donors of the 15 tissues; controlled time of death stays in memory."""
    pb, donors = gtex.load(tissues=TIS, hardy=(1, 2), genes=genes, depth="library_size")
    pb = pb.merge(donors[["tod"]], left_on="donor", right_index=True)
    y = np.log((pb[genes].to_numpy(float) + 0.5) / pb.depth.to_numpy(float)[:, None])
    return fit_waveform(y, pb.tod.to_numpy(float) % 24, genes, group=pb.compartment.to_numpy(), name=tag)


def waveforms():
    pub = cs.gtex()
    g13 = [g for g in CORE if g != "ARNTL"]
    w14, w13 = refit(CORE, "refit14"), refit(g13, "refit13")
    drop = pub.subset(g13)
    tab = [wf.table().reset_index(names="gene").assign(waveform=lab)
           for lab, wf in (("published14", pub), ("refit14", w14), ("refit13", w13), ("published_drop_ARNTL", drop))]
    return pub, w13, drop, pd.concat(tab)


def cmean(a):
    return np.angle(np.exp(1j * np.asarray(a)).mean())


def score(ph, tr):
    """Circular correlation, offset (h, mean error) and median |error| after the offset (h), as A35 `score`."""
    e = np.asarray(ph) - np.asarray(tr)
    m = cmean(e)
    return dict(rho=cs.stats.js_corr(ph, tr), rho_fl=cs.stats.fisher_lee(ph, tr), resultant=float(abs(np.exp(1j * e).mean())), offset_h=float(m / W), med_abs_err_h=float(np.median(np.abs(np.angle(np.exp(1j * (e - m)))) / W)))


def block_centre(A, blk):
    A = A.copy()
    for b in np.unique(blk):
        A[blk == b] -= A[blk == b].mean(0)
    return A


def phase_of(beta):
    return np.angle(beta[:, 0] + 1j * beta[:, 1])


def perm_within(ph, tr, blk, rng, n):
    """Null of rho (JS, Fisher-Lee): time labels shuffled within subject; (n, 2)."""
    idx = [np.flatnonzero(blk == b) for b in np.unique(blk)]
    out = np.empty((n, 2))
    for k in range(n):
        t = tr.copy()
        for i in idx:
            t[i] = tr[rng.permutation(i)]
        out[k] = cs.stats.js_corr(ph, t), cs.stats.fisher_lee(ph, t)
    return out


def sample_states(tab, genes, wf, counts):
    x = tab[genes].reset_index(drop=True)
    L = projector(x, wf, counts=counts)
    beta, _ = states(x, L, wf, depth=tab.depth.to_numpy(float) if counts else None)
    return beta


def windows(beta, hours, blk, tag):
    """A25 estimator on the sample states of a window subset (single-tissue host: frame = GTEx frame)."""
    rows = []
    for name, (lo, hi) in WINDOWS.items():
        for strat in (False, True):
            i = np.flatnonzero((hours >= lo) & (hours <= hi))
            f = beta[i]
            st = blk[i] if strat else None
            S = IT.host_cov(f, st, np.zeros(len(i)))
            e = IT.estimate(f, S, 0.0, st, pd.RangeIndex(len(i)))
            tw = IT.true_window(hours[i])
            t = hours[i]
            sp = sps.spearmanr(e["tangent_day"], t)
            rows.append(dict(dataset=tag, window=name, centred_within_subject=strat, n=len(i), n_times=len(np.unique(t)),
                             spearman_tangent_day=sp[0], p=sp[1], c_day_h=e["c_day"], c_true_h=tw["c_true"],
                             c_err_h=float(IT.cdiff(e["c_day"], tw["c_true"])), w_est_h=e["w"], w_true_h=tw["w_true"], r=e["r"],
                             direction_correct=bool(sp[0] > 0 and abs(IT.cdiff(e["c_day"], tw["c_true"])) < 6)))
    return rows


def score_dataset(tag, tab, time_col, genes, wf, counts, rng, n_perm, summ, subj, phs):
    beta = sample_states(tab, genes, wf, counts)
    blk = tab.block.astype(str).to_numpy()
    hours = tab[time_col].to_numpy(float)
    tr = W * hours
    for arm, b in (("a (single sample)", beta - beta.mean(0)), ("c (within-subject)", block_centre(beta, blk))):
        ph = phase_of(b)
        sc = score(ph, tr)
        null = perm_within(ph, tr, blk, rng, n_perm)
        p = (1 + (null[:, 0] >= sc["rho"]).sum()) / (n_perm + 1)
        p_fl = (1 + (null[:, 1] >= sc["rho_fl"]).sum()) / (n_perm + 1)
        per = pd.Series({k: cs.stats.js_corr(ph[blk == k], tr[blk == k]) for k in np.unique(blk) if (blk == k).sum() >= 3})
        per_fl = pd.Series({k: cs.stats.fisher_lee(ph[blk == k], tr[blk == k]) for k in per.index})
        summ.append(dict(dataset=tag, arm=arm, n_samples=len(tab), n_subjects=len(per), **sc, p_perm=p, p_perm_fl=p_fl, null_sd=null[:, 0].std(), null_sd_fl=null[:, 1].std(),
                         subject_rho_median=per.median(), subject_rho_frac_pos=(per > 0).mean(),
                         subject_rho_fl_median=per_fl.median(), subject_rho_fl_frac_pos=(per_fl > 0).mean(),
                         passes_H113_1_js=bool(arm.startswith("a") and sc["rho"] > 0.3 and p < 0.05),
                         passes_H113_1_fl=bool(arm.startswith("a") and sc["rho_fl"] > 0.3 and p_fl < 0.05)))
        subj += [dict(dataset=tag, arm=arm, subject=k, rho=v, rho_fl=per_fl[k]) for k, v in per.items()]
        if arm.startswith("a"):
            phs.append(pd.DataFrame(dict(dataset=tag, subject=blk, time_h=hours, phase_h=(ph / W) % 24)))
    return beta, blk, hours


def slope_within(age, y, coh):
    num = den = 0.0
    for c in np.unique(coh):
        i = coh == c
        a, v = age[i] - age[i].mean(), y[i] - y[i].mean()
        num += (a * v).sum()
        den += (a ** 2).sum()
    return num / den


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    limit_threads(1)
    out = exploratory_dir(FOLDER, run["quick"])
    n_perm = run["n_perm"] if run["quick"] else N_PERM
    rng = np.random.default_rng(args.seed)
    pub, w13, drop, wt = waveforms()
    wt.to_csv(out / "waveform_refit.csv", index=False)
    summ, subj, phs, wins = [], [], [], []

    # ---- muscle, GSE108539 (exon layer, array input, 13 genes) -------------------------------------------------
    mus, _ = data._builder("muscle_timed_gse108539").build("exon")
    if args.subsample:
        keep = rng.permutation(sorted(mus.block.unique()))[: args.subsample]
        mus = mus[mus.block.isin(keep)].reset_index(drop=True)
    g13 = list(w13.genes)
    for tag, wf in (("muscle GSE108539 (13-gene refit)", w13), ("muscle GSE108539 (sens: ARNTL row dropped)", drop)):
        beta, blk, hours = score_dataset(tag, mus, "clock_h", g13, wf, False, rng, n_perm, summ, subj, phs)
        wins += windows(beta, hours, blk, tag)

    # ---- adipose, GSE104674 (counts, 14 genes, ZT) -----------------------------------------------------------------
    adi, adi_donors = data._builder("adipose_timed_gse104674").build()
    if args.subsample:
        keep = rng.permutation(sorted(adi.block.unique()))[: args.subsample]
        adi = adi[adi.block.isin(keep)].reset_index(drop=True)
    beta_a, blk_a, zt = score_dataset("adipose GSE104674", adi, "zt_h", CORE, pub, True, rng, n_perm, summ, subj, phs)
    wins.append(dict(dataset="adipose GSE104674", window="08-18 h, 09-17 h", note="not applicable: 4 timepoints (ZT 15.5, 0.25, 5.75, 11.25)"))

    summ = pd.DataFrame(summ)
    summ.to_csv(out / "a113_summary.csv", index=False)
    pd.DataFrame(subj).to_csv(out / "a113_subjects.csv", index=False)
    pd.concat(phs).to_csv(out / "a113_phases.csv", index=False)
    pd.DataFrame(wins).to_csv(out / "a113_windows.csv", index=False)

    # ---- A130 -----------------------------------------------------------------------------------------------------
    ph_a = (phase_of(beta_a - beta_a.mean(0)) / W) % 24
    d_a = pd.Series(IT.cdiff(ph_a, zt)).groupby(blk_a).apply(lambda s: cmean(W * s) / W)
    sub_a = pd.DataFrame(dict(cohort="GSE104674", subject=d_a.index, delta_raw_h=d_a.to_numpy(),
                              age=adi_donors.loc[d_a.index, "age"].to_numpy(float), sex=adi_donors.loc[d_a.index, "sex"].to_numpy(),
                              group=adi_donors.loc[d_a.index, "status"].to_numpy()))
    b = data._builder("breast_timed_gse233242")
    pb, donors = b.build()
    if args.subsample:
        keep = rng.permutation(sorted(donors.index))[: max(args.subsample, 10)]
        pb, donors = pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[keep]
    _, e, host, _ = BT.run(pb, donors, n_boot=0, n_perm=1, seed=args.seed)
    sub_b = pd.DataFrame(dict(cohort="V25.5 GSE233242", subject=host.index,
                              delta_raw_h=IT.cdiff(e["donor_day_h"].to_numpy(), donors.resection_h.reindex(host.index).to_numpy()),
                              age=donors.age.reindex(host.index).to_numpy(float), sex="F", group=""))
    sub = pd.concat([sub_a, sub_b], ignore_index=True)
    off = {c: cmean(W * g.delta_raw_h.to_numpy()) / W for c, g in sub.groupby("cohort")}
    sub["delta_h"] = [float(IT.cdiff(v, off[c])) for c, v in zip(sub.cohort, sub.delta_raw_h)]
    sub.to_csv(out / "a130_subjects.csv", index=False)
    age, y, coh = sub.age.to_numpy(float), sub.delta_h.to_numpy(float), sub.cohort.to_numpy()
    obs = slope_within(age, y, coh)
    idx = [np.flatnonzero(coh == c) for c in np.unique(coh)]
    null = np.empty(n_perm)
    for k in range(n_perm):
        a = age.copy()
        for i in idx:
            a[i] = age[rng.permutation(i)]
        null[k] = slope_within(a, y, coh)
    p = (1 + (null <= obs).sum()) / (n_perm + 1)
    feas = bool(len(sub) >= 30 and np.ptp(age) >= 25)
    srows = [dict(scope="pooled (cohort as stratum)", n=len(sub), age_min=age.min(), age_max=age.max(), age_range=np.ptp(age),
                  slope_h_per_year=obs, slope_h_per_decade=10 * obs, p_one_sided=p, null_sd=null.std(), feasible=feas,
                  passes_H130_1=bool(feas and obs < 0 and p < 0.05))]
    for c, g in sub.groupby("cohort"):
        s = sps.linregress(g.age, g.delta_h)
        srows.append(dict(scope=c, n=len(g), age_min=g.age.min(), age_max=g.age.max(), age_range=np.ptp(g.age), slope_h_per_year=s.slope,
                          slope_h_per_decade=10 * s.slope, cohort_offset_h=off[c], spearman_age=sps.spearmanr(g.age, g.delta_h)[0]))
    ss = pd.DataFrame(srows)
    ss.to_csv(out / "a130_summary.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(wt.round(3).to_string(index=False))
    print(summ.round(3).to_string(index=False))
    print(pd.DataFrame(wins).round(3).to_string(index=False))
    print(ss.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
