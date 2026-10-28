"""A123c post hoc diagnostic (descriptive; PI-approved; it explains the H123.5 result and does not change the registered verdict).

Why does the stored A123 clocksync offset estimate agree with truth (iii) at r = 0.679 (H123.5 fail) but better (0.807) with a truth
built from a random half of the out-of-window donors?
 1. residuals.csv: per tissue, estimate - truth (iii), estimate - truth (i), truth (i) - truth (iii) offsets (h, after the centring of
    `compare`), n truth (iii) donors, truth (iii) amplitude and cosinor R2 of the clock genes, and leave-one-tissue-out circular r.
 2. bootstrap_truth.csv: noise ceiling of truth (iii). 200 donor bootstraps of the out-of-window donors (seed 1231), truth rebuilt
    each time with `compute_truth_iii.truth_from`; r(bootstrap truth, truth iii), r(stored estimate, bootstrap truth), r(truth i,
    bootstrap truth). split_half.csv: 200 random balanced half-splits (`split_donors`): r between the two half truths (min 10 donors),
    r(stored estimate, half truth), the rank of the registered split (seed 123), and the Spearman-Brown / attenuation ceiling
    (descriptive estimates only).
 3. construction_checks.csv: whether `truth_from` on all out-of-window donors differs from the half-2 truth in anything but the donors
    (size factors / gene set of `log2_normalised`, samples per tissue, reference tissue, centring set).
The stored estimate and truth (i) are read from output/organ_map.csv and output/truth_i.csv (A123 outputs; the estimate regeneration
was checked against them in compute_truth_iii.py). Controlled access: only per-tissue and summary aggregates are written; donor
identities, times and splits stay in memory.

    python diag_truth_iii.py [--subsample N] [--reps R]     (check mode: N tissues, R replicates, output under diag_check/)
"""

from __future__ import annotations

import multiprocessing as mp
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import compute_truth_iii as C  # noqa: E402  (A123 modules, truth_from, split_donors)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import cli  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.data import in_window, log2_normalised  # noqa: E402
from common.gtex import LABELLED  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

A, B, CLOCK, REF, SEED = C.A, C.B, C.CLOCK, C.REF, C.SEED
BOOT_SEED, N_REP, MIN_HALF = 1231, 200, 10
G = {}  # shared (forked) data: counts, tod, oow, normalised clock-gene tables


def r_off(est: pd.DataFrame, tru: pd.DataFrame) -> float:
    """Circular correlation of the centred offsets over the tissues both have (as `compare`)."""
    te = ev.tissue_errors(est, tru)
    return ev.circ_corr(te.offset_h, te.offset_h_true) if len(te) > 3 else np.nan


def quant(x) -> dict:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    q = np.quantile(x, [.025, .25, .5, .75, .975])
    return dict(n=len(x), mean=x.mean(), sd=x.std(ddof=1), q025=q[0], q25=q[1], median=q[2], q75=q[3], q975=q[4], min=x.min(), max=x.max())


def _truth(donors, own_norm=True, min_n=MIN_HALF):
    """Truth from per-tissue donor lists (duplicates allowed). own_norm=False: size factors from all out-of-window samples."""
    if own_norm:
        t, n, _ = C.truth_from(G["counts"], G["tod"], donors, min_n)
        return t[["amplitude", "offset_h"]]
    Cf = {t: ev.gene_rhythms(G["Lall"][t].loc[d], G["tod"][t].loc[d], CLOCK) for t, d in donors.items() if len(d) >= min_n}
    return ev.csvd_tissues(pd.DataFrame(Cf), REF)[["amplitude", "offset_h"]]


def _boot(i):
    rng = np.random.default_rng(np.random.SeedSequence(BOOT_SEED).spawn(N_REP)[i])
    alld = np.array(sorted(G["all_donors"]))
    pick = alld[rng.integers(0, len(alld), len(alld))]
    mult = pd.Series(pick).value_counts()
    donors = {t: [d for d in G["oow"][t] for _ in range(int(mult.get(d, 0)))] for t in G["tissues"]}
    return i, _truth(donors, min_n=C.MIN_TRUTH)  # the truth (iii) rule: >= 20 donors (with multiplicity)


def _split(i):
    half = C.split_donors(G["oow"], 10_000 + BOOT_SEED + i)
    out = {}
    for h in (1, 2):
        d = {t: [x for x in G["oow"][t] if half[x] == h] for t in G["tissues"]}
        out[h] = (_truth(d), _truth(d, own_norm=False))
    return i, out


def main(argv=None):
    global N_REP
    ap = cli.parser(__doc__)
    ap.add_argument("--reps", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    reps = args.reps or (4 if quick else N_REP)
    N_REP = reps
    out = C.exploratory_dir(C.AMENDMENT, False) / "a123c" / ("diag_check" if quick else "diag")
    out.mkdir(parents=True, exist_ok=True)
    base = C.exploratory_dir(C.AMENDMENT, False)
    t0 = time.time()

    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    o["tissue"] = o.tissue.astype(str)
    a.file.close()
    cnt = A.window_tissues(o)
    tissues = list(cnt.index[:args.subsample]) if quick else list(cnt.index)
    if REF not in tissues:
        tissues = [REF] + tissues[:-1]
    counts, tod_of, oow, Lall = {}, {}, {}, {}
    for t in tissues:
        c, tod = B.S.load([t])[t]
        tod = tod.dropna()
        c = c.loc[tod.index]
        oo = tod.index[~in_window(tod)]
        tod_of[t], oow[t], counts[t] = tod, list(oo), c.loc[oo]
        L = log2_normalised(counts[t])
        Lall[t] = L[[g for g in CLOCK if g in L.columns]]
        print(t, len(oo), f"{time.time() - t0:.0f} s", flush=True)
    G.update(counts=counts, tod=tod_of, oow=oow, tissues=tissues, Lall=Lall, all_donors={d for t in tissues for d in oow[t]})

    om = pd.read_csv(base / "organ_map.csv", index_col=0)
    est = pd.DataFrame({"amplitude": 2 ** om.est_log2_amp, "offset_h": om.est_offset_h}).loc[[t for t in tissues if t in om.index]]
    ti = pd.read_csv(base / "truth_i.csv", index_col=0)[["amplitude", "offset_h"]]
    t3, n3, dropped = C.truth_from(counts, tod_of, oow, C.MIN_TRUTH)
    kept = list(t3.index)
    t3s = t3[["amplitude", "offset_h"]]
    print(f"truth (iii): {len(kept)} tissues; r(est, t3) = {r_off(est, t3s):.4f}; {time.time() - t0:.0f} s", flush=True)
    ref_r = pd.read_csv(base / "a123c" / "tests.csv").set_index("row")
    if not quick:
        assert abs(r_off(est, t3s) - ref_r.loc["H123.5", "circ_corr_offset"]) < 5e-4, "r(est, truth iii) differs from tests.csv"
        td = pd.read_csv(base / "a123c" / "truth_iii.csv", index_col=0)
        tc = ev.centred(t3s, kept)
        assert np.abs(ev.wrap(tc.offset_h - td.offset_h_centred.reindex(kept))).max() < 1e-3, "truth_from differs from truth_iii.csv"

    # ---- 1. per-tissue residuals ------------------------------------------------------------------------------------------
    k = est.dropna().index.intersection(t3s.index)
    e_c, t_c, i_c = ev.centred(est, k), ev.centred(t3s, k), ev.centred(ti, k.intersection(ti.index))
    ti_k = ti.loc[k.intersection(ti.index)]
    i_c = ev.centred(ti_k, k)
    wr = lambda x: pd.Series(ev.wrap(x), index=x.index)  # noqa: E731
    res = pd.DataFrame(index=pd.Index(k, name="tissue"))
    res["n_truth_iii_donors"] = n3.reindex(k).astype(int)
    res["est_offset_h"], res["truth_iii_offset_h"], res["truth_i_offset_h"] = e_c.offset_h, t_c.offset_h, i_c.offset_h
    res["est_minus_truth_iii_h"] = wr(e_c.offset_h - t_c.offset_h)
    res["est_minus_truth_i_h"] = wr(e_c.offset_h - i_c.offset_h)
    res["truth_i_minus_truth_iii_h"] = wr(i_c.offset_h - t_c.offset_h)
    res["truth_iii_amplitude"] = t3.amplitude.reindex(k)
    res["truth_iii_log2_amp_centred"] = t_c.log2_amp
    r2 = {t: ev.cosinor(log2_normalised(counts[t].loc[oow[t]]), tod_of[t].loc[oow[t]], CLOCK).r2 for t in k}
    res["truth_iii_mean_clock_r2"] = pd.Series({t: float(np.mean(v)) for t, v in r2.items()})
    res["truth_iii_median_clock_r2"] = pd.Series({t: float(np.median(v)) for t, v in r2.items()})
    res["n_window_donors"] = cnt.reindex(k).astype(int)
    r_full = r_off(est, t3s)
    r_i = r_off(est, ti)
    res["loto_r_est_truth_iii"] = [r_off(est.drop(index=t), t3s.drop(index=t)) for t in k]
    res["loto_delta_r"] = res.loto_r_est_truth_iii - r_full
    res["loto_r_est_truth_i"] = [r_off(est.drop(index=t), ti.drop(index=t, errors="ignore")) for t in k]
    res["loto_r_truth_i_truth_iii"] = [r_off(ti.drop(index=t, errors="ignore"), t3s.drop(index=t)) for t in k]
    # registered split (seed 123) half-2 truth, min 10: the residual that explains the 0.807
    half = C.split_donors(oow, SEED)
    d2 = {t: [d for d in oow[t] if half[d] == 2] for t in tissues}
    t_h2, n_h2, _ = C.truth_from(counts, tod_of, d2, MIN_HALF)
    r_h2 = r_off(est, t_h2[["amplitude", "offset_h"]])
    kk = k.intersection(t_h2.index)
    h2c = ev.centred(t_h2[["amplitude", "offset_h"]], kk)
    ec2 = ev.centred(est, kk)
    res["n_half2_donors_seed123"] = n_h2.reindex(k)
    res["est_minus_half2_truth_h"] = wr(ec2.offset_h - h2c.offset_h).reindex(k)
    res["half2_minus_truth_iii_h"] = wr(h2c.offset_h - t_c.offset_h.reindex(kk)).reindex(k)
    res["abs_err_truth_iii_h"] = res.est_minus_truth_iii_h.abs()
    res["abs_err_half2_h"] = res.est_minus_half2_truth_h.abs()
    res["sin2_err_truth_iii"] = np.sin(2 * np.pi / 24 * res.est_minus_truth_iii_h) ** 2
    res = res.sort_values("loto_delta_r", ascending=False)
    res.round(4).to_csv(out / "residuals.csv")
    print(f"residuals done {time.time() - t0:.0f} s; r_full {r_full:.4f} r_half2(seed123) {r_h2:.4f}", flush=True)

    # ---- 2. bootstrap noise ceiling + split half -----------------------------------------------------------------------------
    jobs = min(args.jobs, 4)
    ctx = mp.get_context("fork")
    with ctx.Pool(jobs) as pool:
        bt = {}
        for n, (i, t) in enumerate(pool.imap_unordered(_boot, range(reps))):
            bt[i] = t
            if n % 10 == 0:
                print(f"bootstrap {n + 1}/{reps} {time.time() - t0:.0f} s", flush=True)
        sp = {}
        for n, (i, v) in enumerate(pool.imap_unordered(_split, range(reps))):
            sp[i] = v
            if n % 10 == 0:
                print(f"split {n + 1}/{reps} {time.time() - t0:.0f} s", flush=True)
    rb = pd.DataFrame({"r_boot_vs_truth_iii": [r_off(bt[i], t3s) for i in range(reps)],
                       "r_est_vs_boot": [r_off(est, bt[i]) for i in range(reps)],
                       "r_truth_i_vs_boot": [r_off(ti, bt[i]) for i in range(reps)],
                       "n_tissues_boot": [len(bt[i]) for i in range(reps)],
                       "n_common_tissues_boot_vs_truth_iii": [len(bt[i].index.intersection(t3s.index)) for i in range(reps)]})
    rows = []
    for c in ("r_boot_vs_truth_iii", "r_est_vs_boot", "r_truth_i_vs_boot", "n_common_tissues_boot_vs_truth_iii"):
        rows.append(dict(stat=c, **quant(rb[c])))
    ref_est, ref_i = r_full, r_off(ti, t3s)
    # bootstrap-percentile description of the observed r(est, truth iii) against r(est, boot truth)
    rows.append(dict(stat="observed_r_est_truth_iii", mean=ref_est))
    rows.append(dict(stat="observed_r_truth_i_truth_iii", mean=ref_i))
    rows.append(dict(stat="observed_r_est_truth_i", mean=r_i))
    rows.append(dict(stat="frac_boot_r_est_ge_observed_r_est_truth_iii", mean=float((rb.r_est_vs_boot >= ref_est).mean())))
    rows.append(dict(stat="frac_boot_r_est_ge_0.7", mean=float((rb.r_est_vs_boot >= 0.7).mean())))
    bs = pd.DataFrame(rows)
    bs.round(4).to_csv(out / "bootstrap_truth.csv", index=False)

    s = {c: [] for c in ("r_half1_half2", "r_est_half1", "r_est_half2", "r_est_half1_commonnorm", "r_est_half2_commonnorm",
                         "r_half1_half2_commonnorm", "r_truthi_half2", "n_tissues_half2")}
    for i in range(reps):
        (h1, h1n), (h2, h2n) = sp[i][1], sp[i][2]
        s["r_half1_half2"].append(r_off(h1, h2))
        s["r_half1_half2_commonnorm"].append(r_off(h1n, h2n))
        s["r_est_half1"].append(r_off(est, h1))
        s["r_est_half2"].append(r_off(est, h2))
        s["r_est_half1_commonnorm"].append(r_off(est, h1n))
        s["r_est_half2_commonnorm"].append(r_off(est, h2n))
        s["r_truthi_half2"].append(r_off(ti, h2))
        s["n_tissues_half2"].append(len(h2))
    s = pd.DataFrame(s)
    rows = [dict(stat=c, **quant(s[c])) for c in s]
    both = pd.concat([s.r_est_half1, s.r_est_half2])
    rows.append(dict(stat="r_est_either_half", **quant(both)))
    rows.append(dict(stat="registered_split_seed123_r_est_half2_min10", mean=r_h2,
                     note=f"fraction of the {len(both)} random half truths with r(est, half) below it: {(both < r_h2).mean():.3f}"))
    rows.append(dict(stat="frac_r_est_half_ge_registered_0.807", mean=float((both >= r_h2).mean())))
    rows.append(dict(stat="frac_r_est_half_ge_0.7", mean=float((both >= 0.7).mean())))
    rows.append(dict(stat="frac_r_est_half_ge_observed_r_est_truth_iii", mean=float((both >= ref_est).mean())))
    rows.append(dict(stat="frac_r_est_half_below_observed_r_est_truth_iii", mean=float((both < ref_est).mean())))
    # ceiling (descriptive): Spearman-Brown of the split-half r, and the bootstrap r as an alternative reliability proxy
    r_sh = s.r_half1_half2.mean()
    rel_sb = 2 * r_sh / (1 + r_sh)
    rb_m = rb.r_boot_vs_truth_iii.mean()
    for nm, rel in (("spearman_brown_of_mean_split_half_r", rel_sb), ("mean_bootstrap_r_vs_truth_iii", rb_m)):
        rows.append(dict(stat=f"ceiling_reliability[{nm}]", mean=rel))
        rows.append(dict(stat=f"ceiling_max_r_with_perfect_estimate[{nm}]", mean=float(np.sqrt(rel))))
        rows.append(dict(stat=f"disattenuated_r_est_truth_iii[{nm}] (truth noise only)", mean=float(ref_est / np.sqrt(rel))))
        rows.append(dict(stat=f"disattenuated_r_truth_i_truth_iii[{nm}] (truth noise only)", mean=float(ref_i / np.sqrt(rel))))
    rows[-1]["note"] = "descriptive estimates (linear-correlation attenuation applied to circular r); not a test"
    sh = pd.DataFrame(rows)
    sh.round(4).to_csv(out / "split_half.csv", index=False)
    print(f"bootstrap / split done {time.time() - t0:.0f} s", flush=True)

    # ---- 3. construction checks on the registered split (seed 123) ---------------------------------------------------------------
    chk = []
    full_n = n3.reindex(kk)
    chk.append(dict(check="truth_from(all out-of-window) reproduces truth_iii.csv", result="asserted" if not quick else "not run"))
    chk.append(dict(check="reference tissue", full=REF, half2=REF, half2_has_ref=bool(REF in t_h2.index),
                    n_ref_full=int(n3[REF]), n_ref_half2=int(n_h2[REF])))
    chk.append(dict(check="tissue sets", n_full=len(t3s), n_half2=len(t_h2), identical_centring_set=bool(set(t3s.index) == set(t_h2.index)),
                    n_common=len(kk), full_only=";".join(sorted(set(t3s.index) - set(t_h2.index))),
                    half2_only=";".join(sorted(set(t_h2.index) - set(t3s.index)))))
    chk.append(dict(check="donor filters", result="truth_from applies none beyond the donor list: same oow lists, same dropna on tod, "
                                                   "one sample per donor and tissue; half-2 donors are a subset"))
    chk.append(dict(check="samples per tissue (half2 / full)", min=float((n_h2.reindex(kk) / full_n).min()),
                    median=float((n_h2.reindex(kk) / full_n).median()), max=float((n_h2.reindex(kk) / full_n).max())))
    # size factors / gene set: truth from half-2 donors with size factors from all out-of-window samples
    t_h2n = _truth(d2, own_norm=False)
    chk.append(dict(check="size factors: own half-2 vs all out-of-window samples (seed 123 half 2)", r_own=r_h2,
                    r_common=r_off(est, t_h2n),
                    max_abs_offset_diff_h=float(np.abs(ev.wrap(ev.centred(t_h2n, kk).offset_h - h2c.offset_h)).max())))
    gs = {t: int((counts[t] > 0).all().sum()) for t in kk}
    gs2 = {t: int((counts[t].loc[d2[t]] > 0).all().sum()) for t in kk}
    chk.append(dict(check="genes positive in every sample (size-factor gene set), half2 / full", median=float(np.median([gs2[t] / gs[t] for t in kk])),
                    min=float(min(gs2[t] / gs[t] for t in kk)), max=float(max(gs2[t] / gs[t] for t in kk))))
    # estimate overlap: waveform cross-fit uses all out-of-window donors, the half truth a subset of them
    chk.append(dict(check="registered split: truth_from on half 2 with min 10 vs stored tests.csv",
                    r_here=r_h2, r_tests_csv=float(ref_r.loc["stored_clocksync_vs_truth_iii_half2_min10", "circ_corr_offset"])))
    chk.append(dict(check="same 39 tissues, full truth vs half-2 truth: est error (h) circular rms", full=float(np.sqrt((res.est_minus_truth_iii_h ** 2).mean())),
                    half2=float(np.sqrt((res.est_minus_half2_truth_h.dropna() ** 2).mean()))))
    pd.DataFrame(chk).to_csv(out / "construction_checks.csv", index=False)
    rb.round(4).to_csv(out / "bootstrap_truth_replicates_r.csv", index=False)  # only r values (200 rows), no donor data
    s.round(4).to_csv(out / "split_half_replicates_r.csv", index=False)

    print(res[["n_truth_iii_donors", "est_minus_truth_iii_h", "est_minus_truth_i_h", "truth_i_minus_truth_iii_h", "truth_iii_amplitude",
               "truth_iii_mean_clock_r2", "loto_delta_r", "est_minus_half2_truth_h"]].round(3).to_string())
    print(bs.round(3).to_string())
    print(sh.round(3).to_string())
    print(pd.DataFrame(chk).to_string())
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
