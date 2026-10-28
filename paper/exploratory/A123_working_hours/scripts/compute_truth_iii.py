"""A123c: the A123 organ map against a truth that shares no donor with the test (registered in paper/docs/registered_experiments.md).

Truth (iii): per-tissue cosinor of the 14 clock genes on the recorded time of GTEx v10 fast deaths (Hardy 1-2) who died OUTSIDE
08-18 h, rank-1 complex SVD (`ev.csvd_tissues`, reference tissue Muscle - Skeletal, as truth (i)), for the A123 tissues (those with
>= 20 window donors) that have >= 20 such donors; the others are dropped and counted. Expression is log2-normalised on the
out-of-window samples of the tissue (A123 normalised truth (i) on all fast deaths, the estimates on the window samples).
The truth donors and the test (window) donors are disjoint by construction; asserted.
Estimates (not refitted in substance): the A123 clocksync estimate is regenerated in memory with the A123 calls (`B.cross_fit_waveform`
on the donors outside the window, `B.clocksync_estimates`, seed 123; the donor bootstrap only adds intervals, so n_boot = 0) and
checked against `organ_map.csv`.
H123.5 offset circular correlation > 0.7 and tissue-permutation p < 0.05; H123.6 Pearson of log2 relative amplitudes > 0.7 (A123
`compare`, 10,000 permutations; its tissue bootstrap is the function's own 1,000 resamples).
Descriptive: truth (i) vs truth (iii); the A123b methods (CHIRAL, CoFE genome, CoFE clock, CYCLOPS, timed oracle) in the window design
from their cached orderings (`compute_methods.assemble`) vs truth (iii); the split sensitivity (out-of-window donors split at random
into halves, seed 123, balanced within each tissue; waveform from half 1, clocksync refitted on the window donors, scored against a
truth (iii) from half 2 alone, >= 20 half-2 donors per tissue, and a >= 10 variant); the prior-only baseline of the per-donor hour.
Controlled access: only per-tissue aggregates and summary statistics are written; donor times, hours and splits stay in memory.

    python compute_truth_iii.py [--subsample N] [--perm P]
`--subsample N`: check mode, the N tissues with most window donors, 200 permutations, output under output/subsample/a123c/.
Outputs (exploratory/A123_working_hours/output/a123c/): tests.csv, truth_iii.csv, README.md.
"""

from __future__ import annotations

import importlib.util
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
_sp = importlib.util.spec_from_file_location("a123_wh", HERE.parent / "compute_working_hours.py")
A = importlib.util.module_from_spec(_sp)
sys.modules["a123_wh"] = A
_sp.loader.exec_module(A)
_sm = importlib.util.spec_from_file_location("a123b_methods", HERE.parent / "compute_methods.py")
M = importlib.util.module_from_spec(_sm)
sys.modules["a123b_methods"] = M
_sm.loader.exec_module(M)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.data import in_window, log2_normalised  # noqa: E402
from common.gtex import LABELLED  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

B, T, CLOCK, REF, SEED, W = A.B, A.T, A.CLOCK, A.REF, A.SEED, A.W
MIN_TRUTH = 20
AMENDMENT = A.AMENDMENT
METHODS = ("chiral", "cofe_genome", "cofe_clock", "cyclops")


def truth_from(counts: dict, tod: dict, donors: dict, min_n: int):
    """Truth from `donors[t]` (recorded times in `tod[t]`), tissues with >= min_n donors; returns (truth, n per tissue, dropped)."""
    C, n = {}, {}
    for t, d in donors.items():
        n[t] = len(d)
        if len(d) >= min_n:
            C[t] = ev.gene_rhythms(log2_normalised(counts[t].loc[d]), tod[t].loc[d], CLOCK)
    if REF not in C:
        raise RuntimeError(f"reference tissue {REF} has < {min_n} donors")
    return ev.csvd_tissues(pd.DataFrame(C), REF), pd.Series(n), [t for t in donors if t not in C]


def split_donors(oow: dict, seed: int) -> dict:
    """Random halves of the out-of-window donors, balanced within each tissue: tissues in ascending size, each tissue's still
    unassigned donors are shuffled and sent to half 1 until the tissue has round(n/2) donors there, the rest to half 2."""
    rng = np.random.default_rng(seed)
    half = {}
    for t in sorted(oow, key=lambda t: (len(oow[t]), t)):
        d = list(oow[t])
        have1 = sum(half.get(x) == 1 for x in d)
        free = [x for x in d if x not in half]
        free = [free[i] for i in rng.permutation(len(free))]
        k = int(np.clip(round(len(d) / 2) - have1, 0, len(free)))
        for j, x in enumerate(free):
            half[x] = 1 if j < k else 2
    return half


def circ_median_offset(pred_h: np.ndarray, true_h: np.ndarray) -> float:
    """Common offset (h, 0.05 h grid) minimising the median absolute circular error of pred + offset against true."""
    grid = np.arange(0, 24, 0.05)
    err = [np.median(np.abs(T.cdiff(pred_h + g, true_h))) for g in grid]
    return float(grid[int(np.argmin(err))])


def score_row(name, kind, est, tru, n_perm, **extra):
    c = A.compare(est, tru, n_perm, SEED)
    return dict(row=name, kind=kind, **extra, **c)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_perm = args.perm if args.perm is not None else 200 if quick else 10000
    out = exploratory_dir(AMENDMENT, quick) / "a123c"
    out.mkdir(parents=True, exist_ok=True)
    runs = exploratory_dir(AMENDMENT, False) / "a123b" / "runs"
    t0 = time.time()

    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    o["tissue"] = o.tissue.astype(str)
    a.file.close()
    cnt = A.window_tissues(o)
    tissues = list(cnt.index[:args.subsample]) if quick else list(cnt.index)
    if REF not in tissues:
        tissues = [REF] + tissues[:-1]
    print(f"{len(tissues)} A123 tissues (>= {A.MIN_WINDOW} window donors)", flush=True)

    # ---- data: counts of the out-of-window donors (truth), clock table and window log2 (estimates) ----------------------------
    oow_counts, tod_of, oow, rows, LW, TODW = {}, {}, {}, [], {}, {}
    for t in tissues:
        counts, tod = B.S.load([t])[t]
        tod = tod.dropna()
        counts = counts.loc[tod.index]
        win = tod.index[in_window(tod)]
        out_ = tod.index[~in_window(tod)]
        tod_of[t], oow[t] = tod, list(out_)
        oow_counts[t] = counts.loc[out_]
        rows.append(B.clock_table({t: (counts, tod)}))
        Lw = log2_normalised(counts.loc[win])
        LW[t] = Lw[[g for g in CLOCK if g in Lw.columns]]
        TODW[t] = tod.loc[win]
        print(t, f"window {len(win)} out-of-window {len(out_)}", f"{time.time() - t0:.0f} s", flush=True)
    table = pd.concat(rows, ignore_index=True)
    window_donors = sorted({d for t in tissues for d in TODW[t].index})
    tod_all = pd.concat(list(tod_of.values())).groupby(level=0).first()
    truth_donors = {d for t in tissues for d in oow[t]}
    assert not truth_donors & set(window_donors), "truth and test donors overlap"
    assert all(not set(oow[t]) & set(TODW[t].index) for t in tissues)

    # ---- truth (iii) -------------------------------------------------------------------------------------------------------
    t3, n3, dropped = truth_from(oow_counts, tod_of, oow, MIN_TRUTH)
    kept = list(t3.index)
    print(f"truth (iii): {len(kept)} tissues, dropped {len(dropped)}: {dropped}", flush=True)
    truth_i = pd.read_csv(exploratory_dir(AMENDMENT, False) / "truth_i.csv", index_col=0)  # stored A123 truth (i)

    # ---- the stored A123 clocksync estimate, regenerated in memory ----------------------------------------------------------
    wtab = table[table.donor.isin(window_donors) & in_window(table.tod)]
    wf = B.cross_fit_waveform(table, window_donors)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = B.clocksync_estimates(wtab, wf, 0, seed=SEED)
    check = ""
    om = exploratory_dir(AMENDMENT, False) / "organ_map.csv"
    if not quick and om.exists():
        omd = pd.read_csv(om, index_col=0)
        k = omd.index.intersection(est.index)
        ec = ev.centred(est, k)
        d_off = float(np.abs(ev.wrap(ec.offset_h - omd.est_offset_h.loc[k])).max())
        d_amp = float(np.abs(ec.log2_amp - omd.est_log2_amp.loc[k]).max())
        check = f"max |regenerated - stored organ_map.csv| over {len(k)} tissues: offset {d_off:.4f} h (stored has 4 decimals), log2 amplitude {d_amp:.4f}"
        print(check, flush=True)
        assert d_off < 1e-3 and d_amp < 1e-3, check

    tests = []
    main_row = score_row("clocksync_vs_truth_iii", "hypothesis", est, t3, n_perm, method="clocksync", truth="truth_iii")
    tests.append(dict(main_row, row="H123.5", hypothesis_stat="circ_corr_offset",
                      verdict="pass" if main_row["circ_corr_offset"] > 0.7 and main_row["p_offset"] < 0.05 else "fail"))
    tests.append(dict(main_row, row="H123.6", hypothesis_stat="pearson_log2_amp",
                      verdict="pass" if main_row["pearson_log2_amp"] > 0.7 else "fail"))

    # ---- descriptive: truth (i) vs truth (iii) ---------------------------------------------------------------------------------
    tests.append(score_row("truth_i_vs_truth_iii", "descriptive", truth_i, t3, n_perm, method="truth_i", truth="truth_iii"))

    # ---- descriptive: A123b methods, window design, cached orderings ----------------------------------------------------------------
    D = dict(L={"W": LW}, tod=tod_of)
    for m in METHODS:
        e, _, n_fail = M.assemble(m, "W", tissues, D, runs, tod_all)
        if e is None:
            tests.append(dict(row=f"{m}_vs_truth_iii", kind="descriptive", method=m, truth="truth_iii", n_failed_fits=n_fail, estimable=False))
            continue
        tests.append(score_row(f"{m}_vs_truth_iii", "descriptive", e, t3, n_perm, method=m, truth="truth_iii", n_failed_fits=n_fail,
                               n_tissues_estimated=len(e.dropna(subset=["amplitude", "offset_h"])), estimable=True))
    e = M.timed_estimates(LW, TODW)
    tests.append(score_row("timed_vs_truth_iii", "descriptive", e, t3, n_perm, method="timed", truth="truth_iii", estimable=True))

    # ---- descriptive: split sensitivity --------------------------------------------------------------------------------------------
    half = split_donors(oow, SEED)
    h1 = {d for d, v in half.items() if v == 1}
    frac = pd.Series({t: np.mean([half[d] == 1 for d in oow[t]]) for t in tissues if oow[t]})
    wf1 = B.cross_fit_waveform(table, set(table.donor) - h1)  # waveform from the half-1 donors only
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est_s = B.clocksync_estimates(wtab, wf1, 0, seed=SEED)
    split_info = {}
    for mn in (20, 10):
        d2 = {t: [d for d in oow[t] if half[d] == 2] for t in tissues}
        try:
            t3h, n2, drop2 = truth_from(oow_counts, tod_of, d2, mn)
        except RuntimeError as err:
            print("split truth", mn, err)
            continue
        assert not {d for t in tissues for d in d2[t]} & (set(window_donors) | h1)
        split_info[mn] = (len(t3h), drop2)
        tests.append(score_row(f"split_clocksync_vs_truth_iii_half2_min{mn}", "descriptive", est_s, t3h, n_perm, method="clocksync_half1_waveform",
                               truth=f"truth_iii_half2_min{mn}", estimable=True))
        tests.append(score_row(f"stored_clocksync_vs_truth_iii_half2_min{mn}", "descriptive", est, t3h, n_perm, method="clocksync",
                               truth=f"truth_iii_half2_min{mn}", estimable=True))
    print("split done", f"{time.time() - t0:.0f} s", flush=True)

    # ---- descriptive: prior-only baseline for the per-donor hour -------------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(wtab.drop(columns="tod"), anchors=sorted(wtab.compartment.unique()), waveform=wf)
        e_h, _ = T.host_estimate(ds, sorted(wtab.compartment.unique()))
    hr = e_h["donor_day_h"].dropna()  # tangent coordinate -> hour, 12 h flip resolved by the hospital prior
    tr = tod_all.reindex(hr.index).to_numpy(float) % 24
    p = hr.to_numpy(float)
    d_or = circ_median_offset(p, tr)
    c13 = float(T.wrap24(np.angle(np.exp(1j * W * p).mean()) / W))
    d_13 = float(T.wrap24(13.0 - c13))
    err_or = np.abs(T.cdiff(p + d_or, tr))
    err_13 = np.abs(T.cdiff(p + d_13, tr))
    err_base = np.abs(T.cdiff(13.0, tr))
    err_best = np.abs(T.cdiff(np.median(tr), tr))
    tests.append(dict(row="prior_baseline_donor_hour", kind="descriptive", method="clocksync_tangent_hour", truth="recorded_time",
                      n_donors=len(hr), mae_method_offset_min_h=float(np.median(err_or)), mae_method_offset_mean13_h=float(np.median(err_13)),
                      mae_always_13_h=float(np.median(err_base)), mae_best_constant_h=float(np.median(err_best)),
                      offset_min_h=d_or, offset_mean13_h=d_13))
    print("baseline done", f"{time.time() - t0:.0f} s", flush=True)

    # ---- outputs ----------------------------------------------------------------------------------------------------------------------
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "tests.csv", index=False)
    tc = ev.centred(t3, kept)
    td = pd.DataFrame(index=pd.Index(tissues, name="tissue"))
    td["n_window_donors"] = [int(cnt.get(t, np.nan)) for t in tissues]
    td["n_out_of_window_donors"] = [len(oow[t]) for t in tissues]
    td["kept"] = [t in kept for t in tissues]
    td["amplitude"] = t3.amplitude.reindex(td.index)
    td["offset_h"] = t3.offset_h.reindex(td.index)
    td["log2_amp_centred"] = tc.log2_amp.reindex(td.index)
    td["offset_h_centred"] = tc.offset_h.reindex(td.index)
    td.round(4).to_csv(out / "truth_iii.csv")

    pd.set_option("display.width", 250)
    show = ["row", "n_tissues", "circ_corr_offset", "circ_corr_offset_lo", "circ_corr_offset_hi", "p_offset", "pearson_log2_amp",
            "pearson_log2_amp_lo", "pearson_log2_amp_hi", "p_amp", "median_abs_offset_err_h"]
    print(tests[[c for c in show if c in tests]].round(3).to_string())
    print(tests[tests.row == "prior_baseline_donor_hour"].dropna(axis=1).T.to_string())

    nw = len(window_donors)
    readme = f"""# A123c: truth (iii) from out-of-window donors (generated by compute_truth_iii.py{' --subsample (CHECK MODE)' if quick else ''})

Registered in paper/docs/registered_experiments.md (A123c). Aggregates only; no donor rows.

## Counts
- A123 tissues: {len(tissues)}; truth (iii) tissues (>= {MIN_TRUTH} out-of-window fast deaths): {len(kept)}; dropped: {len(dropped)}
  ({', '.join(f'{t} ({n3[t]})' for t in dropped) or 'none'}).
- Window (test) donors: {nw}; out-of-window (truth) donors with a sample in an A123 tissue: {len(truth_donors)}; disjoint (asserted).
- Estimate regeneration check: {check or 'not run in check mode'}.
- Split: half 1 {len(h1)} donors, half 2 {len(truth_donors) - len(h1)}; per-tissue half-1 fraction min {frac.min():.2f}, max {frac.max():.2f}.
  Half-2 truth tissues: {', '.join(f'>= {k}: {v[0]} (dropped {len(v[1])})' for k, v in split_info.items())}.

## Choices beyond the registration
- Truth (iii) log2 normalisation (`log2_normalised`, DESeq2 size factors) is computed on the out-of-window samples of each tissue
  (truth (i): all fast deaths; estimates: window samples). One sample per donor and tissue (as `B.S.load`).
- The A123 estimate is regenerated with the A123 calls and n_boot = 0 (the donor bootstrap affects intervals, not points); checked
  against organ_map.csv. `compare` is the A123 function unchanged: 10,000 permutations, its own 1,000 tissue-bootstrap resamples.
- A123b methods: `compute_methods.assemble` on the cached W orderings; timed oracle = cosinor + SVD on the window donors' recorded times.
- Split: donors split once (seed 123) by `split_donors`; a donor is in one half for all its tissues; waveform from half-1 donors only
  (all other donors excluded), clocksync refitted on the window donors (n_boot = 0). Half-2 truth requires >= 20 half-2 donors per
  tissue (same rule), with a >= 10 variant. The stored A123 estimate is also scored against the half-2 truth.
- Prior baseline: the hour is `donor_day_h` (tangent coordinate with the 12 h flip resolved by the hospital prior, A123 H123.3's
  estimate). Common offset: (a) minimising the median absolute error (0.05 h grid, favours the method), (b) setting the circular
  mean of the hours to 13:00 (uses no recorded time). Baselines: always 13:00; best constant (median recorded time, uses truth).
"""
    (out / "README.md").write_text(readme)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
