"""Figure 1, GTEx benchmark: hidden times of death, `clocksync` against cyclic-ordering methods and the timed oracle.

Analysis (post-hoc design of the figure, not an amendment; the legacy `fig1_method` runs): GTEx fast deaths (Hardy 1-2) with the time
of death hidden. A draw is `n` donors from all fast deaths (`full`) or from those who died 08-18 h (`window`, hospital hours); each
of eight tissues contributes its samples of those donors. Methods: COFE genome-wide (`cofe_genome`) and on the clock genes
(`cofe_clock`), CYCLOPS 2.0 (`cyclops`), CHIRAL (`chiral`), TimeTeller trained on recorded times with 5-fold cross-fitting
(`timeteller`, a supervised reference), `clocksync` (untimed host model, waveform refitted without the drawn donors,
`Dataset.fit` with all tissues as anchors, donor bootstrap), and `timed`, the oracle (cosinor of the clock genes on the drawn
donors' recorded times, rank-1 complex SVD). Ordering methods are scored on their donor times after the best rotation and
reflection (oracle orientation, favourable to them); every method also on the tissue estimands (log2 amplitude and offset, reference
free) against the same SVD on the recorded times of all fast deaths (`truth`), and clocksync against its own timed estimand
(`truth_model`: complex loading of the clock-plane state on recorded time, global waveform) for the coverage of its 95 % intervals.
Readouts per run: tissue errors (log2 amplitude, offset, h), donor-time errors and skill over the best constant (orderings),
coverage and width of the clocksync intervals, runtime, CCD (clock correlation distance) of every tissue draw.
Only aggregates are written (no per-donor orderings or times: controlled access).

    python compute_gtex_benchmark.py [--subsample N] [--methods ...] [--ns 80 40 ...] [--draws K] [--tissues ...]
`--tissue-set working_hours`: the 40 tissues with >= 20 fast deaths who died 08-18 h (A123), outputs in `gtex_benchmark_40/`
(presentation R1.2, post hoc). `cyclops_clock`: CYCLOPS 2.0 on the 14 core clock genes only (expression and seed genes), 3 draws per n (cost, ~4 min per tissue).
`--subsample N`: check mode, one draw of N donors per design, five tissues, cheap methods (clocksync, timed, chiral, cofe_clock).
Full run (legacy design): n in (10, 20, 40, 80), 10 draws (cofe_genome and cyclops 3 at n = 80), methods per `DRAWS_80` / `DRAWS_SWEEP`.
Outputs: paper/results/figure_1/output/method/gtex_benchmark/{tissues,samples,genes,coverage,runtime,ccd,truth,truth_model,gene_truth}.csv
(runs are cached in `runs/`; delete a file to rerun it). Legacy: legacy/paper_code/fig1_method/{run_gtex_benchmark.py,
2_gtex_benchmark.ipynb, 6_gtex_sweep.ipynb} (outputs/paper/fig1_method/gtex/).
"""

from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402
from common import cli  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.cyclops import seeds as cyclops_seeds  # noqa: E402
from common.benchmark.data import cofe_filter, in_window, log2_normalised  # noqa: E402
from common.benchmark.tcga import COFE_CLOCK  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402
from _shared import donor_times, load  # noqa: E402

FOLDER = "method"
REF = "Muscle - Skeletal"
TISSUES = ("Muscle - Skeletal", "Thyroid", "Artery - Tibial", "Nerve - Tibial", "Skin - Sun Exposed (Lower leg)",
           "Skin - Not Sun Exposed (Suprapubic)", "Adipose - Subcutaneous", "Lung")
NS = (10, 20, 40, 80)
DRAWS_80 = {"cofe_genome": 3, "cofe_clock": 10, "cyclops": 3, "cyclops_clock": 3, "chiral": 10, "timeteller": 10,
            "clocksync": 10, "timed": 10}
DRAWS_SWEEP = {"cofe_clock": 10, "cyclops_clock": 3, "chiral": 10, "timeteller": 10, "clocksync": 10, "timed": 10}
ORDERINGS = ("chiral", "timeteller", "cofe_clock", "cofe_genome", "cyclops", "cyclops_clock")
MIN_WINDOW = 20  # `--tissue-set working_hours`: tissues with >= 20 fast deaths who died 08-18 h (as A123)
DESIGNS = ("full", "window")
SEED = 2026
AMP_FAIL, OFFSET_FAIL = 0.5, 2.0
CLOCK = S.CLOCK
LN2 = np.log(2.0)


# ---- data -------------------------------------------------------------------------------------------------------------

def working_hours_tissues() -> tuple[str, ...]:
    """Tissues with >= MIN_WINDOW fast deaths who died 08-18 h (cultured cells excluded), most donors first: A123's 40 tissues."""
    import anndata as ad
    from common.gtex import FAST, LABELLED

    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    a.file.close()
    o["tissue"] = o.tissue.astype(str)
    o = o[~o.tissue.str.startswith("Cells -") & o.DTHHRDY.isin(list(FAST)) & o.tod.notna()].drop_duplicates(["subject", "tissue"])
    n = o[in_window(o.tod.astype(float) % 24)].groupby("tissue").subject.nunique()
    return tuple(n[n >= MIN_WINDOW].sort_values(ascending=False).index)


def draw(data, design: str, n: int, rng) -> np.ndarray:
    tod = donor_times(data)
    pool = tod.index[in_window(tod)] if design == "window" else tod.index
    return np.sort(rng.choice(pool, min(n, len(pool)), replace=False))


def draws(data, n: int, k: int) -> dict:
    rng = np.random.default_rng(SEED if n == 80 else [SEED, n])
    return {(d, i): draw(data, d, n, rng) for d in DESIGNS for i in range(k)}


def clock_table(data, donors=None) -> pd.DataFrame:
    """Long table (donor, compartment, depth, tod, 14 clock-gene counts)."""
    rows = []
    for t, (counts, tod) in data.items():
        c = counts if donors is None else counts.loc[counts.index.intersection(donors)]
        df = c.reindex(columns=CLOCK, fill_value=0).astype(float)
        df.insert(0, "donor", c.index)
        df.insert(1, "compartment", t)
        df.insert(2, "depth", c.sum(axis=1).to_numpy())
        df.insert(3, "tod", tod.loc[c.index].to_numpy())
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


# ---- estimates --------------------------------------------------------------------------------------------------------

def csvd_estimates(data, donors, times: pd.Series) -> pd.DataFrame:
    """Tissue amplitude and offset relative to the reference from donor times: per-tissue cosinor of the clock genes (log2
    normalised), rank-1 complex SVD (`evaluate.csvd_tissues`). local: not in clocksync."""
    t = times.dropna()
    C = {}
    for k, (counts, _) in data.items():
        x = counts.loc[counts.index.intersection(donors).intersection(t.index)]
        if len(x) >= 5:
            C[k] = ev.gene_rhythms(log2_normalised(x), t.loc[x.index], CLOCK)
    return ev.csvd_tissues(pd.DataFrame(C), REF)


def model_truth(table: pd.DataFrame, wf, times: pd.Series) -> pd.DataFrame:
    """clocksync's own estimand when the times are known: per-tissue complex loading of the clock-plane state on the recorded time
    (`Dataset.data` states, local `timed_loading`), amplitude and offset relative to the reference tissue."""
    ds = cs.Dataset(table.drop(columns="tod"), anchors=sorted(table.compartment.unique()), waveform=wf)
    d = ds.data
    z = {}
    for k in ds.compartments:
        i = d.compartment == k
        z[k] = S.timed_loading(d.beta[i], times.reindex(d.resection[i]).to_numpy())
    z = pd.Series(z)
    r = z / z[REF]
    return pd.DataFrame({"amplitude": np.abs(r), "offset_h": np.angle(r) * S.HOURS}).rename_axis("compartment")


def clocksync_estimates(table, wf, n_boot: int, seed: int = 0) -> pd.DataFrame:
    """Untimed host model: `Dataset.fit` with every tissue an anchor (free host), relative amplitude and offset with donor-bootstrap
    percentile intervals; tissues with < 5 donors dropped."""
    cnt = table.compartment.value_counts()
    sub = table[table.compartment.isin(cnt[cnt >= 5].index)].drop(columns="tod")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(sub, anchors=sorted(sub.compartment.unique()), waveform=wf)
        f = ds.fit(n_boot=n_boot, seed=seed, gate=False).table
    est = pd.DataFrame({"amplitude": f.rel_amplitude, "offset_h": f.offset_h, "coupling": f.coupling,
                        "log2_amp_lo": f.log_amp_lo / LN2, "log2_amp_hi": f.log_amp_hi / LN2,
                        "offset_lo": f.offset_lo, "offset_hi": f.offset_hi})
    est["n_boot"] = f.n_boot
    return est.rename_axis("compartment")


def cross_fit_waveform(table, donors):
    """GTEx clock waveform refitted on the fast deaths other than `donors`, tissue as group (`waveform.fit_waveform`), so that a
    draw's estimates do not use a waveform learnt from its own samples."""
    t = table[~table.donor.isin(donors)]
    y = np.log(np.maximum(t[CLOCK].to_numpy(float), 0.5) / t.depth.to_numpy()[:, None])
    return fit_waveform(y, t.tod.to_numpy(), CLOCK, t.compartment.to_numpy(), name="cross-fitted")


# ---- orderings --------------------------------------------------------------------------------------------------------

def run_order(method: str, counts: pd.DataFrame, tmp: Path, seed: int = 0, ncores: int = 4):
    """Predicted phase (h, arbitrary origin and direction) of each sample of one tissue and run metadata."""
    if method in ("cofe_genome", "cofe_clock"):
        from common.benchmark.cofe import run_cofe

        X = cofe_filter(counts)
        if method == "cofe_clock":
            X = X[[g for g in CLOCK if g in X.columns]]
        X = X.loc[:, X.std() > 0]
        try:
            r = run_cofe(X, seed=seed, ncores=ncores)
        except ValueError as e:  # degenerate cross-validation fold at small n: a failed run
            return pd.Series(np.nan, index=counts.index), {"error": str(e), "n_genes": X.shape[1]}
        return r["phase"].rename(None), {"best_s": r["best_s"], "interior_min": r["interior_min"], "n_genes": X.shape[1]}
    if method == "chiral":
        from common.benchmark.chiral import run_chiral

        phase, failed = run_chiral(log2_normalised(counts), seed=seed)
        return phase.rename(None), {"failed_starts": failed}
    if method in ("cyclops", "cyclops_clock"):
        from common.benchmark.cyclops import run_cyclops

        expr = 2 ** log2_normalised(counts) - 1
        seeds = cyclops_seeds()
        if method == "cyclops_clock":  # core clock genes only, as expression and as seed genes
            seeds = [g for g in CLOCK if g in expr.columns]
            expr = expr[seeds]
        r = run_cyclops(expr.T, seeds, tmp, seed=seed, threads=ncores)
        return (r["phase_raw"] * 24 / (2 * np.pi)) % 24, {}
    raise ValueError(method)


def consensus_times(pred: dict[str, pd.Series], ref: str = REF) -> pd.Series:
    """Donor times from per-tissue orderings: each aligned (rotation, reflection) to the reference tissue's ordering on shared
    donors, then averaged on the circle."""
    base = pred[ref]
    aligned = [base]
    for t, p in pred.items():
        if t == ref:
            continue
        shared = p.index.intersection(base.index)
        if len(shared) < 5:
            continue
        a, _, d = ev.align(p.loc[shared], base.loc[shared])
        off = ev.wrap(a - d * p.loc[shared].to_numpy())[0]
        aligned.append((d * p + off) % 24)
    z = pd.concat(aligned, axis=1).apply(lambda r: np.exp(1j * ev.W * r.dropna()).mean(), axis=1)
    return pd.Series(np.angle(z.to_numpy(complex)) / ev.W % 24, index=z.index)


def timeteller_oof(data, tmp: Path, seed: int, cores: int) -> dict:
    """Out-of-fold TimeTeller predictions of every fast death, per tissue (in memory)."""
    from common.benchmark.timeteller import cross_fit

    out = {}
    for t, (counts, tod) in data.items():
        X = log2_normalised(counts)[[g for g in CLOCK if g in counts.columns]]
        out[t] = cross_fit(X, tod, tmp / "timeteller" / t.replace(" ", "_"), seed=seed, cores=cores)
        print("timeteller", t, flush=True)
    return out


# ---- main -------------------------------------------------------------------------------------------------------------

def score_run(method, n, design, i, key, est, truth, truth_model, data, donors, tod, gene_truth, per, meta, rows):
    """Append the scores of one run to `rows` (dict of lists)."""
    errs = ev.tissue_errors(est, truth).reindex(list(data)).assign(**key)
    rows["tissues"].append(errs)
    if method == "clocksync":
        for name, tr in (("csvd", truth), ("model", truth_model)):
            e = ev.tissue_errors(est, tr)
            lo, hi = est.loc[e.index, ["log2_amp_lo", "log2_amp_hi"]].T.to_numpy()
            olo, ohi = est.loc[e.index, ["offset_lo", "offset_hi"]].T.to_numpy()
            d = -e.offset_err_h
            rows["coverage"].append(pd.DataFrame({
                "amp_covered": (e.log2_amp_true >= lo) & (e.log2_amp_true <= hi),
                "offset_covered": (d >= olo - e.offset_h) & (d <= ohi - e.offset_h),
                "amp_width": hi - lo, "offset_width": ohi - olo}, index=e.index).assign(truth=name, **key))
    if per is None:
        return
    for t, p in per.items():
        true = tod.loc[p.index]
        aligned, mae, _ = ev.align(p, true)
        counts = data[t][0].loc[p.index]
        arr = ev.arrangement_error(ev.cosinor(cofe_filter(counts), aligned, COFE_CLOCK), gene_truth.loc[t])
        rows["genes"].append(arr.assign(tissue=t, **key))
        mae_abs = float(np.median(np.abs(ev.wrap(p - true)))) if method == "timeteller" else np.nan
        rows["samples"].append(dict(tissue=t, n_samples=len(p), mae_h=mae, mae_abs_h=mae_abs, skill=ev.skill(mae, true),
                                    circ_corr=abs(ev.circ_corr(p, true)), R_true=ev.spread(true), R_pred=ev.spread(p),
                                    interior_min=meta.get(t, {}).get("interior_min"), **key))
    cons = consensus_times(per)
    _, mae, _ = ev.align(cons, tod.loc[cons.index])
    rows["samples"].append(dict(tissue="consensus", n_samples=len(cons), mae_h=mae, skill=ev.skill(mae, tod.loc[cons.index]),
                                circ_corr=abs(ev.circ_corr(cons, tod.loc[cons.index])), R_true=ev.spread(tod.loc[cons.index]),
                                R_pred=ev.spread(cons), **key))


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--methods", nargs="*", default=None)
    ap.add_argument("--ns", nargs="*", type=int, default=None)
    ap.add_argument("--draws", type=int, default=None, help="draws per design (default: legacy design)")
    ap.add_argument("--tissues", nargs="*", default=None)
    ap.add_argument("--ncores", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--tissue-set", choices=("benchmark", "working_hours"), default="benchmark",
                    help="benchmark: the eight tissues of TISSUES; working_hours: A123's 40 tissues (output folder gtex_benchmark_40)")
    ap.add_argument("--no-ccd", action="store_true", help="skip the clock correlation distance of every tissue draw")
    args = ap.parse_args()
    opt = cli.setup(args)
    quick = opt["quick"]
    out = out_dir(FOLDER, quick) / ("gtex_benchmark_40" if args.tissue_set == "working_hours" else "gtex_benchmark")
    runs = out / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    tissues = (tuple(args.tissues) if args.tissues else working_hours_tissues() if args.tissue_set == "working_hours"
               else TISSUES[:5] if quick else TISSUES)
    ns = tuple(args.ns) if args.ns else (args.subsample,) if quick else NS
    methods = args.methods or (["clocksync", "timed", "chiral", "cofe_clock"] if quick else list(DRAWS_80))
    n_boot = args.n_boot if args.n_boot is not None else 10 if quick else 200
    ccd_perm = 20 if quick else 1000

    from common.benchmark.ccd import ccd

    t0 = time.time()
    data = load(tissues)
    tod = donor_times(data)
    print({t: len(c) for t, (c, _) in data.items()}, f"{time.time() - t0:.0f} s", flush=True)
    table = clock_table(data)
    wf_all = cs.gtex()
    truth = csvd_estimates(data, tod.index, tod)
    truth_model = model_truth(table, wf_all, tod)
    truth.to_csv(out / "truth.csv")
    truth_model.to_csv(out / "truth_model.csv")
    gene_truth = pd.concat({t: ev.cosinor(cofe_filter(c), tt, COFE_CLOCK) for t, (c, tt) in data.items()}, names=["tissue", "gene"])
    gene_truth.to_csv(out / "gene_truth.csv")

    def n_draws(m, n):
        k = (DRAWS_80 if n == 80 else DRAWS_SWEEP).get(m, 0)
        return min(k, args.draws) if args.draws else (1 if quick else k)

    rows = {k: [] for k in ("tissues", "samples", "genes", "coverage", "runtime", "ccd")}
    oof = None
    for n in ns:
        D = draws(data, n, max((DRAWS_80 if n == 80 else DRAWS_SWEEP).values()))  # the generator advances over all legacy draws
        for m in methods:
            if m == "timeteller" and oof is None and n_draws(m, n):
                oof = timeteller_oof(data, runs, SEED, args.ncores)
            for design in DESIGNS:
                for i in range(n_draws(m, n)):
                    donors = D[design, i]
                    key = dict(method=m, n=n, design=design, draw=i)
                    f = runs / f"{m}_n{n}_{design}_{i}.pkl"
                    if f.exists():
                        r = pd.read_pickle(f)
                    else:
                        r = {k: [] for k in rows}
                        ts = time.time()
                        per, meta = None, {}
                        if m == "timed":
                            est = csvd_estimates(data, donors, tod)
                        elif m == "clocksync":
                            tab = clock_table(data, donors)
                            wf = cross_fit_waveform(table, donors)
                            est = clocksync_estimates(tab, wf, 0)
                            t_fit = time.time() - ts
                            est = clocksync_estimates(tab, wf, n_boot, seed=i) if n_boot else est
                            r["runtime"].append(dict(runtime_s=round(t_fit, 2), bootstrap_s=round(time.time() - ts - t_fit), **key))
                        else:
                            per = {}
                            for t, (counts, _) in data.items():
                                c = counts.loc[counts.index.intersection(donors)]
                                if len(c) < 5:
                                    continue
                                if m == "timeteller":
                                    o = oof[t].loc[c.index]
                                    ph, mt = o.time_pred, {"theta": float(o.theta.median())}
                                else:
                                    ph, mt = run_order(m, c, runs / "tmp" / f"{m}_n{n}_{design}_{i}_{t.replace(' ', '_')}", seed=i, ncores=args.ncores)
                                per[t] = ph.dropna()
                                meta[t] = mt
                            r["runtime"].append(dict(runtime_s=round(time.time() - ts), **key))
                            if REF not in per or per[REF].empty:
                                est = pd.DataFrame(columns=["amplitude", "offset_h"])
                            else:
                                est = csvd_estimates(data, donors, consensus_times(per))
                        score_run(m, n, design, i, key, est, truth, truth_model, data, donors, tod, gene_truth, per, meta, r)
                        pd.to_pickle(r, f)
                        print(m, f"n{n}", design, i, f"{time.time() - ts:.0f} s", flush=True)
                    for k in rows:
                        rows[k] += r[k]
        for design in () if args.no_ccd else DESIGNS:  # clock correlation distance of every tissue draw
            for (d, i), donors in D.items():
                if d != design:
                    continue
                for t, (counts, _) in data.items():
                    c = counts.loc[counts.index.intersection(donors)]
                    if len(c) >= 5:
                        v, p = ccd(log2_normalised(c), n_perm=ccd_perm, seed=i)
                        rows["ccd"].append(dict(tissue=t, design=design, draw=i, n_donors=len(c), n=n, ccd=v, p=p))

    for name, fr in rows.items():
        if not fr:
            continue
        df = pd.concat(fr) if isinstance(fr[0], pd.DataFrame) else pd.DataFrame(fr)
        if name in ("tissues", "coverage", "genes"):
            df = df.rename_axis("tissue" if name != "genes" else "gene").reset_index()
        df.to_csv(out / f"{name}.csv", index=False)
    T = pd.concat(rows["tissues"]).rename_axis("tissue").reset_index()
    T = T[T.tissue != REF]
    T["failed"] = T.log2_amp_err.isna() | (T.log2_amp_err.abs() > AMP_FAIL) | (T.offset_err_h.abs() > OFFSET_FAIL)

    def summ(g):
        ok = g.dropna(subset=["log2_amp_err"])
        return pd.Series({"amp_err": ok.log2_amp_err.abs().median(), "offset_err_h": ok.offset_err_h.abs().median(),
                          "r_amp": np.corrcoef(ok.log2_amp, ok.log2_amp_true)[0, 1] if len(ok) > 2 else np.nan,
                          "r_offset": np.corrcoef(ok.offset_h_true + ok.offset_err_h, ok.offset_h_true)[0, 1] if len(ok) > 2 else np.nan,
                          "failure_rate": g.failed.mean(), "n_estimates": len(ok)})

    S_ = T.groupby(["method", "design", "n"]).apply(summ).reset_index()
    S_.to_csv(out / "summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(S_.round(2).to_string(index=False))
    if rows["coverage"]:
        C = pd.concat(rows["coverage"])
        print(C.groupby(["truth", "design", "n"])[["amp_covered", "offset_covered", "amp_width", "offset_width"]]
              .agg({"amp_covered": "mean", "offset_covered": "mean", "amp_width": "median", "offset_width": "median"}).round(2))


if __name__ == "__main__":
    main()
