"""A123b: every benchmark method on the 40-tissue GTEx working-hours organ map (registered in paper/docs/registered_experiments.md).

Designs: (W) the A123 design, fast deaths (Hardy 1-2) who died 08-18 h, tissues with >= 20 such donors, window samples only;
(F, descriptive) the same tissues with all fast deaths. Methods as in `results/figure_1/scripts/method/compute_gtex_benchmark.py`
(`run_order`, default settings): `cofe_genome`, `cofe_clock`, `cyclops`, `chiral` (each orders every tissue separately),
`timeteller` (supervised, 5-fold out-of-fold over all fast deaths; design W uses the window donors' predictions), `timed`
(oracle: per-tissue cosinor on the design's recorded times, rank-1 complex SVD) and `clocksync` (as A123; design W with the waveform
cross-fitted on the fast deaths outside the window, design F with the GTEx waveform; point estimates, no donor bootstrap).
Orderings -> organ map as A123 did for CHIRAL: every tissue aligned (rotation, reflection) to the reference tissue (Muscle -
Skeletal, else the tissue with most donors whose fit succeeded), averaged on the circle, per-tissue cosinor of the clock genes on
the consensus times, rank-1 complex SVD. A failed fit (error or all phases missing) drops the tissue; the failure rate is reported.
A method with estimates in < 20 tissues is not estimable (check mode: < 2).
Readouts: A123 `compare` (offset circular correlation, MAE, log2 amplitude Pearson; truth (i) and (ii)); for the orderings the donor
order (|circular correlation| of the consensus with recorded time after the best rotation and reflection, MAE, median per-tissue
value). H123.4 (design W, truth (i)): clocksync minus each of CoFE genome, CoFE clock genes, CYCLOPS in offset circular
correlation, paired tissue bootstrap (2,000 resamples, the same uniform draws for the three), one-sided p, Holm over the three.
Controlled access: per-donor phases live only in the cache `runs/` under the git-ignored output dir; only tissue-level aggregates
and summary statistics are written elsewhere.

    python compute_methods.py [--design W|F|both] [--methods ...] [--jobs J] [--threads T] [--subsample N] [--boot B] [--perm P]
`--jobs` x `--threads` = workers (<= 12 for now). Each (design, method, tissue) fit is cached in output/a123b/runs/ (delete a file
to refit); a rerun resumes. `--subsample N`: check mode, the N tissues with most window donors, 200 permutations, 200 bootstraps,
output under output/subsample/a123b/.
Outputs (exploratory/A123_working_hours/output/a123b/): tests.csv, donor_order.csv, h123_4.csv, runtime.csv.
"""

from __future__ import annotations

import concurrent.futures as cf
import importlib.util
import multiprocessing as mp
import shutil
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
_spec = importlib.util.spec_from_file_location("a123_wh", HERE.parent / "compute_working_hours.py")
A = importlib.util.module_from_spec(_spec)
sys.modules["a123_wh"] = A
_spec.loader.exec_module(A)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.data import in_window, log2_normalised  # noqa: E402
from common.gtex import LABELLED  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

B, CLOCK, REF, SEED = A.B, A.CLOCK, A.REF, A.SEED
ORDER = ("cofe_genome", "cofe_clock", "cyclops", "chiral")
ALL_METHODS = ORDER + ("timeteller", "timed", "clocksync")
COST = {"cyclops": 0, "cofe_genome": 1, "timeteller": 2, "chiral": 3, "cofe_clock": 4}  # launch order (longest first)
MIN_EST, N_H4, TT_SEED = 20, 2000, B.SEED
AMENDMENT = A.AMENDMENT


def slug(t: str) -> str:
    return t.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "-")


def cache_path(runs: Path, design: str, method: str, t: str) -> Path:
    return runs / (f"all_{method}_{slug(t)}.pkl" if method == "timeteller" else f"{design}_{method}_{slug(t)}.pkl")


# ---- one tissue fit (worker) --------------------------------------------------------------------------------------------

def fit_task(design: str, method: str, t: str, runs: str, threads: int) -> dict:
    """Fit `method` on tissue `t` of `design`, cache and return the result (phases per sample, runtime, failure)."""
    from clocksync.config import limit_threads

    limit_threads(threads)
    runs = Path(runs)
    f = cache_path(runs, design, method, t)
    if f.exists():
        return dict(key=(design, method, t), cached=True, **{k: v for k, v in pd.read_pickle(f).items() if k in ("runtime_s", "failed")})
    ts = time.time()
    counts, tod = B.S.load([t])[t]
    tod = tod.dropna()
    counts = counts.loc[tod.index]
    tmp = runs / "tmp" / f.stem
    err, phase, n_in = "", None, len(counts)
    try:
        if method == "timeteller":
            from common.benchmark.timeteller import cross_fit

            X = log2_normalised(counts)[[g for g in CLOCK if g in counts.columns]]
            phase = cross_fit(X, tod, tmp, seed=TT_SEED, cores=threads).time_pred
        else:
            c = counts.loc[tod.index[in_window(tod)]] if design == "W" else counts
            n_in = len(c)
            phase, meta = B.run_order(method, c, tmp, seed=0, ncores=threads)
            err = str(meta.get("error", ""))
    except Exception as e:  # a failed fit is data (registered failed-fit rule)
        err = f"{type(e).__name__}: {e}"[:300]
    shutil.rmtree(tmp, ignore_errors=True)
    phase = phase.dropna() if phase is not None else pd.Series(dtype=float)
    r = dict(phase=phase, runtime_s=time.time() - ts, failed=phase.empty, error=err, n_samples=len(phase),
             n_samples_in=n_in)
    pd.to_pickle(r, f)
    try:  # joblib's reusable pool can keep a pool worker from exiting
        from joblib.externals.loky import get_reusable_executor

        get_reusable_executor().shutdown(wait=True)
    except Exception:
        pass
    return dict(key=(design, method, t), cached=False, runtime_s=r["runtime_s"], failed=r["failed"])


# ---- data of the main process (clock genes only) ------------------------------------------------------------------------

def load_clock_data(tissues, quick):
    """Per tissue, in memory: recorded times of all fast deaths, log2-normalised clock genes of both designs' samples, the clock
    table (clocksync), truth (i) and (ii) coefficients."""
    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    o["tissue"] = o.tissue.astype(str)
    meta = pd.read_csv(A.META, sep=" ")
    meta["sample"] = meta.fullID.str.replace(".", "-", regex=False)
    ct = meta.set_index("sample").CT.dropna()
    D = dict(tod={}, L={"W": {}, "F": {}}, rows=[], Ci={}, Cii={}, n_all={}, n_win={})
    for t in tissues:
        counts, tod = B.S.load([t])[t]
        tod = tod.dropna()
        counts = counts.loc[tod.index]
        win = tod.index[in_window(tod)]
        D["tod"][t] = tod
        D["rows"].append(B.clock_table({t: (counts, tod)}))
        L = log2_normalised(counts)
        D["L"]["F"][t] = L[[g for g in CLOCK if g in L.columns]]
        Lw = log2_normalised(counts.loc[win])  # normalised on the window samples alone, as A123
        D["L"]["W"][t] = Lw[[g for g in CLOCK if g in Lw.columns]]
        if len(counts) >= 5:
            D["Ci"][t] = ev.gene_rhythms(L, tod, CLOCK)
        c2 = A.truth_ii_coef(a, o, ct, t)
        if c2 is not None:
            D["Cii"][t] = c2
        D["n_all"][t], D["n_win"][t] = len(counts), len(win)
        print("loaded", t, len(counts), len(win), flush=True)
    a.file.close()
    return D


def timed_estimates(Lt: dict, tod: dict) -> pd.DataFrame:
    C = {t: ev.gene_rhythms(L, tod[t].loc[L.index], CLOCK) for t, L in Lt.items() if len(L) >= 5}
    return ev.csvd_tissues(pd.DataFrame(C), REF)


def clocksync_est(design: str, table: pd.DataFrame, runs: Path) -> tuple[pd.DataFrame, float]:
    f = runs / f"{design}_clocksync_ALL.pkl"
    if f.exists():
        r = pd.read_pickle(f)
        return r["est"], r["runtime_s"]
    ts = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if design == "W":
            donors = sorted(table.donor[in_window(table.tod)].unique())
            wtab = table[table.donor.isin(donors) & in_window(table.tod)]
            est = B.clocksync_estimates(wtab, B.cross_fit_waveform(table, donors), 0, seed=SEED)
        else:
            est = B.clocksync_estimates(table, cs.gtex(), 0, seed=SEED)
    r = dict(est=est, runtime_s=time.time() - ts)
    pd.to_pickle(r, f)
    return est, r["runtime_s"]


def assemble(method, design, tissues, D, runs, tod_all):
    """Tissue estimates (or None), donor-order row (or None), failure counts, from the cached fits."""
    per, n_fail = {}, 0
    for t in tissues:
        r = pd.read_pickle(cache_path(runs, design, method, t))
        p = r["phase"]
        if design == "W" and method == "timeteller":
            p = p.loc[p.index.intersection(D["L"]["W"][t].index)]
        if r["failed"] or p.empty:
            n_fail += 1
        else:
            per[t] = p
    if not per:
        return None, None, n_fail
    ref = REF if REF in per else max(per, key=lambda t: len(per[t]))
    cons = B.consensus_times(per, ref)
    Cc = {}
    for t in per:
        ix = D["L"][design][t].index.intersection(cons.index)
        if len(ix) >= 5:
            Cc[t] = ev.gene_rhythms(D["L"][design][t].loc[ix], cons.loc[ix], CLOCK)
    est = ev.csvd_tissues(pd.DataFrame(Cc), ref) if ref in Cc else None
    trw = tod_all.reindex(cons.index)
    _, mae, _ = ev.align(cons, trw)
    pt = [abs(ev.circ_corr(p, D["tod"][t].reindex(p.index))) for t, p in per.items() if len(p) >= 5]
    don = dict(consensus_n=len(cons), consensus_circ_corr=abs(ev.circ_corr(cons, trw)), consensus_mae_h=mae,
               median_tissue_circ_corr=float(np.median(pt)) if pt else np.nan)
    return est, don, n_fail


# ---- H123.4 ---------------------------------------------------------------------------------------------------------------

def h123_4(ests: dict, truth: pd.DataFrame, n_boot: int, seed: int) -> pd.DataFrame:
    """clocksync minus each untimed method in offset circular correlation (design W, truth (i)), paired tissue bootstrap with the
    same uniform draws for every method, one-sided p = (1 + #{diff <= 0}) / (1 + B), Holm over the three (not estimable: p = 1)."""
    cs_est = ests["clocksync"]
    rng = np.random.default_rng(seed)
    U = rng.random((n_boot, len(truth)))
    rows = []
    for m in ("cofe_genome", "cofe_clock", "cyclops"):
        est = ests.get(m)
        k = (est.dropna(subset=["amplitude", "offset_h"]).index.intersection(truth.index).intersection(cs_est.index)
             if est is not None else pd.Index([]))
        if len(k) < MIN_EST:
            rows.append(dict(method=m, n_tissues=len(k), estimable=False))
            continue
        te_c, te_m = ev.tissue_errors(cs_est.loc[k], truth.loc[k]), ev.tissue_errors(est.loc[k], truth.loc[k])
        oc, om, ot = te_c.offset_h.to_numpy(), te_m.offset_h.to_numpy(), te_c.offset_h_true.to_numpy()
        r_c, r_m = ev.circ_corr(oc, ot), ev.circ_corr(om, ot)
        idx = np.floor(U[:, :len(k)] * len(k)).astype(int)
        d = np.array([ev.circ_corr(oc[i], ot[i]) - ev.circ_corr(om[i], ot[i]) for i in idx])
        rows.append(dict(method=m, n_tissues=len(k), estimable=True, circ_corr_clocksync=r_c, circ_corr_method=r_m,
                         diff=r_c - r_m, diff_lo=np.nanquantile(d, .025), diff_hi=np.nanquantile(d, .975),
                         p_one_sided=(1 + np.sum(d <= 0)) / (1 + n_boot), n_boot=n_boot))
    h = pd.DataFrame(rows)
    p = h.get("p_one_sided", pd.Series(dtype=float)).reindex(h.index).fillna(1.0).to_numpy()
    order = np.argsort(p)
    adj, run = np.empty(len(p)), 0.0
    for j, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - j) * p[i]))
        adj[i] = run
    h["p_holm"] = adj
    h["pass"] = np.where(h.estimable & (h.p_holm < 0.05), "pass", np.where(h.estimable, "fail", "not estimable"))
    return h


# ---- main -----------------------------------------------------------------------------------------------------------------

def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--design", choices=("W", "F", "both"), default="both")
    ap.add_argument("--methods", nargs="*", default=list(ALL_METHODS), choices=ALL_METHODS)
    ap.add_argument("--threads", type=int, default=4, help="threads per fit (workers = jobs x threads)")
    ap.add_argument("--boot", type=int, default=None, help="H123.4 resamples")
    ap.add_argument("--perm", type=int, default=None)
    ap.add_argument("--tissues", nargs="*", default=None, help="restrict to these tissues (timing checks)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_boot = args.boot if args.boot is not None else 200 if quick else N_H4
    n_perm = args.perm if args.perm is not None else 200 if quick else 10000
    global MIN_EST
    MIN_EST = 2 if quick else 20
    designs = ["W", "F"] if args.design == "both" else [args.design]
    out = exploratory_dir(AMENDMENT, quick) / "a123b"
    runs = out / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    print(f"workers = {args.jobs} x {args.threads} = {args.jobs * args.threads}; designs {designs}; methods {args.methods}", flush=True)
    t0 = time.time()

    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    o["tissue"] = o.tissue.astype(str)
    a.file.close()
    cnt = A.window_tissues(o)
    tissues = list(cnt.index[:args.subsample]) if quick else list(cnt.index)
    if quick and REF not in tissues:
        tissues = [REF] + tissues[:-1]
    if args.tissues:
        tissues = [t for t in cnt.index if t in args.tissues]
    print(f"{len(tissues)} tissues with >= {A.MIN_WINDOW} window donors", flush=True)
    D = load_clock_data(tissues, quick)
    tod_all = pd.concat(list(D["tod"].values())).groupby(level=0).first()
    truths = {"truth_i": ev.csvd_tissues(pd.DataFrame(D["Ci"]), REF),
              "truth_ii": ev.csvd_tissues(pd.DataFrame(D["Cii"]), REF) if REF in D["Cii"] else None}
    table = pd.concat(D["rows"], ignore_index=True)

    # ---- tissue fits in parallel (cache; longest first, design W first) ----------------------------------------------------
    todo = []
    for d in designs:
        for m in args.methods:
            if m in ORDER or m == "timeteller":
                for t in tissues:
                    k = (("F" if m == "timeteller" else d), m, t)
                    if not cache_path(runs, d, m, t).exists() and k not in todo:
                        todo.append(k)
    todo.sort(key=lambda k: (k[0] != "W", COST[k[1]], -D["n_all" if k[0] == "F" else "n_win"][k[2]]))
    print(f"{len(todo)} fits to run", flush=True)
    done, last = 0, time.time()
    if todo:
        with cf.ProcessPoolExecutor(args.jobs, mp_context=mp.get_context("spawn")) as ex:
            pend = {ex.submit(fit_task, d, m, t, str(runs), args.threads): (d, m, t) for d, m, t in todo}
            while pend:
                fin, _ = cf.wait(pend, timeout=600, return_when=cf.FIRST_COMPLETED)
                for fu in fin:
                    r = fu.result()
                    pend.pop(fu)
                    done += 1
                    print(f"[{time.time() - t0:.0f} s] {done}/{len(todo)} {r['key']} {r['runtime_s']:.0f} s failed={r['failed']}", flush=True)
                if not fin or time.time() - last > 600:
                    print(f"[{time.time() - t0:.0f} s] heartbeat: {done}/{len(todo)} done, {len(pend)} pending", flush=True)
                last = time.time()
            for p in list(getattr(ex, "_processes", {}).values()):  # idle pool workers of this run; all fits are cached
                p.terminate()

    # ---- estimates, scores -----------------------------------------------------------------------------------------------------
    tests, donor, runtime, ests = [], [], [], {}
    for d in designs:
        for m in args.methods:
            if m == "clocksync":
                est, rt = clocksync_est(d, table, runs)
                runtime.append(dict(design=d, method=m, tissue="ALL", runtime_s=round(rt, 1), failed=False))
                n_fail, don = 0, None
            elif m == "timed":
                est, don, n_fail = timed_estimates(D["L"][d], D["tod"]), None, 0
            else:
                est, don, n_fail = assemble(m, d, tissues, D, runs, tod_all)
                for t in tissues:
                    r = pd.read_pickle(cache_path(runs, d, m, t))
                    runtime.append(dict(design=d, method=m, tissue=t, n_samples=r["n_samples_in"], runtime_s=round(r["runtime_s"], 1),
                                        failed=r["failed"], error=r["error"]))
            ests[d, m] = est
            if don is not None:
                donor.append(dict(method=m, design=d, **don))
            nk = len(est.dropna(subset=["amplitude", "offset_h"])) if est is not None else 0
            for tn, tr in truths.items():
                base = dict(method=m, design=d, truth=tn, n_tissues_attempted=len(tissues), n_failed_fits=n_fail,
                            failure_rate=n_fail / len(tissues), n_tissues_estimated=nk)
                if tr is None or est is None or len(est.dropna(subset=["amplitude", "offset_h"]).index.intersection(tr.index)) < MIN_EST:
                    tests.append(dict(**base, estimable=False))
                    continue
                tests.append(dict(**base, estimable=True, **A.compare(est, tr, n_perm, SEED)))
    tests = pd.DataFrame(tests)
    tests["H123.1"] = np.where(tests.estimable & (tests.get("circ_corr_offset", 0) > 0.7) & (tests.get("p_offset", 1) < 0.05), "pass", "fail")
    tests["H123.2"] = np.where(tests.estimable & (tests.get("pearson_log2_amp", 0) > 0.7), "pass", "fail")
    tests.to_csv(out / "tests.csv", index=False)
    pd.DataFrame(donor).to_csv(out / "donor_order.csv", index=False)
    pd.DataFrame(runtime).to_csv(out / "runtime.csv", index=False)
    if "W" in designs and all(("W", m) in ests for m in ("clocksync", "cofe_genome", "cofe_clock", "cyclops")):
        h = h123_4({m: ests["W", m] for m in ("clocksync", "cofe_genome", "cofe_clock", "cyclops")}, truths["truth_i"], n_boot, SEED)
        h.to_csv(out / "h123_4.csv", index=False)
        print(h.round(4).to_string())
    pd.set_option("display.width", 250)
    cols = ["method", "design", "truth", "estimable", "failure_rate", "n_tissues", "circ_corr_offset", "median_abs_offset_err_h", "pearson_log2_amp"]
    print(tests[[c for c in cols if c in tests]].round(3).to_string())
    print(pd.DataFrame(donor).round(3).to_string())
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
