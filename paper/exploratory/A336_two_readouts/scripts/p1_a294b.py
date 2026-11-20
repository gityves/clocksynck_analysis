"""A336 part 1 wrapper: A294b H294b.1 / H294b.2 (own-variation floor of the 27 peripheral GTEx tissues) on the 14 genes,
the core loop (A329, existing), the trio contribution, the 364 gene drops and the trio alone.

Original: A294_floor/scripts/score_a294.py; core: A329_core_departure/{scripts/score_a329.py, output/tissues.csv, tests.csv}.
Fit as A294b (fast deaths, 27 peripheral tissues, every tissue an anchor, leave-one-out host; `A104.per_tissue`).
- 14 genes reproduced first: the clock fit is rerun and compared with A294 tissues.csv (own, S, coupling), and
  tests.csv (median log ratio, bootstrap interval, H294b.2 slope and interval, replaying score_a294's rng sequence).
- Core: A329 tissues.csv (log ratios against the 100 matched sets restricted to the core positions); the clock core fit of
  the drop loop must equal it (assert).
- Trio contribution: 200 tissue bootstraps (the statistic's unit; rng [336, 1]) of the paired per-tissue log ratios.
- Gene-drop null: clock fit only, no companions: the matched floor cannot be recomputed, so the drop statistics are the
  median own variance 2 psi of the clock (H294b.1) and the Theil-Sen slope of log own on log S (H294b.2, the descriptive
  slope of tests.csv), direction +1 (deviation noted in the output).
- Trio alone: the trio fit (own variance, no floor), if the model-free gate of the trio on the same rows passes.
Aggregates only (GTEx controlled access). Long run: `RUN_CAP_GB=30 paper/common/run.sh <this> --jobs 4`.

    python p1_a294b.py [--jobs 4] [--check]
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool

from p1_util import *  # noqa: F401,F403
import p1_util as U
import a336_common as C
from scipy.stats import theilslopes
sys.path.insert(0, str(U.PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads
from common import gtex

META = ["donor", "compartment", "depth"]
S = {}
A294 = U.EXP / "A294_floor" / "output"
A329 = U.EXP / "A329_core_departure" / "output"
SCRIPT = str(U.HERE.parent / "p1_a294b.py")


def load():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    return pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)


def fit(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    lab, w = S["sets"][i]
    pb = S["pb"]
    comps = sorted(pb.compartment.unique())
    try:
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=lab)
        pt = A104.per_tissue(ds)
        return pt.assign(label=lab, S=pt.z2 * pt.V, own=2 * pt.psi).rename_axis("tissue").reset_index()
    except Exception as e:  # noqa: BLE001
        return pd.DataFrame([dict(label=lab, error=f"{type(e).__name__}: {e}"[:200])])


def slope_stat(own, Sx, floor=None):
    """Theil-Sen slope of log excess (own - floor, positive tissues) or of log own on log S."""
    y = own - floor if floor is not None else own
    ok = y > 0
    return float(theilslopes(np.log(y[ok]), np.log(Sx[ok]))[0]) if ok.sum() > 2 else np.nan


def slope_boot(own, Sx, floor, rng, n=2000):
    """score_a294's H294.2 bootstrap on the positive-excess tissues."""
    ok = (own - floor) > 0
    xs, ys = np.log(Sx[ok]), np.log((own - floor)[ok])
    bs = []
    for _ in range(n):
        i = rng.integers(0, ok.sum(), ok.sum())
        bs.append(theilslopes(ys[i], xs[i])[0] if len(np.unique(xs[i])) > 2 else np.nan)
    return np.nanquantile(bs, [0.05, 0.95])


def signature(lo, hi):
    return ("phase/amplitude (slope 1)" if lo <= 1 <= hi and lo > 0 else
            "additive (slope 0)" if lo <= 0 <= hi and hi < 1 else "undecided")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    chk = a.check
    t0 = time.time()
    pb = load()
    limit_threads(1)
    comps = sorted(pb.compartment.unique())
    tw = C.trio()
    dl, tl = C.drops()
    if chk:
        dl = dl[:2] + [x for x in dl if x[0] == tl]
    sets = [("clock14", U.WF), ("trio", tw)] + [(lab, w) for lab, w in dl]
    keep = sorted(set(U.WF.genes))
    S.update(pb=pb[META + keep].copy(), sets=sets)
    del pb
    print(f"{len(sets)} fits {time.time() - t0:.0f}s", flush=True)

    # ---- 14 genes first ------------------------------------------------------------------------------------------
    r14 = fit(0)
    assert "error" not in r14, r14
    ref = pd.read_csv(A294 / "tissues.csv", index_col=0)
    c14 = r14.set_index("tissue")
    bad = [c for c in ("S", "own", "coupling") if not np.allclose(c14[c].reindex(ref.index), ref[c], atol=1e-6, rtol=0)]
    lr14 = np.log(c14.own.reindex(ref.index) / ref.own_matched_median)
    tests = pd.read_csv(A294 / "tests.csv").set_index("test")
    rng = np.random.default_rng(294)
    lo14, hi14 = np.quantile([np.median(rng.choice(lr14.to_numpy(), 27)) for _ in range(2000)], [0.05, 0.95])
    ex = ref.own - ref.own_matched_median
    s14 = slope_stat(c14.own.reindex(ref.index).to_numpy(), c14.S.reindex(ref.index).to_numpy(), ref.own_matched_median.to_numpy())
    blo14, bhi14 = slope_boot(c14.own.reindex(ref.index).to_numpy(), c14.S.reindex(ref.index).to_numpy(), ref.own_matched_median.to_numpy(), rng)
    h1, h2 = tests.loc["H294.1"], tests.loc["H294.2"]
    rep = {"median_log_ratio": (float(np.median(lr14)), h1.median_log_ratio), "lo90": (lo14, h1.lo90), "hi90": (hi14, h1.hi90),
           "slope": (s14, h2.slope), "slope_lo90": (blo14, h2.lo90), "slope_hi90": (bhi14, h2.hi90)}
    badr = {k: v for k, v in rep.items() if not np.isclose(*v, atol=1e-6, rtol=0)}
    if (bad or badr) and not chk:
        print("14-gene NOT reproduced", bad, badr, flush=True)
        rows = [dict(id=i, figure_panel="Fig 4e", value_14_stored=v, value_14=w_, reproduced=False,
                     notes=f"14-gene values not reproduced (1e-6): columns {bad}, stats {badr}", script=SCRIPT)
                for i, v, w_ in (("A294b_H294b.1", h1.median_log_ratio, rep["median_log_ratio"][0]),
                                 ("A294b_H294b.2", h2.slope, rep["slope"][0]))]
        C.write_rows("a294b", rows)
        return
    print("14-gene reproduction OK", flush=True)

    # ---- fits of the trio and the drops --------------------------------------------------------------------------
    res = []
    with Pool(a.jobs) as pool:
        for k, r in enumerate(pool.imap(fit, range(1, len(sets)), chunksize=1)):
            res.append(r)
            if k % 10 == 0:
                print(f"{k + 1}/{len(sets) - 1} fits {time.time() - t0:.0f}s", flush=True)
    res = pd.concat(res, ignore_index=True)
    err = res.get("error", pd.Series(dtype=str)).notna().sum()
    print(f"fits done {time.time() - t0:.0f}s, errors {err}", flush=True)
    tr = res[res.label == "trio"]
    dr = res[res.label.isin([lab for lab, _ in dl])]

    # ---- core (A329, existing) and its reproduction by the trio-drop fit -----------------------------------------
    t29 = pd.read_csv(A329 / "tissues.csv", index_col=0)
    cdrop = dr[dr.label == tl].set_index("tissue").reindex(t29.index)
    assert np.allclose(cdrop.own, t29.own_core, atol=1e-6, rtol=0) and np.allclose(cdrop.S, t29.S_core, atol=1e-6, rtol=0)
    lr_f, lr_c = t29.log_ratio_full.to_numpy(), t29.log_ratio_core.to_numpy()
    rngc = np.random.default_rng(329)
    lo_c, hi_c = np.quantile([np.median(rngc.choice(lr_c, 27)) for _ in range(2000)], [0.05, 0.95])
    own_c, S_c, fl_c = t29.own_core.to_numpy(), t29.S_core.to_numpy(), t29.own_matched_median_core.to_numpy()
    s_c = slope_stat(own_c, S_c, fl_c)
    blo_c, bhi_c = slope_boot(own_c, S_c, fl_c, np.random.default_rng([336, 2]))

    # ---- trio contribution (200 tissue bootstraps of the paired per-tissue quantities) ---------------------------
    rng = np.random.default_rng(U.SEED_BOOT)
    own_f, S_f, fl_f = t29.own_full.to_numpy(), t29.S_full.to_numpy(), t29.own_matched_median_full.to_numpy()
    d1, d2 = [], []
    for _ in range(3 if chk else U.N_BOOT_UNIT):
        i = rng.integers(0, 27, 27)
        d1.append(np.median(lr_f[i]) - np.median(lr_c[i]))
        d2.append(slope_stat(own_f[i], S_f[i], fl_f[i]) - slope_stat(own_c[i], S_c[i], fl_c[i]))

    # ---- gene drops (clock only) ---------------------------------------------------------------------------------
    wfs = dict(dl)
    v1, v2 = {}, {}
    for lab, g in dr.groupby("label"):
        v1[lab] = float(g.own.median()) if len(g) == 27 else np.nan
        v2[lab] = slope_stat(g.own.to_numpy(), g.S.to_numpy()) if len(g) == 27 else np.nan
    # ---- trio alone ----------------------------------------------------------------------------------------------
    pbt = None
    gp = None
    try:
        pbt = load()
        ds = cs.Dataset(pbt[META + list(tw.genes)], anchors=comps, waveform=tw, name="trio")
        g = ds.gate(n_perm=1000 if not chk else 50, seed=336)
        gp = (bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), float(g.p_pairing), float(g.p_acrophase))
    except Exception as e:  # noqa: BLE001
        gp = ("not computable", np.nan, np.nan)
        print("gate error", e, flush=True)
    trio_ok = len(tr) == 27 and "error" not in tr
    own_t = float(tr.own.median()) if trio_ok else np.nan
    slope_t = slope_stat(tr.own.to_numpy(), tr.S.to_numpy()) if trio_ok else np.nan

    if chk:
        dp = C.drop_percentiles(v1, tl, wfs, +1)
        assert len(d1) == 3 and not np.isnan(lo_c)
        print("check ok")
        return
    out = []
    sig14, sig_c = signature(blo14, bhi14), signature(blo_c, bhi_c)
    spec = [("A294b_H294b.1", "median over 27 tissues of log(2psi_clock / median 2psi of 100 matched sets)", float(np.median(lr14)),
             h1.median_log_ratio, float(np.median(lr_c)), (lo_c > 0 and np.median(lr_c) > 0), "median log ratio > 0 AND tissue-bootstrap 90% lower bound > 0",
             d1, v1, own_t, f"core median {np.median(lr_c):.4f} (90% {lo_c:.4f} to {hi_c:.4f}) vs 14 genes {np.median(lr14):.4f} ({lo14:.4f} to {hi14:.4f}); "
             "gene-drop statistic = median own variance 2psi of the clock (no floor without companions), direction +1"),
            ("A294b_H294b.2", "Theil-Sen slope of log excess own variance on log S_clock (27 tissues)", s14, h2.slope, s_c,
             sig_c == sig14, f"core signature (slope interval vs 0 and 1, as score_a294) equals the 14-gene one ({sig14})", d2, v2, slope_t,
             f"core slope {s_c:.4f} (90% {blo_c:.4f} to {bhi_c:.4f}, core interval from rng [336,2]) signature {sig_c}; 14 genes {s14:.4f} ({blo14:.4f} to {bhi14:.4f}) {sig14}; "
             "gene-drop statistic = slope of log own on log S (no floor), direction +1")]
    for i, name, v14, vstored, vcore, ok, crit, dd, vals, vtrio, note in spec:
        dp = C.drop_percentiles(vals, tl, wfs, +1)
        st = "holds" if ok else "falls"
        lo, hi = U.ci90(dd)
        out.append(dict(id=i, figure_panel="Fig 4e", statistic=name, direction=1, value_14_stored=vstored, value_14=v14, reproduced=True,
                        value_core=vcore, core_criterion=crit, core_status=st, diff_14_minus_core=v14 - vcore, diff_lo90=lo, diff_hi90=hi,
                        diff_draws="200 tissue bootstraps of paired per-tissue values, seed [336,1]", trio_gate_passes=gp[0],
                        trio_alone_value=vtrio if gp[0] is True else np.nan, reading=C.reading(st, dp["trio_drop_below_q10"]),
                        blind="core not blind (A335 / A329)", script=SCRIPT,
                        notes=note + f"; gate p_pairing {gp[1]:.4f}, p_acrophase {gp[2]:.4f}; fit errors {err}",
                        **{x: dp[x] for x in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")}))
    df = C.write_rows("a294b", out)
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
