"""Baboon degraded to working-hours windows (A40): are tissue offsets, order and window centre recoverable from a few
consecutive animals?

Pre-registered as A40 (2026-10-16; docs/PLAN_internal_desync.md (archive)). Data: 12 baboons killed every 2 h (ZT0-ZT22), 63
tissues as compartments and as anchors (`data/baboon`), human GTEx waveform. Windows: cyclic runs of 3, 4, 5 (8 h,
primary) and 6 animals from each of the 12 starts; the A36 estimator refitted on the window's animals only; "hospital-like"
window ZT2-ZT10 (5 animals) by name. Readouts per window: (1) untimed tissue offsets (arg z relative to the circular mean)
against the full-day timed offsets (-arg of the tissue's gain against the GTEx waveform, cosinor on ZT): Fisher-Lee and
Jammalamadaka-SenGupta circular correlation and median absolute error after a common shift; (2) ordering: Spearman of the
tangent coordinate of the A25 estimator (direction by curvature, no prior) with ZT, and whether the flip is correct (rho > 0);
(3) window centre: curvature centre against the true centre, in the full-day species frame F (the circular mean of host
phase minus ZT of the full-day fit) and without it. Primary (8 h windows): (1) median Fisher-Lee r > 0.5 over the 12
windows; (2) median |Spearman| > 0.7 and flip correct in >= 9/12.

    python compute_baboon_windows.py [--subsample N] [--seed S] [--jobs 4] [--n-perm N]
`--subsample N` runs the first N start positions of each window size.
Outputs (paper/results/figure_1/output/windows/): baboon_windows.csv, baboon_windows_summary.csv, baboon_tissue_timed.csv
Legacy: legacy/desync/validation/baboon/scripts/a40_window.py -> output/a40_{windows,summary}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _internal_time as T  # noqa: E402
import _shared as sh  # noqa: E402

FOLDER = "windows"
SIZES = (3, 4, 5, 6)


def offsets_vs_timed(off_h: pd.Series, timed: pd.Series, n_perm: int, seed: int = 0) -> dict:
    ok = off_h.notna() & np.isfinite(off_h.to_numpy(float))
    if ok.sum() < 6:
        return dict(n=int(ok.sum()))
    r = cs.timed_agreement(off_h[ok].to_numpy(float) * sh.W, timed[ok].to_numpy(float), n_perm=n_perm, seed=seed)
    return dict(r_fl=r.r_fl, p_fl=r.p_fl, r_js=r.r_js, p_js=r.p_js, mae_h=r.median_abs_err_h, n=int(r.n))


def fit_window(pb, kept, timed, animals, zt_unwrapped, frame_F, label, n_perm, seed):
    limit_threads(1)
    row = dict(window=label, m=len(animals), start_zt=zt_unwrapped[0], animals=",".join(animals),
               c_true=float(np.mod(np.mean(zt_unwrapped), 24)))
    sub = pb[pb.donor.isin(animals)].reset_index(drop=True)
    try:
        ds = cs.Dataset(sub, anchors=kept, min_donors=min(len(animals), 5), name=label)
        e, ft = T.host_estimate(ds, kept)
    except Exception as ex:  # recorded as infeasible, as legacy
        return row | dict(status=f"infeasible: {type(ex).__name__}: {str(ex)[:80]}")
    row |= dict(status="ok", kappa=e["kappa"], r_axis=e["r"], w_est=e["w"], shift_h=e["shift_h"], c_curv=e["c_curv"])
    ztmap = dict(zip(animals, zt_unwrapped))
    t = e["tangent"]
    rs = sps.spearmanr(t.to_numpy(), [ztmap[a] for a in t.index]).statistic
    row |= dict(spearman=rs, abs_spearman=abs(rs), flip_ok=bool(rs > 0))
    row["c_est_frame"] = float(np.mod(e["c_curv"] - frame_F, 24))
    row["err_centre_frame"] = float(T.cdiff(row["c_est_frame"], row["c_true"]))
    row["err_centre_noframe"] = float(T.cdiff(e["c_curv"], row["c_true"]))
    row["flip_ok_centre"] = bool(abs(row["err_centre_frame"]) < 6)
    if len(animals) >= 5:  # tissue offsets need >= 5 rows per tissue (`Dataset.fit`)
        off = ft.table.offset_h.reindex(kept)
        row |= {f"loo_{k}": v for k, v in offsets_vs_timed(off, timed, n_perm, seed).items()}
    else:
        row["loo_status"] = "infeasible: tissue offsets need >= 5 animals"
    return row


def summarise(D: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for m, g in D.groupby("m"):
        ok = g[g.status == "ok"]
        r = dict(m=m, span_h=2 * (m - 1), n_windows=len(g), n_ok=len(ok), flips_ok=int(ok.flip_ok.sum()) if len(ok) else 0,
                 flips_ok_centre=int(ok.flip_ok_centre.sum()) if len(ok) else 0)
        for c in ("loo_r_fl", "loo_r_js", "loo_mae_h", "spearman", "abs_spearman", "err_centre_frame", "err_centre_noframe"):
            if c in ok and ok[c].notna().any():
                x = ok[c].abs() if c.startswith("err") else ok[c]
                suf = "_abs" if c.startswith("err") else ""
                r[f"{c}{suf}_med"], r[f"{c}{suf}_min"], r[f"{c}{suf}_max"] = x.median(), x.min(), x.max()
        rows.append(r)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None, help="permutations of the offset correlations (default 2000; 100 with --subsample)")
    ap.add_argument("--sizes", type=int, nargs="+", default=list(SIZES), help="window sizes (animals)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = args.n_perm if args.n_perm is not None else (100 if run["quick"] else 2000)
    pb, donors = sh.baboon()
    kept = sorted(pb.compartment.unique())
    animals = sorted(pb.donor.unique())
    assert len(animals) == 12
    tt = sh.timed_tissue_gains(pb, donors)
    tt.to_csv(out / "baboon_tissue_timed.csv")
    timed = tt.timed_phase.reindex(kept)
    # species frame F from the full-day fit (mean-anchor frame)
    ds = cs.Dataset(pb, anchors=kept, name="baboon_full")
    e, ft = T.host_estimate(ds, kept)
    f = ft.host[["x", "y"]].to_numpy()
    f = f - f.mean(0)
    phi = np.angle(f[:, 0] + 1j * f[:, 1]) + e["shift_h"] * sh.W
    zt = donors.zt.reindex(ft.host.index).to_numpy(float)
    F = float(np.mod(np.angle(np.exp(1j * (phi - zt * sh.W)).mean()) / sh.W, 24))
    full = offsets_vs_timed(ft.table.offset_h.reindex(kept), timed, n_perm)
    print(f"full-day frame F = {F:.2f} h; c_curv {e['c_curv']:.2f} h, r_axis {e['r']:.2f}; full-day offsets vs timed: "
          f"{ {k: round(v, 3) for k, v in full.items()} }", flush=True)
    jobs = []
    for m in args.sizes:
        for s in range(12 if args.subsample is None else min(12, max(2, args.subsample // 6))):
            an = [animals[(s + i) % 12] for i in range(m)]
            ztu = [2 * s + 2 * i for i in range(m)]
            jobs.append((an, ztu, f"m{m}_ZT{2 * s:02d}-ZT{(2 * s + 2 * (m - 1)) % 24:02d}"))
    rows = Parallel(n_jobs=args.jobs)(delayed(fit_window)(pb, kept, timed, an, ztu, F, lab, n_perm, args.seed)
                                      for an, ztu, lab in jobs)
    D = pd.DataFrame(rows)
    D.insert(1, "full_day_F_h", F)
    D.to_csv(out / "baboon_windows.csv", index=False)
    S = summarise(D)
    S.to_csv(out / "baboon_windows_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(S.round(3).T.to_string())
    ok = D[(D.m == 5) & (D.status == "ok")]
    if len(ok) and "loo_r_fl" in ok:
        p1, p2s, p2f = ok.loo_r_fl.median(), ok.abs_spearman.median(), int(ok.flip_ok.sum())
        print(f"PRIMARY (1) m=5 median Fisher-Lee r = {p1:.3f} (> 0.5): {'PASSES' if p1 > 0.5 else 'FAILS'}")
        print(f"PRIMARY (2) m=5 median |Spearman| = {p2s:.3f} (> 0.7), flips correct {p2f}/{len(ok)} (>= 9 of 12): "
              f"{'PASSES' if p2s > 0.7 and p2f >= 9 else 'FAILS'}")


if __name__ == "__main__":
    main()
