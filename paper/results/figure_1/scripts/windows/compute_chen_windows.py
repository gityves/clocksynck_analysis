"""V25.3: frame-transfer test of the A25 internal-time estimator in Chen et al. 2016 prefrontal cortex (BA11, BA47; ZT of death known).

Pre-registered as A25 (2026-10-14; docs/PLAN_internal_desync.md (archive)): "V25.3, Chen 2016 prefrontal cortex (2 regions, 146 donors, microarray, ZT of death):
truth in clock time = ZT + 6.5 h (mean local sunrise); windows as V25.2 (centres 0, 3, ..., 21 h; widths 6 / 9 / 12 / 18 h; plus all donors), >= 20 donors
per region; this is the frame-transfer test (independent of GTEx). Acceptance: median |curvature-centre error| <= 3 h over windows with r < 0.6." Run note: V25.3 fails.
Host: the two regions (a 2-compartment host, `compute_breast_timed.two_host_estimate`: local, not in clocksync). Array input: `Dataset` takes counts, so the log2
intensities of the 14 clock genes (`GSE71620`, as `organs/compute_chen_age.py`) are written as pseudo-counts 1e6 exp(ln2 (y - median_region(y))) with
a constant depth of 1e8: every gene has the same median count (equal gene weights, as the legacy projection of centred log expression) and the Poisson variance is
negligible (the legacy states had V = 0). Estimands: axis centre mod 12 h, curvature centre c_curv, hospital-prior centre, radial / tangential ratio r,
width w, per-donor internal time vs recorded time (median absolute error, circular correlation), flip agreement. Null: none (descriptive; acceptance rule above).

    python compute_chen_windows.py [--subsample N]
Outputs: paper/results/figure_1/output/windows/a25_chen.csv, a25_chen_acceptance.csv
Legacy: legacy/desync/common/a25_internal_time.py (chen_step) -> common_output/a25_chen.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

import _internal_time as IT  # noqa: E402
from compute_breast_timed import two_host_estimate  # noqa: E402

FOLDER = "windows"
RAW = DATA / "brain_tod" / "GSE71620"
KEY = ["age", "sex", "race", "pmi", "ph", "tod_raw"]
SUNRISE = 6.5  # mean local sunrise (Pittsburgh): clock time = ZT + 6.5 h
CENTRES, WIDTHS = range(0, 24, 3), (6, 9, 12, 18)
REGIONS = ["BA11", "BA47"]


def load() -> tuple[pd.DataFrame, pd.Series, cs.Waveform]:
    """Pseudo-count pseudobulk (donor, region) of the 14 clock genes and the ZT of death per donor (donor ids rebuilt as in compute_chen_age)."""
    s = pd.read_csv(RAW / "samples.csv")
    s = s[s.tod_raw.notna()].copy()
    s["donor"] = s[KEY].astype(str).agg("|".join, axis=1)
    s["donor"] = s.donor.map({d: f"d{i:03d}" for i, d in enumerate(sorted(s.donor.unique()))})
    e = pd.read_csv(RAW / "expression_clock.csv").set_index("gsm").loc[s.gsm]
    e.columns = [canonical(c.split("|")[0]) for c in e.columns]
    wf = cs.gtex()
    genes = [g for g in wf.genes if g in e.columns]
    y = e[genes].T.groupby(level=0).mean().T[genes].reset_index(drop=True)
    s = s.reset_index(drop=True)
    y = y - y.groupby(s.region.to_numpy()).transform("median")
    cnt = 1e6 * np.exp(np.log(2.0) * y)
    pb = pd.concat([s[["donor"]].assign(compartment=s.region.astype(str), depth=1e8), cnt], axis=1)
    return pb, s.groupby("donor").tod_raw.first(), wf.subset(genes)


def window(ds_all, keep, t_true, row, min_donors=20, min_comp_donors=20):
    dd = ds_all.subset(donors=keep)
    cnt = dd.table.groupby("compartment").resection.nunique()
    row = dict(row, n=len(dd.donor_ids), n_anchors=int((cnt.reindex(REGIONS) >= min_comp_donors).sum()))
    if row["n"] < min_donors or row["n_anchors"] < 2:
        return row | dict(status="too few")
    try:
        e, host = two_host_estimate(dd, REGIONS)
    except Exception as ex:
        return row | dict(status=f"error: {type(ex).__name__}: {ex}")
    t = t_true.reindex(e["donor_h"].index)
    tw = IT.true_window(t_true.reindex(host.index).to_numpy())
    est = e["donor_h"].to_numpy()
    sa = np.sin(IT.W * est - np.angle(np.exp(1j * IT.W * est).mean()))
    sb = np.sin(IT.W * t.to_numpy() - np.angle(np.exp(1j * IT.W * t.to_numpy()).mean()))
    return row | {k: v for k, v in e.items() if k not in ("donor_h", "tangent", "donor_day_h", "tangent_day", "sigma_f")} | tw | dict(
        status="ok", err_axis12=IT.cdiff(2 * e["c_axis12"], 2 * np.mod(tw["c_true"], 12.0)) / 2,
        err_curv=IT.cdiff(e["c_curv"], tw["c_true"]), err_day=IT.cdiff(e["c_day"], tw["c_true"]),
        donor_mae=float(np.median(np.abs(IT.cdiff(est, t.to_numpy())))),
        donor_ccorr=float((sa * sb).sum() / np.sqrt((sa ** 2).sum() * (sb ** 2).sum())),
        flip_ok=bool(abs(IT.cdiff(e["c_curv"], tw["c_true"])) < 6))


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    limit_threads(1)
    out = out_dir(FOLDER, run["quick"])
    pb, zt, wf = load()
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[: args.subsample]
        pb, zt = pb[pb.donor.isin(keep)].reset_index(drop=True), zt.loc[keep]
    ds = cs.Dataset(pb, anchors=REGIONS, waveform=wf, name="chen_brain")
    t = pd.Series(IT.wrap24(zt.to_numpy(float) + SUNRISE), index=zt.index)
    sets = [("all", np.nan, np.nan, t.index.to_numpy())]
    for c in CENTRES:
        for w in WIDTHS:
            sets.append((f"{c:g}h+-{w / 2:g}", c, w, t.index[np.abs(IT.cdiff(t.to_numpy(), c)) <= w / 2].to_numpy()))
    mind = 8 if args.subsample else 20
    rows = [window(ds, k, t, dict(data="chen_brain", window=n, c_set=c, w_set=w), mind, mind) for n, c, w, k in sets]
    D = pd.DataFrame(rows)
    D["flip_ok"] = np.abs(D.get("err_curv", np.nan)) < 6
    D.to_csv(out / "a25_chen.csv", index=False)
    ok = D[(D.status == "ok") & (D.r < IT.R_MAX)]
    acc = pd.DataFrame([dict(windows_ok=int((D.status == "ok").sum()), windows_r_lt_0_6=len(ok),
                             median_abs_err_curv_h=float(ok.err_curv.abs().median()) if len(ok) else np.nan,
                             median_abs_err_axis12_h=float(ok.err_axis12.abs().median()) if len(ok) else np.nan,
                             flip_ok=int(ok.flip_ok.sum()),
                             accept=bool(len(ok) and ok.err_curv.abs().median() <= 3))])
    acc.to_csv(out / "a25_chen_acceptance.csv", index=False)
    pd.set_option("display.width", 250)
    print(D[D.status == "ok"][["window", "n", "c_true", "c_axis12", "c_curv", "err_axis12", "err_curv", "r", "r_true", "w", "w_true", "donor_mae", "flip_ok"]]
          .round(2).to_string(index=False))
    print(acc.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
