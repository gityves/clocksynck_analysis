"""Slide data (results, R1.2): CHIRAL's donor internal phase (DIP) against clocksync's untimed internal time, GTEx fast deaths.

Descriptive, post hoc (presentation figure, not an amendment). Donors and tissues as A123 / A123b: GTEx v10 fast deaths
(Hardy 1-2), the 40 tissues with >= 20 donors who died 08-18 h, every tissue an anchor (free host).
- `full`: clocksync on all fast deaths, GTEx waveform (as A123b design F); internal time = host phase (`donor_h`, A25 estimator).
- `window`: clocksync on the donors who died 08-18 h only, waveform cross-fitted on the fast deaths outside the window (as A123);
  internal time = host phase with the 12 h flip resolved by the hospital prior (`donor_day_h`).
DIP: CHIRAL's published per-sample CT (naef-lab/CHIRAL commit 1f78a9a, `sample_metadata.txt`), one value per donor (circular mean
over its samples). CHIRAL learnt it from all GTEx donors, whole day. clocksync's origin and direction are arbitrary: its internal
time is aligned to the DIP by the rotation and reflection minimising the median absolute difference (`evaluate.align`).

Controlled access: per-donor phases stay in memory. Written: per-bin aggregates (bins of 3 h of DIP, bins with < 5 donors blanked),
a 2D histogram (3 h x 3 h cells, cells with < 5 donors blanked) and summary statistics.

    paper/common/run.sh paper/presentation/scripts/gtex_dip_compute.py
Outputs: paper/presentation/output/gtex_dip/{summary,bins,hist2d}.csv
"""
from __future__ import annotations

import importlib.util
import sys
import time
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("a123_wh", PAPER / "exploratory/A123_working_hours/scripts/compute_working_hours.py")
A = importlib.util.module_from_spec(_spec)
sys.modules["a123_wh"] = A
_spec.loader.exec_module(A)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from common.gtex import LABELLED  # noqa: E402

B, T = A.B, A.T
OUT = PAPER / "presentation" / "output" / "gtex_dip"
BIN_H, MIN_CELL = 3.0, 5
W = 2 * np.pi / 24


def dip_per_donor() -> pd.Series:
    meta = pd.read_csv(A.META, sep=" ").dropna(subset=["CT"])
    z = np.exp(1j * W * meta.CT.to_numpy(float))
    donor = meta.fullID.str.split(".").str[:2].str.join("-")  # GTEX.XXXXX.<sample> -> GTEX-XXXXX
    d = pd.Series(z, index=donor.to_numpy()).groupby(level=0).mean()
    return pd.Series(np.mod(np.angle(d.to_numpy()) / W, 24), index=d.index, name="dip")


def internal_time(table: pd.DataFrame, wf, key: str) -> tuple[pd.Series, dict]:
    anchors = sorted(table.compartment.unique())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(table.drop(columns="tod"), anchors=anchors, waveform=wf)
        e, _ = T.host_estimate(ds, anchors)
    return e[key].dropna(), dict(window_width_h=e["w"], centre_day_h=e["c_day"])


def summarise(design: str, it: pd.Series, dip: pd.Series, tod: pd.Series, extra: dict):
    k = it.index.intersection(dip.index)
    x, d, t = it.loc[k].to_numpy(), dip.loc[k].to_numpy(), tod.reindex(k).to_numpy()
    xa, mae, sign = ev.align(x, d)
    diff = ev.wrap(xa - d)
    row = dict(design=design, n_donors=len(k), circ_corr_dip=ev.circ_corr(xa, d), median_abs_diff_h=mae, reflected=sign < 0,
               circ_corr_it_tod=abs(ev.circ_corr(x, t)), circ_corr_dip_tod=abs(ev.circ_corr(d, t)),
               dip_mae_tod_h=ev.align(d, t)[1], it_mae_tod_h=ev.align(x, t)[1], **extra)
    edges = np.arange(0, 24 + BIN_H, BIN_H)
    b = np.digitize(d, edges[1:-1])
    bins = []
    for j in range(len(edges) - 1):
        m = b == j
        n = int(m.sum())
        r = dict(design=design, dip_lo=edges[j], dip_hi=edges[j + 1], n=n)
        if n >= MIN_CELL:
            z = np.exp(1j * W * xa[m]).mean()
            r.update(it_mean_h=np.mod(np.angle(z) / W, 24), it_circ_sd_h=np.sqrt(-2 * np.log(abs(z))) / W,
                     diff_median_h=float(np.median(diff[m])), diff_q25_h=float(np.quantile(diff[m], .25)),
                     diff_q75_h=float(np.quantile(diff[m], .75)))
        bins.append(r)
    h, _, _ = np.histogram2d(d, xa, bins=[edges, edges])
    hist = pd.DataFrame(h, index=pd.Index(edges[:-1], name="dip_lo"), columns=edges[:-1]).stack().rename("n").reset_index()
    hist.columns = ["dip_lo", "it_lo", "n"]
    hist["n"] = hist.n.where(hist.n >= MIN_CELL)
    return row, bins, hist.assign(design=design)


def main():
    cs.config.limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    a.file.close()
    o["tissue"] = o.tissue.astype(str)
    tissues = list(A.window_tissues(o).index)
    rows = []
    for t in tissues:
        counts, tod = B.S.load([t])[t]
        tod = tod.dropna()
        rows.append(B.clock_table({t: (counts.loc[tod.index], tod)}))
        print(f"loaded {t} ({time.time() - t0:.0f} s)", flush=True)
    table = pd.concat(rows, ignore_index=True)
    tod = table.groupby("donor").tod.first()
    dip = dip_per_donor()

    win = sorted(table.donor[in_window(table.tod)].unique())
    wtab = table[table.donor.isin(win) & in_window(table.tod)]
    designs = {"full": (table, cs.gtex(), "donor_h"),
               "window": (wtab, B.cross_fit_waveform(table, win), "donor_day_h")}
    summ, bins, hist = [], [], []
    for name, (tab, wf, key) in designs.items():
        it, extra = internal_time(tab, wf, key)
        s, b, h = summarise(name, it, dip, tod, extra)
        summ.append(s), bins.extend(b), hist.append(h)
        print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()}, f"{time.time() - t0:.0f} s", flush=True)
    pd.DataFrame(summ).to_csv(OUT / "summary.csv", index=False)
    pd.DataFrame(bins).to_csv(OUT / "bins.csv", index=False)
    pd.concat(hist).to_csv(OUT / "hist2d.csv", index=False)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
