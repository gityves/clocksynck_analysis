"""Diagnostic (post hoc, R1.2 DIP figure): where and why CHIRAL's donor internal phase (DIP) and clocksync's internal time differ.

Same donors, tissues and fits as `gtex_dip_compute.py`, design `full` (clocksync on all fast deaths). Per-donor values stay in
memory; written are aggregates over groups of >= 5 donors and summary statistics:
- `by_dip_bin.csv`: per 3 h bin of DIP: recorded time of death (circular mean, resultant length), clocksync host strength
  (radius of the donor's posterior host mean, relative to the cohort median), tissues per donor, CHIRAL samples per donor,
  and the median absolute error of each phase against the recorded time (each aligned to it by `evaluate.align`).
- `by_tod_bin.csv`: the same errors per 3 h bin of recorded time of death.
- `by_strength.csv`: both phases against the recorded time within tertiles of clocksync host strength.
- `agreement.csv`: |circular r| and median absolute error of each phase against the recorded time, with paired donor bootstrap
  95 % intervals of the differences (1000 replicates); fraction of donors with an error > 6 h.

    paper/common/run.sh paper/presentation/scripts/gtex_dip_diag.py
Outputs: paper/presentation/output/gtex_dip/diag/*.csv
"""
from __future__ import annotations

import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd

import gtex_dip_compute as G

ev, A, B, T = G.ev, G.A, G.B, G.T
OUT = G.OUT / "diag"
MIN_N, N_BOOT = G.MIN_CELL, 1000


def errors(x, t):
    """Absolute error (h) of phase x against recorded time t after the best rotation and reflection."""
    xa, _, _ = ev.align(x, t, step=0.1)
    return np.abs(ev.wrap(xa - t))


def agg(g: pd.DataFrame) -> dict:
    n = len(g)
    if n < MIN_N:
        return dict(n=n)
    z = np.exp(1j * G.W * g.tod.to_numpy()).mean()
    return dict(n=n, tod_mean_h=np.mod(np.angle(z) / G.W, 24), tod_R=abs(z), host_strength_median=g.strength.median(),
                tissues_median=g.n_tissues.median(), chiral_samples_median=g.n_chiral.median(),
                it_err_median_h=g.err_it.median(), dip_err_median_h=g.err_dip.median(),
                it_err_gt6h=(g.err_it > 6).mean(), dip_err_gt6h=(g.err_dip > 6).mean())


def stats(x, d, t) -> dict:
    return dict(r_it=abs(ev.circ_corr(x, t)), r_dip=abs(ev.circ_corr(d, t)),
                mae_it=np.median(errors(x, t)), mae_dip=np.median(errors(d, t)))


def main():
    G.cs.config.limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    a = ad.read_h5ad(G.LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    a.file.close()
    o["tissue"] = o.tissue.astype(str)
    rows = []
    for tis in A.window_tissues(o).index:
        counts, tod = B.S.load([tis])[tis]
        tod = tod.dropna()
        rows.append(B.clock_table({tis: (counts.loc[tod.index], tod)}))
    table = pd.concat(rows, ignore_index=True)
    print(f"loaded ({time.time() - t0:.0f} s)", flush=True)
    anchors = sorted(table.compartment.unique())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = G.cs.Dataset(table.drop(columns="tod"), anchors=anchors, waveform=G.cs.gtex())
        e, ft = T.host_estimate(ds, anchors)
    host = ft.host[["x", "y"]]
    radius = np.hypot(*(host - host.mean()).to_numpy().T)
    meta = pd.read_csv(A.META, sep=" ").dropna(subset=["CT"])
    n_chiral = meta.fullID.str.split(".").str[:2].str.join("-").value_counts()

    dip = G.dip_per_donor()
    tod = table.groupby("donor").tod.first()
    d = pd.DataFrame(dict(it=e["donor_h"])).dropna()
    d = d.loc[d.index.intersection(dip.index)]
    d["dip"], d["tod"] = dip.loc[d.index], tod.loc[d.index]
    d["strength"] = pd.Series(radius, index=host.index).reindex(d.index)
    d["strength"] /= d.strength.median()
    d["n_tissues"] = table.groupby("donor").compartment.nunique().reindex(d.index)
    d["n_chiral"] = n_chiral.reindex(d.index).fillna(0)
    d["err_it"], d["err_dip"] = errors(d.it, d.tod), errors(d.dip, d.tod)

    edges = np.arange(0, 25, 3.0)
    for col, name in (("dip", "by_dip_bin"), ("tod", "by_tod_bin")):
        b = pd.cut(d[col], edges, right=False)
        pd.DataFrame([dict(lo=k.left, hi=k.right, **agg(g)) for k, g in d.groupby(b, observed=False)]).to_csv(
            OUT / f"{name}.csv", index=False)
    terc = pd.qcut(d.strength, 3, labels=["weak", "middle", "strong"])
    rows = []
    for k, g in d.groupby(terc, observed=True):
        rows.append(dict(tertile=k, n=len(g), strength_median=g.strength.median(), **stats(g.it, g.dip, g.tod),
                         it_err_gt6h=(g.err_it > 6).mean(), dip_err_gt6h=(g.err_dip > 6).mean()))
    pd.DataFrame(rows).to_csv(OUT / "by_strength.csv", index=False)

    point = stats(d.it, d.dip, d.tod)
    rng = np.random.default_rng(0)
    reps = []
    for _ in range(N_BOOT):
        g = d.iloc[rng.integers(0, len(d), len(d))]
        s = stats(g.it.to_numpy(), g.dip.to_numpy(), g.tod.to_numpy())
        reps.append((s["r_dip"] - s["r_it"], s["mae_dip"] - s["mae_it"]))
    reps = np.array(reps)
    out = dict(n=len(d), **point, r_dip_minus_it=point["r_dip"] - point["r_it"],
               r_diff_lo=np.quantile(reps[:, 0], .025), r_diff_hi=np.quantile(reps[:, 0], .975),
               mae_dip_minus_it=point["mae_dip"] - point["mae_it"],
               mae_diff_lo=np.quantile(reps[:, 1], .025), mae_diff_hi=np.quantile(reps[:, 1], .975),
               it_err_gt6h=(d.err_it > 6).mean(), dip_err_gt6h=(d.err_dip > 6).mean(),
               corr_strength_err_it=d[["strength", "err_it"]].corr("spearman").iloc[0, 1],
               corr_strength_err_dip=d[["strength", "err_dip"]].corr("spearman").iloc[0, 1])
    pd.Series(out).to_csv(OUT / "agreement.csv", header=["value"])
    print(pd.Series(out).round(3).to_string(), f"\ndone in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
