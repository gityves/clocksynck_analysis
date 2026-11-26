"""A400b (descriptive check of A400, registered 2026-11-26): is the mouse shared clock loudness a time-of-day, sampling
or technical effect?

A400 leverage (20 halvings of the 16 peripheral tissues; radii residualised on two harmonics of both half-host phases,
log tissue counts and the arm), with in addition: V1 harmonics 1-4 of the recorded ZT; V2 ZT x arm cell indicators;
V3 the mouse's and each half's median log depth; V4 V2 + V3. Clock with a mice-permuted null (within arm); 100 companion
sets adjusted the same way. Descriptive: R^2 of the half-host radii on the ZT x arm cells.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A400_mouse_person/scripts/a400b.py --jobs 6
Outputs (`../output/`): b_leverages.csv, b_summary.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
for p in (PAPER, HERE.parent, EXP / "A385_sheep_person" / "scripts", EXP / "A394_mouse_organisation" / "scripts"):
    sys.path.insert(0, str(p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a385b as R  # noqa: E402
import a394 as M  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
VARIANTS = ("V0 unadjusted", "V1 ZT harmonics 1-4", "V2 ZT x arm cells", "V3 depth", "V4 cells + depth")
N_SETS, N_PERM = 50, 1000
S = {}


def states(w):
    """Per halving: mice, states, log counts, log median depth of each half."""
    pb = S["mouse"]
    ds = cs.Dataset(pb[M.KEYS + list(w.genes)], anchors=S["per"], waveform=w, name=w.name, min_donors=5)
    halvings = R.halves(S["per"], [400, 1])
    out = []
    for (sa, sb), (j, a, b, cov) in zip(halvings, R.states(ds, halvings)):
        dep = [np.log(pb[pb.compartment.isin(s)].groupby("donor").depth.median().reindex(j).to_numpy(float))
               for s in (sa, sb)]
        out.append((j, a, b, cov, np.column_stack(dep)))
    return out


def covariates(j, cov, dep, variant):
    """(extra covariate matrix, group labels for indicators)."""
    zt = S["zt"].reindex(j).to_numpy(float) * 2 * np.pi / 24
    harm = np.column_stack([g(k * zt) for k in range(1, 5) for g in (np.cos, np.sin)])
    dd = np.column_stack([S["depth"].reindex(j).to_numpy(float), dep])
    cell = S["cell"].reindex(j).to_numpy()
    arm = S["arm"].reindex(j).to_numpy()
    return {"V0 unadjusted": (cov, arm), "V1 ZT harmonics 1-4": (np.column_stack([cov, harm]), arm),
            "V2 ZT x arm cells": (cov, cell), "V3 depth": (np.column_stack([cov, dd]), arm),
            "V4 cells + depth": (np.column_stack([cov, dd]), cell)}[variant]


def leverage(st, variant, rng=None):
    L = []
    for j, a, b, cov, dep in st:
        X, grp = covariates(j, cov, dep, variant)
        if rng is not None:
            arm = S["arm"].reindex(j).to_numpy()
            b = b.copy()
            for p in np.unique(arm):
                m = np.flatnonzero(arm == p)
                b[m] = b[rng.permutation(m)]
        L.append(R.pooled_leverage(a, b, X, grp)[0])
    return float(np.nanmedian(L))


def radius_r2(st):
    """Median over halvings of the R^2 of each half-host radius on the ZT x arm cells."""
    out = []
    for j, a, b, _, _ in st:
        cell = pd.get_dummies(pd.Series(S["cell"].reindex(j).to_numpy())).to_numpy(float)
        for z in (a, b):
            r = np.abs(z)
            res = r - cell @ np.linalg.lstsq(cell, r, rcond=None)[0]
            out.append(1 - res.var() / r.var())
    return float(np.median(out))


def task(job):
    warnings.filterwarnings("ignore")
    kind, i, rd, w = job
    try:
        st = states(w)
        out = [dict(kind=kind, set=i, readout=rd, variant=v, leverage=leverage(st, v)) for v in VARIANTS]
        if kind == "clock":
            for o in out:
                rng = np.random.default_rng([400, 4])
                null = np.array([leverage(st, o["variant"], rng) for _ in range(N_PERM)])
                o.update(p_perm=float((1 + np.sum(null >= o["leverage"])) / (1 + N_PERM)), null_q95=float(np.quantile(null, 0.95)),
                         radius_r2_cells=radius_r2(st))
        return out
    except Exception as e:  # counted, not dropped silently
        return [dict(kind=kind, set=i, readout=rd, variant=None, error=f"{type(e).__name__}: {e}"[:200])]


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    pb, per, _ = M.mouse_table()
    pb = pb[pb.compartment.isin(per)].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    lm = M.log_mean(pb)
    sets = {k: pseudo_waveforms(lm, N_SETS, np.random.default_rng([394, t]), R.WF, kind=k)
            for k, t in (("matched", 2), ("responsive", 3))}
    jobs = [("clock", 0, "core", R.CORE), ("clock", 0, "14 genes", R.WF)]
    jobs += [(k, j, rd, w if rd == "14 genes" else R.core_of(w, f"core{j}"))
             for k, ws in sets.items() for j, w in enumerate(ws) for rd in ("14 genes", "core")]
    genes = sorted({g for *_, w in jobs for g in w.genes})
    donors = pb.groupby("donor").stratum.first()
    zt = pd.Series(donors.index.str.extract(r"_ZT(\d+)_", expand=False).astype(float), index=donors.index)
    S.update(mouse=pb[M.KEYS + genes].copy(), per=per, arm=donors, zt=zt,
             cell=pd.Series(donors.to_numpy() + "_" + zt.astype(int).astype(str).to_numpy(), index=donors.index),
             depth=np.log(pb.groupby("donor").depth.median()))
    print(f"{len(jobs)} waveforms; {donors.size} mice, {S['cell'].nunique()} cells ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows += r
            if (k + 1) % 10 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    C = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    C.to_csv(OUT / "b_leverages.csv", index=False)
    summ = []
    for rd in ("core", "14 genes"):
        for v in VARIANTS:
            c = C[(C.readout == rd) & (C.variant == v)]
            clk = c[c.kind == "clock"].iloc[0]
            q = {k: float(np.nanquantile(c[c.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}
            comp = c[c.kind != "clock"].leverage.dropna()
            summ.append(dict(readout=rd, variant=v, clock=clk.leverage, p_perm=clk.p_perm, null_q95=clk.null_q95,
                             matched_q95=q["matched"], responsive_q95=q["responsive"],
                             pct_companions_below=float(100 * np.mean(comp < clk.leverage)),
                             above_both=bool(clk.leverage > max(q.values())), radius_r2_cells=clk.radius_r2_cells))
    Sm = pd.DataFrame(summ)
    Sm.to_csv(OUT / "b_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(Sm.round(4).to_string(index=False), f"\nfailed: {int(C.variant.isna().sum())} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
