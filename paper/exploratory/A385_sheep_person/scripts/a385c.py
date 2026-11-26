"""A385c (descriptive check of A385b, registered 2026-11-26): does the glucocorticoid-target state explain the sheep
shared clock loudness?

A385b D1 leverage (A385 halvings, within project with >= 30 animals, Fisher-combined), radii residualised on the phase
harmonics and log tissue counts (A385) plus: (A) the animal's GC and IEG scores; (B) both halves' GC and IEG scores.
Scores as A386 (non-clock genes; mean log fraction centred within tissue x project, standardised within tissue).
Clock (14 genes, core) with a sheep-permuted null; 100 companion sets adjusted the same way.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A385_sheep_person/scripts/a385c.py --jobs 6
Outputs (`../output/`): c_leverages.csv, c_summary.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a385_a386 as M  # noqa: E402
import a385b as R  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
VARIANTS = ("unadjusted", "A: animal GC + IEG", "B: half GC + IEG")
N_SETS, N_PERM = 50, 1000
S = {}


def row_scores(pb, d, genes):
    """Per (animal, tissue) row of the Dataset: A386 score (centred within tissue x project, standardised within tissue)."""
    key = pd.Series(np.arange(len(pb)), index=pd.MultiIndex.from_arrays([pb.donor, pb.compartment]))
    rows = key.loc[list(zip(np.asarray(d.resection).astype(str), d.compartment))].to_numpy()
    y = np.log((pb[genes].to_numpy(float)[rows] + 0.5) / pb.depth.to_numpy(float)[rows, None]).mean(1)
    comp = np.asarray(d.compartment)
    grp = pd.Series(comp.astype(str)) + "|" + pd.Series(np.asarray(d.group).astype(str))
    y = y - pd.Series(y).groupby(grp.to_numpy()).transform("mean").to_numpy()
    z = pd.Series(y).groupby(comp).transform(lambda v: (v - v.mean()) / v.std(ddof=1)).to_numpy()
    return pd.DataFrame({"donor": np.asarray(d.resection).astype(str), "comp": comp, "score": z})


def states(w):
    """Per halving: ids, states, log counts, animal GC / IEG, half GC / IEG (A, B)."""
    pb = S["pb"]
    ds = cs.Dataset(pb[M.META + list(w.genes)], anchors=S["per"], waveform=w, mixed_strata=True, name=w.name,
                    min_donors=5)
    d = ds.data
    sc = row_scores(pb, d, S["gc"]).rename(columns={"score": "gc"})
    sc["ieg"] = row_scores(pb, d, S["ieg"]).score.to_numpy()
    animal = sc.groupby("donor")[["gc", "ieg"]].mean()
    out = []
    for (sa, sb), (j, a, b, cov) in zip(M.halvings(S["per"]), R.states(ds, M.halvings(S["per"]))):
        half = [sc[sc.comp.isin(s)].groupby("donor")[["gc", "ieg"]].mean().reindex(j).to_numpy(float) for s in (sa, sb)]
        out.append((j, a, b, cov, animal.reindex(j).to_numpy(float), np.column_stack(half)))
    return out


def d1(st, variant, rng=None):
    """A385 D1 leverage with the variant's extra covariates; rng: permute B within project (null)."""
    L = []
    for j, a, b, cov, an, hf in st:
        extra = {"unadjusted": None, "A: animal GC + IEG": an, "B: half GC + IEG": hf}[variant]
        X = cov if extra is None else np.column_stack([cov, extra])
        ok = np.isfinite(X).all(1)
        pr = S["project"].reindex(j).to_numpy()
        z = w = 0.0
        for p in np.unique(pr):
            m = np.flatnonzero((pr == p) & ok)
            if len(m) < M.MIN_PROJECT:
                continue
            bm = b[m] if rng is None else b[m][rng.permutation(len(m))]
            z += (len(m) - 3) * np.arctanh(M.leverage(a[m], bm, X[m]))
            w += len(m) - 3
        L.append(np.tanh(z / w) if w > 0 else np.nan)
    return float(np.nanmedian(L))


def task(job):
    warnings.filterwarnings("ignore")
    kind, i, rd, w = job
    try:
        st = states(w)
        out = [dict(kind=kind, set=i, readout=rd, variant=v, leverage=d1(st, v)) for v in VARIANTS]
        if kind == "clock":
            for v in VARIANTS:
                rng = np.random.default_rng([385, 5])
                null = np.array([d1(st, v, rng) for _ in range(N_PERM)])
                lev = [o for o in out if o["variant"] == v][0]["leverage"]
                [o for o in out if o["variant"] == v][0].update(p_perm=float((1 + np.sum(null >= lev)) / (1 + N_PERM)))
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
    pb, per, _ = M.load()
    clock = {canonical(g) for g in M.WF.genes}
    gc = [g for g in dict.fromkeys(GC_TARGETS) if g in pb.columns and g not in {"NFIL3", "BHLHE40"} and canonical(g) not in clock]
    ieg = [g for g in dict.fromkeys(IMMEDIATE_EARLY) if g in pb.columns and canonical(g) not in clock]
    lm = M.log_mean(pb)
    sets = {k: pseudo_waveforms(lm, N_SETS, np.random.default_rng([385, t]), M.WF, kind=k)
            for k, t in (("matched", 2), ("responsive", 3))}
    jobs = [("clock", 0, "14 genes", M.WF), ("clock", 0, "core", M.CORE)]
    jobs += [(k, j, rd, w if rd == "14 genes" else R.core_of(w, f"core{j}"))
             for k, ws in sets.items() for j, w in enumerate(ws) for rd in ("14 genes", "core")]
    genes = sorted({g for *_, w in jobs for g in w.genes} | set(gc) | set(ieg))
    S.update(pb=pb[M.META + genes].copy(), per=per, gc=gc, ieg=ieg, project=pb.groupby("donor").stratum.first())
    del pb
    print(f"{len(jobs)} waveforms; GC {len(gc)} IEG {len(ieg)} genes ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows += r
            if (k + 1) % 10 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    C = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    C.to_csv(OUT / "c_leverages.csv", index=False)
    summ = []
    for rd in ("14 genes", "core"):
        for v in VARIANTS:
            c = C[(C.readout == rd) & (C.variant == v)]
            clk = c[c.kind == "clock"].iloc[0]
            q = {k: float(np.nanquantile(c[c.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}
            comp = c[c.kind != "clock"].leverage
            summ.append(dict(readout=rd, variant=v, clock=clk.leverage, p_perm=clk.p_perm, matched_q95=q["matched"],
                             responsive_q95=q["responsive"], companion_median=float(comp.median()),
                             pct_companions_below=float(100 * np.mean(comp.dropna() < clk.leverage)),
                             above_both=bool(clk.leverage > max(q.values()))))
    Sm = pd.DataFrame(summ)
    Sm.to_csv(OUT / "c_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(Sm.round(4).to_string(index=False), f"\nfailed: {int(C.variant.isna().sum())} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
