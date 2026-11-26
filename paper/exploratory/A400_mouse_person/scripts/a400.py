"""A400 (mouse person loudness; gate calibrated on GTEx thinned to the mouse design), registered 2026-11-26.

Statistic (A385b): per halving of the anchors, the half-host radii of the animals in both halves, residualised on two
phase harmonics of both phases, log tissue count per half and the arm indicator; correlation; median over 20 halvings.
Mouse: Deota 2023 (A394 table), 16 peripheral tissues, arm as stratum. GTEx: fast deaths thinned to the mouse design.

Stages:
  power  GTEx thinned (clock 40 replicates; A336 companions 100 + 100, one replicate each); mouse companions (the A394
         sets, 200 + 200); P_a, P_b and the gate. No mouse clock value.
  check  2 sets per arm, 2 replicates; nothing written.
  score  H400.1 (core) and the 14-gene readout, only if the gate passed.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A400_mouse_person/scripts/a400.py --stage power --jobs 6
Outputs (`../output/`): companions.csv, gtex_thin.csv, power.csv; score: tests.csv. Aggregates only (GTEx: leverages).
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
for p in (PAPER, EXP / "A385_sheep_person" / "scripts", EXP / "A394_mouse_organisation" / "scripts",
          EXP / "A336_two_readouts" / "scripts"):
    sys.path.insert(0, str(p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_person_gtex as G  # noqa: E402
import a385b as R  # noqa: E402
import a394 as M  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
WF, CORE = R.WF, R.CORE
N_REP, N_GTEX_SETS, N_MOUSE_SETS, N_PERM = 40, 100, 200, 1000
S = {}


def leverage(st, arm, rng=None):
    """Median over halvings of the pooled leverage with the arm indicator; rng: permute B within arm (null)."""
    L = []
    for j, a, b, cov in st:
        pr = arm.reindex(j).to_numpy()
        if rng is not None:
            b = b.copy()
            for p in np.unique(pr):
                m = np.flatnonzero(pr == p)
                b[m] = b[rng.permutation(m)]
        L.append(R.pooled_leverage(a, b, cov, pr)[0])
    return float(np.nanmedian(L))


def mouse_states(w):
    pb = S["mouse"]
    ds = cs.Dataset(pb[M.KEYS + list(w.genes)], anchors=S["per"], waveform=w, name=w.name, min_donors=5)
    return R.states(ds, R.halves(S["per"], [400, 1]))


def assign(r):
    rng = np.random.default_rng([400, 2, r])
    tis = S["gtex_tissues"]
    free, rows = sorted(tis.index), []
    for aid, k in S["k"].sort_values(ascending=False).items():
        cand = [d for d in free if len(tis[d]) >= k]
        d = cand[rng.integers(len(cand))]
        free.remove(d)
        rows += [(aid, d, t) for t in rng.choice(tis[d], int(k), replace=False)]
    return pd.DataFrame(rows, columns=["animal", "gtex", "compartment"])


def gtex_states(w, r):
    a = assign(r)
    g = S["gtex"].merge(a, left_on=["donor", "compartment"], right_on=["gtex", "compartment"])
    g = g.drop(columns=["donor", "gtex"]).rename(columns={"animal": "donor"})
    ds = cs.Dataset(g[G.META + list(w.genes)], anchors=S["comps"], waveform=w, name=w.name)
    return R.states(ds, R.halves(S["comps"], [331, 1]))


def task(job):
    warnings.filterwarnings("ignore")
    arm_, kind, i, rd, w = job
    base = dict(arm=arm_, kind=kind, set=i, readout=rd)
    try:
        if arm_ == "mouse":
            lev = leverage(mouse_states(w), S["arm"])
        else:
            r = i if kind == "clock" else 1000 + i + (0 if kind == "matched" else N_GTEX_SETS)
            lev = leverage(gtex_states(w, r), S["arm"])
    except Exception as e:  # counted, not dropped silently
        return dict(base, leverage=np.nan, error=f"{type(e).__name__}: {e}"[:200])
    return dict(base, leverage=lev)


def setup(n_mouse, n_gtex, n_rep):
    t0 = time.time()
    pb, per, _ = M.mouse_table()
    pb = pb[pb.compartment.isin(per)].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    lm = M.log_mean(pb)
    ms = {kind: pseudo_waveforms(lm, n_mouse, np.random.default_rng([394, tag]), WF, kind=kind)
          for kind, tag in (("matched", 2), ("responsive", 3))}
    jm = [(kind, j, rd, w if rd == "14 genes" else R.core_of(w, f"core{j}"))
          for kind, ws in ms.items() for j, w in enumerate(ws) for rd in ("14 genes", "core")]
    genes = sorted({g for *_, w in jm for g in w.genes} | set(WF.genes))
    S["mouse"] = pb[M.KEYS + genes].copy()
    S["per"] = per
    S["arm"] = pb.groupby("donor").stratum.first()
    S["k"] = pb.groupby("donor").compartment.nunique()
    gt, _ = G.load()
    gl = [c for c in gt.columns if c not in G.META]
    x = gt[gl].to_numpy(np.float32)
    glm = pd.Series(np.log((x / gt.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=gl)[x.mean(0) >= 1]
    del x
    R.S["gtex_lm"] = glm[~glm.index.duplicated()]
    jg = R.gtex_sets(n_gtex)
    genes = sorted({g for *_, w in jg for g in w.genes} | set(WF.genes))
    gt["donor"] = gt.donor.astype(str)
    S["gtex"] = gt[G.META + genes].copy()
    del gt
    S["comps"] = sorted(S["gtex"].compartment.unique())
    S["gtex_tissues"] = S["gtex"].groupby("donor").compartment.apply(lambda v: sorted(v.unique()))
    print(f"mouse {S['arm'].size} animals, {len(per)} tissues; GTEx {S['gtex_tissues'].size} donors "
          f"({time.time() - t0:.0f}s)", flush=True)
    jobs = [("gtex", "clock", r, rd, w) for r in range(n_rep) for rd, w in (("14 genes", WF), ("core", CORE))]
    jobs += [("gtex", k, j, rd, w) for k, j, rd, w in jg]
    jobs += [("mouse", k, j, rd, w) for k, j, rd, w in jm]
    return jobs


def run(jobs, n_jobs):
    rows, t0 = [], time.time()
    with Pool(n_jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows.append(r)
            if (k + 1) % 25 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def q95(v):
    return float(np.nanquantile(v, 0.95))


def power(n_jobs, check):
    limit_threads(1)
    t0 = time.time()
    C = run(setup(2 if check else N_MOUSE_SETS, 2 if check else N_GTEX_SETS, 2 if check else N_REP), n_jobs)
    if check:
        print(C.groupby(["arm", "kind", "readout"]).agg(rows=("set", "size"), nan=("leverage", lambda v: v.isna().sum()))
              .to_string(), f"({time.time() - t0:.0f}s)")
        return
    C[C.arm == "mouse"].to_csv(OUT / "companions.csv", index=False)
    C[C.arm == "gtex"].to_csv(OUT / "gtex_thin.csv", index=False)
    rows = []
    for rd in ("14 genes", "core"):
        sel = lambda arm_, kind: C[(C.arm == arm_) & (C.kind == kind) & (C.readout == rd)].leverage  # noqa: E731
        clk = sel("gtex", "clock").to_numpy(float)
        g_thr = max(q95(sel("gtex", "matched")), q95(sel("gtex", "responsive")))
        m_thr = max(q95(sel("mouse", "matched")), q95(sel("mouse", "responsive")))
        rows.append(dict(readout=rd, gtex_thin_clock_median=float(np.nanmedian(clk)),
                         gtex_thin_clock_q10=float(np.nanquantile(clk, 0.1)),
                         gtex_thin_matched_median=float(sel("gtex", "matched").median()), gtex_thin_companion_q95=g_thr,
                         mouse_matched_median=float(sel("mouse", "matched").median()),
                         mouse_matched_q95=q95(sel("mouse", "matched")), mouse_responsive_q95=q95(sel("mouse", "responsive")),
                         mouse_companion_q95=m_thr, P_a=float(np.mean(clk > g_thr)), P_b=float(np.mean(clk > m_thr))))
    P = pd.DataFrame(rows)
    P["min_P"] = P[["P_a", "P_b"]].min(axis=1)
    P["gate"] = bool(P[P.readout == "core"].min_P.iloc[0] >= 0.8)
    P["failed_jobs"] = int(C.leverage.isna().sum())
    P.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(P.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)")


def score():
    limit_threads(1)
    P = pd.read_csv(OUT / "power.csv")
    if not bool(P.gate.iloc[0]):
        print("gate failed in the power stage: not scored (registered)")
        return
    pb, per, _ = M.mouse_table()
    pb = pb[pb.compartment.isin(per)].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    S.update(mouse=pb[M.KEYS + list(WF.genes)].copy(), per=per, arm=pb.groupby("donor").stratum.first())
    rows = []
    for rd, w in (("core", CORE), ("14 genes", WF)):
        st = mouse_states(w)
        lev = leverage(st, S["arm"])
        rng = np.random.default_rng([400, 3])
        null = np.array([leverage(st, S["arm"], rng) for _ in range(N_PERM)])
        p = P[P.readout == rd].iloc[0]
        pv = float((1 + np.sum(null >= lev)) / (1 + N_PERM))
        rows.append(dict(test="H400.1" if rd == "core" else "14 genes (reported)", leverage=lev, p_perm=pv,
                         null_q95=q95(null), mouse_matched_q95=p.mouse_matched_q95,
                         mouse_responsive_q95=p.mouse_responsive_q95,
                         passes=bool(lev > 0 and pv < 0.05 and lev > p.mouse_companion_q95)))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tests.csv", index=False)
    print(T.round(4).to_string(index=False))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "check", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    score() if a.stage == "score" else power(a.jobs, a.stage == "check")


if __name__ == "__main__":
    main()
