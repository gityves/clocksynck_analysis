"""A385b (sheep person loudness, redesigned; gate calibrated on GTEx thinned to the sheep design), registered 2026-11-26.

Statistic (A385 / A331 leverage): per halving, the half-host radii of the animals in both halves, residualised on two
phase harmonics of both phases and on log tissue count per half; correlation; median over 20 halvings.
Designs: D1 (A385: within project, projects with >= 30 animals, Fisher-combined, weights n - 3); D2_8 / D2_12 (animals
with >= 8 / >= 12 tissues, pooled after residualising also on project indicators).

Stages:
  power  GTEx fast deaths thinned to the sheep design (clock 40 replicates; A336 companions, 100 + 100, one replicate
         each), the unthinned A336 clock leverage (reproduction check), sheep companions (the A385 sets, 200 + 200)
         under every design; powers P_a, P_b, design choice and gate. No sheep clock value.
  check  2 sets per arm and 2 replicates; nothing written.
  score  H385b.1 (core) and the 14-gene readout under the chosen design, only if the gate passed.

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A385_sheep_person/scripts/a385b.py --stage power --jobs 6
Outputs (`../output/`): b_sheep_companions.csv, b_gtex_thin.csv, b_power.csv; score: b_tests.csv. Aggregates only
(GTEx: leverages only).
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
sys.path.insert(0, str(PAPER / "exploratory" / "A336_two_readouts" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_person_gtex as G  # noqa: E402
import a385_a386 as M  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
WF, CORE, CORE_IDX = M.WF, M.CORE, M.CORE_IDX
DESIGNS = ("D1", "D2_8", "D2_12")
MIN_K = {"D2_8": 8, "D2_12": 12}
N_REP, N_GTEX_SETS, N_SHEEP_SETS, N_SPLITS, N_PERM = 40, 100, 200, 20, 1000
S = {}


def core_of(w, name):
    return Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, name)


def halves(others, seed):
    """The halvings of `elasticity.split_hosts` (same rng use), as tissue lists."""
    rng = np.random.default_rng(seed)
    others, out = list(others), []
    for _ in range(N_SPLITS):
        p = list(rng.permutation(others))
        out.append((sorted(p[: len(p) // 2]), sorted(p[len(p) // 2:])))
    return out


def states(ds, halvings):
    """Per halving: animal ids, complex half-host states and tissue counts per half."""
    d = ds.data
    tis = pd.DataFrame({"donor": np.asarray(d.resection).astype(str), "comp": np.asarray(d.compartment)})
    out = []
    for sa, sb in halvings:
        hA, hB = ds.host(sa), ds.host(sb)
        hA.index, hB.index = hA.index.astype(str), hB.index.astype(str)
        j = hA.index.intersection(hB.index)
        cnt = [tis[tis.comp.isin(s)].groupby("donor").size().reindex(j).to_numpy(float) for s in (sa, sb)]
        out.append((np.asarray(j), (hA.x + 1j * hA.y)[j].to_numpy(), (hB.x + 1j * hB.y)[j].to_numpy(),
                    np.log(np.column_stack(cnt))))
    return out


def pooled_leverage(a, b, cov, proj):
    ra, rb = np.abs(a), np.abs(b)
    P = pd.get_dummies(pd.Series(proj), drop_first=True).to_numpy(float)
    H = np.column_stack([np.ones(len(a)), P, cov] + [g(k * p) for p in (np.angle(a), np.angle(b)) for k in (1, 2)
                                                     for g in (np.cos, np.sin)])
    res = lambda v: v - H @ np.linalg.lstsq(H, v, rcond=None)[0]  # noqa: E731
    return float(np.corrcoef(res(ra), res(rb))[0, 1]), len(a) - H.shape[1] - 1


def design_stats(st, project, ntis, rng=None):
    """design -> (median leverage over halvings, median weight). rng: permute B within project (null)."""
    out = {}
    for dsg in DESIGNS:
        L, N = [], []
        for j, a, b, cov in st:
            pr = project.reindex(j).to_numpy()
            if rng is not None:
                b = b.copy()
                for p in np.unique(pr):
                    m = np.flatnonzero(pr == p)
                    b[m] = b[rng.permutation(m)]
            if dsg == "D1":
                z = w = 0.0
                for p in np.unique(pr):
                    m = np.flatnonzero(pr == p)
                    if len(m) >= M.MIN_PROJECT:
                        z += (len(m) - 3) * np.arctanh(M.leverage(a[m], b[m], cov[m]))
                        w += len(m) - 3
                L.append(np.tanh(z / w) if w > 0 else np.nan)
                N.append(w)
            else:
                m = np.flatnonzero(ntis.reindex(j).to_numpy() >= MIN_K[dsg])
                lv, df = pooled_leverage(a[m], b[m], cov[m], pr[m])
                L.append(lv)
                N.append(df)
        out[dsg] = (float(np.nanmedian(L)), float(np.median(N)))
    return out


# ---- sheep ---------------------------------------------------------------------------------------------------

def sheep_states(w):
    ds = cs.Dataset(S["sheep"][M.META + list(w.genes)], anchors=S["per"], waveform=w, mixed_strata=True,
                    name=w.name, min_donors=5)
    return states(ds, M.halvings(S["per"]))


def sheep_sets(n):
    lm = M.log_mean(S["sheep_full"])
    sets = {kind: pseudo_waveforms(lm, n, np.random.default_rng([385, tag]), WF, kind=kind)
            for kind, tag in (("matched", 2), ("responsive", 3))}
    return [(kind, j, rd, w if rd == "14 genes" else core_of(w, f"core{j}"))
            for kind, ws in sets.items() for j, w in enumerate(ws) for rd in ("14 genes", "core")]


# ---- GTEx thinned --------------------------------------------------------------------------------------------

def assign(r):
    """Replicate r: sheep animal -> (GTEx donor, tissues kept)."""
    rng = np.random.default_rng([3852, r])
    animals = S["animals"].sort_values("k", ascending=False)
    tis = S["gtex_tissues"]
    free = sorted(tis.index)  # a list: deterministic order (set order depends on hash randomisation)
    rows = []
    for aid, x in animals.iterrows():
        cand = [d for d in free if len(tis[d]) >= x.k]
        if not cand:
            top = max(len(tis[d]) for d in free)
            cand = [d for d in free if len(tis[d]) == top]
        d = cand[rng.integers(len(cand))]
        free.remove(d)
        keep = rng.choice(tis[d], min(int(x.k), len(tis[d])), replace=False)
        rows += [(aid, d, t) for t in keep]
    return pd.DataFrame(rows, columns=["animal", "gtex", "compartment"])


def gtex_thin(w, r):
    a = assign(r)
    g = S["gtex"].merge(a, left_on=["donor", "compartment"], right_on=["gtex", "compartment"])
    g = g.drop(columns=["donor", "gtex"]).rename(columns={"animal": "donor"})
    ds = cs.Dataset(g[G.META + list(w.genes)], anchors=S["comps"], waveform=w, name=w.name)
    st = states(ds, halves(S["comps"], [331, 1]))
    ntis = g.groupby("donor").compartment.nunique()
    return st, ntis


def gtex_sets(n):
    lm = S["gtex_lm"]
    sets = {kind: pseudo_waveforms(lm, G.N_SETS, np.random.default_rng([336, 2, tag]), WF, kind=kind)[:n]
            for kind, tag in (("matched", 1), ("responsive", 2))}
    return [(kind, j, rd, w if rd == "14 genes" else core_of(w, f"core{j}"))
            for kind, ws in sets.items() for j, w in enumerate(ws) for rd in ("14 genes", "core")]


# ---- jobs ----------------------------------------------------------------------------------------------------

def task(job):
    warnings.filterwarnings("ignore")
    arm, kind, i, rd, w = job
    base = dict(arm=arm, kind=kind, set=i, readout=rd)
    try:
        if arm == "sheep":
            st = sheep_states(w)
            res = design_stats(st, S["project"], S["ntis_sheep"])
        elif arm == "gtex_full":
            ds = cs.Dataset(S["gtex"][G.META + list(w.genes)], anchors=S["comps"], waveform=w, name=w.name)
            st = states(ds, halves(S["comps"], [331, 1]))
            lev = float(np.median([E.radius_covariance(a, b)["leverage"] for _, a, b, _ in st]))
            return [dict(base, design="A336", leverage=lev, weight=np.nan)]
        else:
            r = i if kind == "clock" else 1000 + i + (0 if kind == "matched" else N_GTEX_SETS)
            st, ntis = gtex_thin(w, r)
            res = design_stats(st, S["project"], ntis)
    except Exception as e:  # counted, not dropped silently
        return [dict(base, design=None, error=f"{type(e).__name__}: {e}"[:200])]
    return [dict(base, design=k, leverage=v[0], weight=v[1]) for k, v in res.items()]


def setup(n_sheep, n_gtex, n_rep):
    t0 = time.time()
    sheep, per, md = M.load()
    S.update(sheep_full=sheep, per=per)
    project = sheep.groupby("donor").stratum.first()
    ntis = sheep.groupby("donor").compartment.nunique()
    big = project.map(project.value_counts()) >= M.MIN_PROJECT
    used = big | (ntis >= min(MIN_K.values()))
    S.update(project=project, ntis_sheep=ntis, animals=pd.DataFrame({"k": ntis[used]}))
    js = sheep_sets(n_sheep)
    genes = sorted({g for *_, w in js for g in w.genes} | set(WF.genes))
    S["sheep"] = sheep[M.META + genes].copy()
    del S["sheep_full"]
    gt, _ = G.load()
    gl = [c for c in gt.columns if c not in G.META]
    x = gt[gl].to_numpy(np.float32)
    lm = pd.Series(np.log((x / gt.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=gl)[x.mean(0) >= 1]
    del x
    S["gtex_lm"] = lm[~lm.index.duplicated()]
    jg = gtex_sets(n_gtex)
    genes = sorted({g for *_, w in jg for g in w.genes} | set(WF.genes))
    gt["donor"] = gt.donor.astype(str)
    S["gtex"] = gt[G.META + genes].copy()
    del gt
    S["comps"] = sorted(S["gtex"].compartment.unique())
    S["gtex_tissues"] = S["gtex"].groupby("donor").compartment.apply(lambda v: sorted(v.unique()))
    print(f"sheep {project.size} animals, {len(per)} tissues, {int(used.sum())} animals mimicked; GTEx "
          f"{S['gtex'].donor.nunique()} donors, {len(S['comps'])} tissues ({time.time() - t0:.0f}s)", flush=True)
    jobs = [("gtex_full", "clock", 0, rd, w) for rd, w in (("14 genes", WF), ("core", CORE))]
    jobs += [("gtex", "clock", r, rd, w) for r in range(n_rep) for rd, w in (("14 genes", WF), ("core", CORE))]
    jobs += [("gtex", k, j, rd, w) for k, j, rd, w in jg]
    jobs += [("sheep", k, j, rd, w) for k, j, rd, w in js]
    return jobs


def run(jobs, n_jobs):
    rows, t0 = [], time.time()
    with Pool(n_jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows += r
            if (k + 1) % 25 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def q95(v):
    return float(np.nanquantile(v, 0.95))


def power(n_jobs, check):
    limit_threads(1)
    t0 = time.time()
    jobs = setup(2 if check else N_SHEEP_SETS, 2 if check else N_GTEX_SETS, 2 if check else N_REP)
    R = run(jobs, n_jobs)
    if check:
        print(R.groupby(["arm", "kind", "readout"]).agg(rows=("set", "size"), err=("design", lambda v: v.isna().sum()))
              .to_string(), f"({time.time() - t0:.0f}s)")
        return
    R[R.arm == "sheep"].to_csv(OUT / "b_sheep_companions.csv", index=False)
    R[R.arm != "sheep"].to_csv(OUT / "b_gtex_thin.csv", index=False)
    stored = pd.read_csv(OUT / "companions.csv")
    rows = []
    for rd in ("14 genes", "core"):
        full = R[(R.arm == "gtex_full") & (R.readout == rd)].leverage
        for dsg in DESIGNS:
            sel = lambda arm, kind: R[(R.arm == arm) & (R.kind == kind) & (R.readout == rd) & (R.design == dsg)]  # noqa: E731
            clk = sel("gtex", "clock").leverage.to_numpy(float)
            g_thr = max(q95(sel("gtex", "matched").leverage), q95(sel("gtex", "responsive").leverage))
            s_m, s_r = sel("sheep", "matched"), sel("sheep", "responsive")
            s_thr = max(q95(s_m.leverage), q95(s_r.leverage))
            row = dict(readout=rd, design=dsg, gtex_unthinned=float(full.iloc[0]) if len(full) else np.nan,
                       gtex_thin_clock_median=float(np.nanmedian(clk)), gtex_thin_clock_q10=float(np.nanquantile(clk, 0.1)),
                       gtex_thin_matched_median=float(sel("gtex", "matched").leverage.median()),
                       gtex_thin_companion_q95=g_thr, sheep_matched_median=float(s_m.leverage.median()),
                       sheep_matched_q95=q95(s_m.leverage), sheep_responsive_q95=q95(s_r.leverage),
                       sheep_companion_q95=s_thr, weight=float(sel("sheep", "matched").weight.median()),
                       P_a=float(np.mean(clk > g_thr)), P_b=float(np.mean(clk > s_thr)))
            if dsg == "D1":
                st = stored[stored.label == rd].set_index(["kind", "set"]).leverage
                mine = s_m.set_index(["kind", "set"]).leverage
                row["d1_reproduces_a385"] = bool(np.allclose(mine.reindex(st.loc[["matched"]].index).to_numpy(float),
                                                             st.loc[["matched"]].to_numpy(float), atol=1e-6))
            rows.append(row)
    P = pd.DataFrame(rows)
    P["min_P"] = P[["P_a", "P_b"]].min(axis=1)
    c = P[P.readout == "core"].copy()
    c["order"] = c.design.map({d: i for i, d in enumerate(DESIGNS)})
    best = c.sort_values(["min_P", "order"], ascending=[False, True]).iloc[0]
    P["chosen"] = P.design == best.design
    P["gate"] = bool(best.min_P >= 0.8)
    P["failed_jobs"] = int(R.design.isna().sum())
    P.to_csv(OUT / "b_power.csv", index=False)
    pd.set_option("display.width", 250)
    print(P.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)")


def score(n_jobs, deviation=False):
    limit_threads(1)
    t0 = time.time()
    P = pd.read_csv(OUT / "b_power.csv")
    if not bool(P.gate.iloc[0]) and not deviation:
        print("gate failed in the power stage: not scored (registered)")
        return
    dsg = P[P.chosen].design.iloc[0]
    sheep, per, _ = M.load()
    S.update(sheep=sheep[M.META + list(WF.genes)].copy(), per=per,
             project=sheep.groupby("donor").stratum.first(), ntis_sheep=sheep.groupby("donor").compartment.nunique())
    rows = []
    for rd, w in (("core", CORE), ("14 genes", WF)):
        st = sheep_states(w)
        lev = design_stats(st, S["project"], S["ntis_sheep"])[dsg][0]
        rng = np.random.default_rng([385, 4])
        null = np.array([design_stats(st, S["project"], S["ntis_sheep"], rng)[dsg][0] for _ in range(N_PERM)])
        p = P[(P.readout == rd) & (P.design == dsg)].iloc[0]
        pv = float((1 + np.sum(null >= lev)) / (1 + N_PERM))
        rows.append(dict(test="H385b.1" if rd == "core" else "14 genes (reported)", design=dsg, leverage=lev,
                         p_perm=pv, null_q95=q95(null), sheep_matched_q95=p.sheep_matched_q95,
                         sheep_responsive_q95=p.sheep_responsive_q95,
                         passes=bool(lev > 0 and pv < 0.05 and lev > p.sheep_companion_q95)))
        print(rd, "done", f"({time.time() - t0:.0f}s)", flush=True)
    T = pd.DataFrame(rows)
    T["deviation"] = bool(deviation)
    T.to_csv(OUT / ("b_tests_deviation.csv" if deviation else "b_tests.csv"), index=False)
    print(T.round(4).to_string(index=False))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "check", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--deviation", action="store_true",
                    help="score despite the failed gate (user decision 2026-11-26; registered deviation)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    score(a.jobs, a.deviation) if a.stage == "score" else power(a.jobs, a.stage == "check")


if __name__ == "__main__":
    main()
