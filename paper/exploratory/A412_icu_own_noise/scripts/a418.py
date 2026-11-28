"""A418 (registered in `paper/docs/registered_experiments.md`): is the fixed own noise of A412 specific to the clock, or
shared by rhythmic output programmes of the same organs?

A412 exactly (GTEx fast deaths vs ventilated with a living brain, 12 body tissues, matched design, host anchored on the A164
tissues, A412 estimator `dpsi_k = log psi_fast - log psi_vent`, `dS_k = log S_fast - log S_vent`), with 200 sets of 14
rhythmic output genes instead of non-rhythmic companions.
Pool: non-clock genes rhythmic by time of death in GTEx fast deaths (cosinor of the log fraction on 24 h time of death per
body tissue, Fisher-combined over tissues, BH q < 0.05; mean count >= 1 over all rows as A263; clock-like and A263
responsive genes excluded). Sets: each clock gene replaced by a pool gene of the same pooled peak-phase bin (6 x 4 h) and
expression tercile (median log fraction), seed [418, 0]; each set's waveform is fitted on the fast deaths by time of death
(`clocksync.waveform.fit_waveform`, group = tissue) and read in both groups as the clock.

  --stage power  pool, sets, per set ONLY median dS (shared-rhythm change) and whether kappa is finite (dpsi is computed
                 inside the worker to test finiteness and discarded: never written or printed). Gate: >= 100 finite sets and
                 median over sets of median dS >= 0.2.
  --stage score  kappa_j = median_k dpsi_jk / median_k dS_jk (14 genes, matched design) per set; clock raw kappa recomputed
                 in the same run (and the A412 companion-excess kappa beside it); H418.1 = clock raw kappa < q05 of sets.
                 Descriptive: set band, by phase bin (genes of that bin only, >= 3), core loop.
Aggregates only (GTEx time of death and donors are controlled access).

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A412_icu_own_noise/scripts/a418.py --stage power|score [--jobs 4]
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import a412 as K  # noqa: E402  (A412 estimator, design, bootstrap)
import clocksync as cs  # noqa: E402
from a349_common import gtex  # noqa: E402
from a349_common import H  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402
from clocksync._core.plane import FLOOR  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.waveform import Waveform, fit_waveform  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402

OUT = HERE.parents[1] / "output" / "a418"
SEED = 418
WF = K.WF
TRIO = K.TRIO
N_SETS = 200
CORE_POS = [i for i, g in enumerate(WF.genes) if g not in TRIO]
CLOCK_RELATED = ["ARNTL2", "CLOCK", "BHLHE40", "BHLHE41", "NFIL3", "HLF", "CRY1", "CRY2", "PER1", "PER2", "PER3",
                 "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "DBP", "TEF", "CIART", "CSNK1D", "CSNK1E", "FBXL3", "FBXL21"]
S = {}


def bh(p):
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    o = np.argsort(p[ok])
    r = p[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1)
    q[np.flatnonzero(ok)[o]] = np.minimum(np.minimum.accumulate(r[::-1])[::-1], 1)
    return q


# ---- data, pool, sets ------------------------------------------------------------------------------------------

def load():
    """Genome-wide pseudobulk of the A412 design (same donors and tissues as `a412.setup`), the A412 designs, fast-death
    time of death. Returns pb (all genes), tod (donor -> h), anc, designs, fast ids."""
    T, fast, pm = K.C.tissues()
    ids_f = T.donors[fast].astype(str)
    # designs: reuse a412.setup's logic through its pieces (matched / pooled), without its companion gene list
    import a412 as M
    tmp = OUT / "_clock_only_sets.csv"  # a412.setup reads the companion table only for its gene list: clock genes here
    pd.DataFrame(dict(kind=["matched"], genes=[";".join(WF.genes)])).to_csv(tmp, index=False)
    M.SETS = tmp
    _, _, _, designs, _ = M.setup()
    tmp.unlink()
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    allg = list(a.var_names)
    a.file.close()
    pb, donors = gtex.load(tissues=list(K.SH.TISSUES), hardy=(0, 1, 2), genes=allg, depth="library_size", require_tod=True)
    keep = set(designs["pooled"]["fast"]) | set(designs["pooled"]["ventilated"])
    pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    anc = [t for t in K.SH.TISSUES if t in set(pb.compartment)]
    return pb, donors.tod.astype(float) % 24, anc, designs, set(map(str, ids_f))


def cosinor(pb, fast_ids, tod, genes):
    """Per body tissue and gene on fast deaths: cosinor of the log fraction on 24 h time of death -> amp, phase (h), p, median."""
    amp, ph, pv, med = {}, {}, {}, {}
    for k in H.BODY:
        t = pb[(pb.compartment == k) & pb.donor.isin(fast_ids)]
        if len(t) < 10:
            continue
        y = np.log(np.maximum(t[genes].to_numpy(np.float64), FLOOR) / t.depth.to_numpy(float)[:, None])
        h = tod.reindex(t.donor).to_numpy(float) * 2 * np.pi / 24
        X = np.column_stack([np.ones(len(h)), np.cos(h), np.sin(h)])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        rss = ((y - X @ beta) ** 2).sum(0)
        tss = ((y - y.mean(0)) ** 2).sum(0)
        df = len(h) - 3
        with np.errstate(divide="ignore", invalid="ignore"):
            F = ((tss - rss) / 2) / (rss / df)
        pv[k] = np.where(rss > 1e-12, stats.f.sf(F, 2, df), np.nan)
        amp[k] = np.hypot(beta[1], beta[2])
        ph[k] = np.mod(np.arctan2(beta[2], beta[1]), 2 * np.pi) * 24 / (2 * np.pi)
        med[k] = np.median(y, axis=0)
    f = lambda d: pd.DataFrame(d, index=genes)  # noqa: E731
    return f(amp), f(ph), f(pv), f(med)


def pool_table(pb, fast_ids, tod):
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "depth", "stratum")]
    x = pb[genes].to_numpy(np.float32)
    expressed = pd.Series(x.mean(0) >= 1, index=genes)  # A263 expression filter
    del x
    genes = [g for g in genes if expressed[g] or g in WF.genes]
    amp, ph, pv, med = cosinor(pb, fast_ids, tod, genes)
    ok = pv.notna().sum(axis=1)
    fis = -2 * np.log(np.clip(pv.to_numpy(float), 1e-300, 1)).sum(axis=1)
    p_comb = pd.Series(np.where(ok > 0, stats.chi2.sf(fis, 2 * ok.clip(lower=1)), np.nan), index=genes)
    q_comb = pd.Series(bh(p_comb.to_numpy()), index=genes)
    z = (amp * np.exp(1j * ph * 2 * np.pi / 24)).sum(1)
    phase = np.mod(np.angle(z) * 24 / (2 * np.pi), 24)
    T = pd.DataFrame(dict(phase_h=phase, bin=(phase // 4).astype(int) % 6, median_logfrac=med.median(1),
                          n_tissues=ok, p_combined=p_comb, q_combined=q_comb))
    excl = set(WF.genes) | {canonical(g) for g in CLOCK_RELATED} | {canonical(g) for g in CLOCK_LIKE} \
        | {canonical(g) for g in RESPONSIVE_GENES}
    pool = T[(T.q_combined < 0.05) & ~T.index.isin(excl) & (T.index == [canonical(g) for g in T.index])
             & (T.index.isin([g for g in genes if expressed[g]]))].copy()
    cuts = np.quantile(pool.median_logfrac, [1 / 3, 2 / 3])
    pool["tercile"] = np.searchsorted(cuts, pool.median_logfrac)
    ck = T.loc[list(WF.genes)].copy()
    ck["tercile"] = np.searchsorted(cuts, ck.median_logfrac)
    merged = [int(b) for b in range(6) if (pool[pool.bin == b].groupby("tercile").size().reindex(range(3), fill_value=0) < 10).any()]
    for b in merged:  # fewer than 10 pool genes in a (bin, tercile) cell: terciles of that bin merged
        pool.loc[pool.bin == b, "tercile"] = 0
        ck.loc[ck.bin == b, "tercile"] = 0
    return pool, ck, dict(n_combined_q05=int((T.q_combined < 0.05).sum()), n_tested=len(genes), merged_bins=merged)


def draw_sets(pool, ck, n=N_SETS):
    rng = np.random.default_rng([SEED, 0])
    cells = {key: list(g.index) for key, g in pool.groupby(["bin", "tercile"])}
    rows = []
    for i in range(n):
        used = set()
        for pos, g in enumerate(WF.genes):
            avail = [x for x in cells.get((ck.bin[g], ck.tercile[g]), []) if x not in used]
            if not avail:
                raise RuntimeError(f"empty cell for {g}: bin {ck.bin[g]}, tercile {ck.tercile[g]}")
            c = avail[rng.integers(len(avail))]
            used.add(c)
            rows.append(dict(set=i, position=pos, clock_gene=g, gene=c, bin=int(ck.bin[g]), core=pos in CORE_POS))
    return pd.DataFrame(rows)


def set_waveform(genes):
    """Waveform of a set fitted on fast deaths by time of death (12 body tissues, group = tissue)."""
    r = S["fit_rows"]
    y = np.log(np.maximum(r[genes].to_numpy(np.float64), FLOOR) / r.depth.to_numpy(float)[:, None])
    return fit_waveform(y, S["tod"].reindex(r.donor).to_numpy(float), genes, group=r.compartment.to_numpy())


# ---- per set reads ---------------------------------------------------------------------------------------------

def read(genes, wf, drop=None, name="a418"):
    """dtable (dpsi, dS per tissue) of a gene set with waveform wf on the matched design."""
    cols = ["donor", "compartment", "depth"] + list(genes)
    ds = cs.Dataset(S["pb"][cols], anchors=S["anc"], waveform=wf, name=name, drop_genes=drop)
    d, Hs = K.restrict(ds, S["anc"], S["host"])
    return K.dtable(K.estimates(d, Hs, S["design"]))


def finite_kappa(t):
    mp, ms = np.nanmedian(t.dpsi), np.nanmedian(t.dS)
    return bool(np.isfinite(mp) and np.isfinite(ms) and ms != 0), float(ms)


def task_power(i):
    t0 = time.time()
    genes = list(S["sets"][S["sets"].set == i].sort_values("position").gene)
    try:
        t = read(genes, set_waveform(genes), name=f"a418 set {i}")
        fin, ms = finite_kappa(t)  # dpsi enters only the finiteness flag
        return dict(set=i, finite=fin, median_dS=ms, n_tissues=int(t.dS.notna().sum()), seconds=time.time() - t0)
    except Exception as e:
        return dict(set=i, finite=False, median_dS=np.nan, error=f"{type(e).__name__}: {e}"[:200], seconds=time.time() - t0)


def kappa_of(t):
    mp, ms = np.nanmedian(t.dpsi), np.nanmedian(t.dS)
    return float(mp / ms) if np.isfinite(mp) and np.isfinite(ms) and ms != 0 else np.nan


def task_score(i):
    genes = list(S["sets"][S["sets"].set == i].sort_values("position").gene)
    bins = S["sets"][S["sets"].set == i].sort_values("position").bin.to_numpy()
    out = dict(set=i)
    try:
        wf = set_waveform(genes)
        t = read(genes, wf, name=f"a418 set {i}")
        out.update(kappa=kappa_of(t), median_dpsi=float(np.nanmedian(t.dpsi)), median_dS=float(np.nanmedian(t.dS)))
        core = [genes[p] for p in CORE_POS]
        out["kappa_core"] = kappa_of(read(core, wf.subset(core), name=f"a418 set {i} core"))
        for b in range(6):
            g = [genes[p] for p in range(len(genes)) if bins[p] == b]
            out[f"kappa_bin{b}"] = kappa_of(read(g, wf.subset(g), name=f"a418 set {i} bin{b}")) if len(g) >= 3 else np.nan
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"[:200]
    return out


def run_pool(fn, n, jobs, t0, tag):
    rows = []
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(fn, range(n))):
            rows.append(r)
            if (k + 1) % 5 == 0 or k + 1 == n:
                print(f"{tag} {k + 1}/{n} sets {time.time() - t0:.0f}s", flush=True)
    return pd.DataFrame(rows).sort_values("set").reset_index(drop=True)


def prepare(a, t0):
    pb, tod, anc, designs, fast_ids = load()
    print(f"data {len(pb)} rows x {pb.shape[1]} cols ({time.time() - t0:.0f}s)", flush=True)
    pool, ck, info = pool_table(pb, fast_ids, tod)
    sets = draw_sets(pool, ck)
    need = sorted(set(sets.gene) | set(WF.genes))
    pb = pb[["donor", "compartment", "depth"] + need].copy()
    m = designs["matched"]
    S.update(pb=pb, tod=tod, anc=anc, sets=sets, design=m, host=np.union1d(m["fast"], m["ventilated"]),
             fit_rows=pb[pb.compartment.isin(H.BODY) & pb.donor.isin(fast_ids)])
    return pool, ck, sets, info


def power(a):
    t0 = time.time()
    pool, ck, sets, info = prepare(a, t0)
    pool.reset_index(names="gene").to_csv(OUT / "pool.csv", index=False)
    sets.to_csv(OUT / "sets.csv", index=False)
    print(f"pool {len(pool)} genes ({info}); per bin {pool.bin.value_counts().sort_index().to_dict()}; "
          f"clock bins {dict(ck.bin)} terciles {dict(ck.tercile)}", flush=True)
    n = 3 if a.check else N_SETS
    R = run_pool(task_power, n, a.jobs, t0, "power")
    ok = R[R.finite]
    P = pd.DataFrame([dict(pool_size=len(pool), sets=len(R), sets_finite=len(ok), median_set_dS=float(ok.median_dS.median()),
                           q05_dS=float(ok.median_dS.quantile(0.05)), q95_dS=float(ok.median_dS.quantile(0.95)),
                           errors=int(R.get("error", pd.Series(dtype=object)).notna().sum()),
                           gate_finite=bool(len(ok) >= 100), gate_dS=bool(ok.median_dS.median() >= 0.2),
                           runtime_s=round(time.time() - t0))])
    P["gate_pass"] = P.gate_finite & P.gate_dS
    P.to_csv(OUT / ("power_check.csv" if a.check else "power.csv"), index=False)
    print(P.T.to_string(), flush=True)
    if "error" in R:
        print("errors:", R.error.dropna().unique()[:5].tolist())


def score(a):
    t0 = time.time()
    P = pd.read_csv(OUT / "power.csv").iloc[0]
    if not bool(P.gate_pass):
        sys.exit("gate not passed: score stage not run")
    pool, ck, sets, info = prepare(a, t0)
    chk = pd.read_csv(OUT / "sets.csv")
    assert chk.gene.tolist() == sets.gene.tolist(), "sets differ from the power stage"
    R = run_pool(task_score, 3 if a.check else N_SETS, a.jobs, t0, "score")
    R.to_csv(OUT / ("sets_kappa_check.csv" if a.check else "sets_kappa.csv"), index=False)
    # clock, same run
    wfc = WF
    clock = {}
    for lab, drop in (("14 genes", None), ("core", TRIO)):
        t = read(list(WF.genes), wfc, drop=drop, name=f"a418 clock {lab}")
        clock[lab] = dict(raw=kappa_of(t), median_dpsi=float(np.nanmedian(t.dpsi)), median_dS=float(np.nanmedian(t.dS)))
    t14 = read(list(WF.genes), wfc, name="a418 clock excess")
    ct = pd.read_csv(K.OUT / "companions.csv")
    cm = K.comp_matrix(ct, "matched", list(t14.index))
    excess = K.kappa(t14.dpsi, t14.dS, np.nanmedian(cm.to_numpy(float), axis=0))
    ok = R.dropna(subset=["kappa"])
    k = ok.kappa
    q05, q50, q95 = k.quantile([0.05, 0.5, 0.95])
    ck14 = clock["14 genes"]["raw"]
    rows = [dict(row="set kappa 14 genes", n=len(k), q05=q05, median=q50, q95=q95, sd=float(k.std(ddof=1)),
                 mean=float(k.mean()), sets_below_zero=float((k < 0).mean())),
            dict(row="clock raw kappa 14 genes", value=ck14, median_dpsi=clock["14 genes"]["median_dpsi"],
                 median_dS=clock["14 genes"]["median_dS"], pct_sets_at_or_below_clock=100 * float((k <= ck14).mean())),
            dict(row="clock excess kappa (A412 companions)", value=excess),
            dict(row="H418.1 clock raw kappa < q05 of sets", value=float(ck14 < q05)),
            dict(row="clock raw kappa core loop", value=clock["core"]["raw"], median_dpsi=clock["core"]["median_dpsi"],
                 median_dS=clock["core"]["median_dS"],
                 pct_sets_at_or_below_clock=100 * float((ok.kappa_core.dropna() <= clock["core"]["raw"]).mean()))]
    kc = ok.kappa_core.dropna()
    rows.append(dict(row="set kappa core loop", n=len(kc), q05=kc.quantile(0.05), median=kc.median(), q95=kc.quantile(0.95)))
    for b in range(6):
        g = [x for x, bb in zip(WF.genes, ck.bin[list(WF.genes)]) if bb == b]
        kb = ok[f"kappa_bin{b}"].dropna()
        cb = np.nan
        if len(g) >= 3:
            cb = kappa_of(read(g, WF.subset(g), name=f"a418 clock bin{b}"))
        rows.append(dict(row=f"phase bin {b} ({4 * b}-{4 * b + 4} h)", n=len(kb), clock_genes_in_bin=len(g),
                         q05=kb.quantile(0.05) if len(kb) else np.nan, median=kb.median() if len(kb) else np.nan,
                         q95=kb.quantile(0.95) if len(kb) else np.nan, value=cb))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / ("tests_check.csv" if a.check else "tests.csv"), index=False)
    print(T.round(3).to_string(index=False))
    print(f"errors {int(R.get('error', pd.Series(dtype=object)).notna().sum())}, done {time.time() - t0:.0f}s", flush=True)


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true", help="3 sets")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    limit_threads(4)
    (power if a.stage == "power" else score)(a)


if __name__ == "__main__":
    main()
