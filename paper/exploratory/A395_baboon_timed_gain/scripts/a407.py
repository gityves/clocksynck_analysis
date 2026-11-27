"""A407 (is the timed own-noise law of A395 specific to the clock? rhythmic non-clock companions), registered 2026-11-27.

A395 pipeline (baboon, Mure 2018, peripheral tissues, ZT known; `timed_loadings("zt")`, precondition by 1,000 ZT
permutations, animal correction, b_t = OLS slope of log psi on log |a|^2, phase-scatter slope) applied to 200 gene sets
that replace each clock gene by a rhythmic pool gene of the same pooled-phase bin (six 4 h bins) and median-expression
tercile. Pool: genome-wide baboon table of A263 (`build_genomewide.py`); per peripheral tissue a cosinor of the log
fraction on ZT, Fisher-combined (Note 1) BH q < 0.05, median log fraction >= the lowest median of the 14 clock genes in
that tissue in >= half of the tissues; the 14 clock genes and the clock-related list are excluded. Each set's waveform is
fitted on the baboon (`clocksync.waveform.fit_waveform`, hours = ZT, group = tissue).
(Note 1: pool by Fisher-combined per-tissue cosinor p over the 35 tissues, BH over genes q < 0.05; terciles merged in a bin
with a cell < 10 genes. `FLOOR` is imported from `clocksync._core.plane` (private; no public export yet).)

Stages:
  power  pool, sets, and per set: passing organs (14 positions), b_t, scatter slope (14 positions and the 11 core
         positions, descriptive). The b_t and slopes are written to companions_unread.csv and not printed; only the
         number of sets with >= 15 passing organs is reported. The clock under a baboon-fitted waveform is not computed.
  score  clock re-read with the baboon-fitted 14-gene waveform (and core), H407.1 against the q05 of the companion b_t,
         companion summaries, the A395 GTEx-waveform value beside it.

    paper/common/run.sh paper/exploratory/A395_baboon_timed_gain/scripts/a407.py --stage power --jobs 4
Outputs (`../output/a407/`): pool.csv, sets.csv, companions_unread.csv, power.csv; tests.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import a395 as B  # noqa: E402  (A395 pipeline: timed, organ_table, readouts, slope, boot)
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402
from clocksync._core.plane import FLOOR  # noqa: E402

A = B.A
OUT = HERE.parents[1] / "output" / "a407"
SEED = 407
KEYS = B.KEYS
WF = A.WF
TRIO = B.TRIO
CORE_POS = [i for i, g in enumerate(WF.genes) if g not in TRIO]
N_SETS, N_PERM, MIN_ORGANS = 200, B.N_PERM, B.MIN_ORGANS
CLOCK_RELATED = ["ARNTL2", "CLOCK", "BHLHE40", "BHLHE41", "NFIL3", "HLF", "CRY1", "CRY2", "PER1", "PER2", "PER3",
                 "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "DBP", "TEF", "CIART", "CSNK1D", "CSNK1E", "FBXL3", "FBXL21"]
S = {}


# ---- pool ------------------------------------------------------------------------------------------------------

def load_genome():
    """All rows of the peripheral tissues with every gene (counts, depth), the peripheral tissue list, ZT of each row."""
    pb = A.table("baboon")
    pb["donor"] = pb.donor.astype(str)
    _, per, _ = A.groups("baboon", sorted(pb.compartment.unique()))
    _, donors = B.data.load("baboon")
    zt = donors.zt.astype(float)
    zt.index = zt.index.astype(str)
    pb = pb[pb.compartment.isin(per)].reset_index(drop=True)
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "depth", "stratum")]
    return pb, per, zt, genes


def bh(p):
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    o = np.argsort(p[ok])
    r = p[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1)
    q[np.flatnonzero(ok)[o]] = np.minimum(np.minimum.accumulate(r[::-1])[::-1], 1)
    return q


def cosinor(pb, per, zt, genes):
    """Per tissue and gene: cosinor of the log fraction on ZT (24 h): amplitude, acrophase (h), BH q, median log fraction."""
    amp, ph, q, med, pv = {}, {}, {}, {}, {}
    for k in per:
        t = pb[pb.compartment == k]
        x = t[genes].to_numpy(np.float64)
        y = np.log(np.maximum(x, FLOOR) / t.depth.to_numpy(float)[:, None])
        h = zt.reindex(t.donor).to_numpy(float) * 2 * np.pi / 24
        X = np.column_stack([np.ones(len(h)), np.cos(h), np.sin(h)])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        rss = ((y - X @ beta) ** 2).sum(0)
        tss = ((y - y.mean(0)) ** 2).sum(0)
        df = len(h) - 3
        with np.errstate(divide="ignore", invalid="ignore"):
            F = ((tss - rss) / 2) / (rss / df)
        p = np.where(rss > 1e-12, stats.f.sf(F, 2, df), np.nan)
        pv[k] = p
        amp[k] = np.hypot(beta[1], beta[2])
        ph[k] = np.mod(np.arctan2(beta[2], beta[1]), 2 * np.pi) * 24 / (2 * np.pi)
        q[k], med[k] = bh(p), np.median(y, axis=0)
    f = lambda d: pd.DataFrame(d, index=genes)
    return f(amp), f(ph), f(q), f(med), f(pv)


def pool_table(pb, per, zt, genes):
    amp, ph, q, med, pv = cosinor(pb, per, zt, genes)
    clock = list(WF.genes)
    lowest = med.loc[clock].min(0)  # per tissue: lowest median log fraction of the 14 clock genes
    n_rhy = (q < 0.05).sum(axis=1)  # tissues with per-tissue q < 0.05 (descriptive; the pool rule is Note 1)
    ok = pv.notna().sum(axis=1)
    fis = -2 * np.log(np.clip(pv.to_numpy(float), 1e-300, 1)).sum(axis=1)
    p_comb = pd.Series(np.where(ok > 0, stats.chi2.sf(fis, 2 * ok.clip(lower=1)), np.nan), index=pv.index)
    q_comb = pd.Series(bh(p_comb.to_numpy()), index=pv.index)
    n_exp = (med >= lowest).sum(axis=1)
    z = (amp * np.exp(1j * ph * 2 * np.pi / 24)).sum(1)  # amplitude-weighted circular mean of the tissue phases
    phase = np.mod(np.angle(z) * 24 / (2 * np.pi), 24)
    pooled_med = med.median(1)
    T = pd.DataFrame(dict(phase_h=phase, bin=(phase // 4).astype(int) % 6, median_logfrac=pooled_med,
                          n_rhythmic=n_rhy, n_expressed=n_exp, p_combined=p_comb,
                          q_combined=q_comb))
    excl = set(clock) | {canonical(g) for g in CLOCK_RELATED}
    n_exp_min = int(np.ceil(0.5 * len(per)))
    pool = T[(T.q_combined < 0.05) & (n_exp >= n_exp_min) & ~T.index.isin(excl)
             & (T.index == [canonical(g) for g in T.index])].copy()
    cuts = np.quantile(pool.median_logfrac, [1 / 3, 2 / 3])
    pool["tercile"] = np.searchsorted(cuts, pool.median_logfrac)
    ck = T.loc[clock].copy()
    ck["tercile"] = np.searchsorted(cuts, ck.median_logfrac)
    merged = [int(b) for b in range(6) if (pool[pool.bin == b].groupby("tercile").size().reindex(range(3), fill_value=0) < 10).any()]
    for b in merged:  # fewer than 10 pool genes in a (bin, tercile) cell: the terciles of that bin are merged
        pool.loc[pool.bin == b, "tercile"] = 0
        ck.loc[ck.bin == b, "tercile"] = 0
    return pool, ck, dict(n_combined_q05=int((T.q_combined < 0.05).sum()), n_expressed_ge_half=int((n_exp >= n_exp_min).sum()),
                          merged_bins=merged)


def draw_sets(pool, ck, n=N_SETS):
    """n sets, each clock gene replaced by a pool gene of its (bin, tercile) cell, without replacement within a set."""
    rng = np.random.default_rng([SEED, 1])
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
            rows.append(dict(set=i, position=pos, clock_gene=g, gene=c, core=pos in CORE_POS))
    return pd.DataFrame(rows)


# ---- per set ---------------------------------------------------------------------------------------------------

def set_waveform(pb, per, zt, genes):
    """Waveform of a gene set fitted on the baboon peripheral rows (hours = ZT, group = tissue)."""
    y = np.log(np.maximum(pb[genes].to_numpy(np.float64), FLOOR) / pb.depth.to_numpy(float)[:, None])
    return fit_waveform(y, zt.reindex(pb.donor).to_numpy(float), genes, group=pb.compartment.to_numpy())


def perms(zt):
    rng = np.random.default_rng([SEED, 2])
    return [pd.Series(rng.permutation(zt.to_numpy()), index=zt.index) for _ in range(N_PERM)]


def reads(pb, per, zt, genes, wf, label):
    """Passing organs and the three A395 slopes of one gene set under waveform wf (A395 pipeline)."""
    ds = cs.Dataset(pb[KEYS + list(genes)], anchors=per, waveform=wf, name=label, min_donors=5)
    obs = B.timed(ds, zt).abs_amplitude.reindex(per)
    null = np.array([B.timed(ds, p).abs_amplitude.reindex(per).to_numpy(float) for p in S["perms"]])
    passing = pd.Series(obs.to_numpy(float) > np.nanquantile(null, 0.95, axis=0), index=per)
    organs = list(passing.index[passing])
    out = dict(passing=len(organs))
    if len(organs) >= 3:
        T = B.organ_table(ds, zt, organs)
        r = B.readouts(T)
        out.update(used=len(r[0][1]), b_t=B.slope(r[0][1], r[0][2]), b_t_uncorrected=B.slope(r[1][1], r[1][2]),
                   scatter_slope=B.slope(r[2][1], r[2][2]))
    return out


def task(i):
    t0 = time.time()
    s = S["sets"][S["sets"].set == i].sort_values("position")
    genes = list(s.gene)
    pb = S["pb"][KEYS + genes]
    try:
        wf = set_waveform(pb, S["per"], S["zt"], genes)
        row = dict(set=i, **reads(pb, S["per"], S["zt"], genes, wf, f"a407 set {i}"))
        core = [genes[p] for p in CORE_POS]
        c = reads(pb[KEYS + core], S["per"], S["zt"], core, wf.subset(core), f"a407 set {i} core")
        row.update({f"core_{k}": v for k, v in c.items()})
    except Exception as e:
        row = dict(set=i, error=f"{type(e).__name__}: {e}"[:200])
    row["seconds"] = time.time() - t0
    return row


def power(jobs, n_sets):
    t0 = time.time()
    pb, per, zt, genes = load_genome()
    pool, ck, info = pool_table(pb, per, zt, genes)
    pool.reset_index(names="gene").to_csv(OUT / "pool.csv", index=False)
    sets = draw_sets(pool, ck)
    sets.to_csv(OUT / "sets.csv", index=False)
    print(f"pool {len(pool)} genes from {len(genes)} ({info}); clock-gene bins {dict(ck.bin)}; "
          f"peripheral tissues {len(per)}, rows {len(pb)} ({time.time() - t0:.0f}s)", flush=True)
    print("clock genes combined p:", {g: float(f"{ck.p_combined[g]:.3g}") for g in ck.index},
          "q:", {g: float(f"{ck.q_combined[g]:.3g}") for g in ck.index}, flush=True)
    print("per-bin pool counts", pool.bin.value_counts().sort_index().to_dict(),
          "per bin x tercile", pool.groupby(["bin", "tercile"]).size().to_dict(), flush=True)
    S.update(pb=pb, per=per, zt=zt, sets=sets[sets.set < n_sets], perms=perms(zt))
    rows = []
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, range(n_sets))):
            rows.append(r)
            if (k + 1) % 5 == 0 or k + 1 == n_sets:
                print(f"{k + 1}/{n_sets} sets {time.time() - t0:.0f}s (last {r['seconds']:.0f}s)", flush=True)
    R = pd.DataFrame(rows).sort_values("set")
    sfx = "_check" if n_sets < N_SETS else ""
    R.to_csv(OUT / f"companions_unread{sfx}.csv", index=False)  # b_t and slopes: not read before the score stage
    n_ok = int((R.passing >= MIN_ORGANS).sum())
    P = pd.DataFrame([dict(sets=len(R), sets_ge15_organs=n_ok, needed=len(R) // 2, errors=int(R.get("error", pd.Series(dtype=object)).notna().sum()),
                           pool_size=len(pool), **{f"pool_bin{b}": int(n) for b, n in pool.bin.value_counts().sort_index().items()},
                           runtime_s=round(time.time() - t0))])
    P.to_csv(OUT / f"power{sfx}.csv", index=False)
    print(P.T.to_string(), flush=True)
    if "error" in R:
        print("errors:", R.error.dropna().unique()[:5].tolist(), flush=True)
    print(f"passing organs per set: min {R.passing.min()}, median {R.passing.median()}, max {R.passing.max()}")


# ---- score -----------------------------------------------------------------------------------------------------

def score():
    pb, per, zt, genes = load_genome()
    S.update(perms=perms(zt))
    R = pd.read_csv(OUT / "companions_unread.csv")
    P = pd.read_csv(OUT / "power.csv").iloc[0]
    ok = R[R.passing >= MIN_ORGANS].dropna(subset=["b_t"])
    rows = []
    clock_tab = {}
    for lab, g in (("14 genes", list(WF.genes)), ("core", [x for x in WF.genes if x not in TRIO])):
        d = pb[KEYS + g]
        wf = set_waveform(d, per, zt, g)
        ds = cs.Dataset(d, anchors=per, waveform=wf, name=f"a407 clock {lab}", min_donors=5)
        obs = B.timed(ds, zt).abs_amplitude.reindex(per)
        null = np.array([B.timed(ds, p).abs_amplitude.reindex(per).to_numpy(float) for p in S["perms"]])
        organs = list(per[i] for i in np.flatnonzero(obs.to_numpy(float) > np.nanquantile(null, 0.95, axis=0)))
        T = B.organ_table(ds, zt, organs)
        clock_tab[lab] = T
        for j, (name, x, y) in enumerate(B.readouts(T)):
            bs = B.boot(x, y, np.random.default_rng([SEED, 3, j]))
            rows.append(dict(readout=lab, waveform="baboon-fitted", slope=name, value=B.slope(x, y),
                             lo90=np.quantile(bs, 0.05), hi90=np.quantile(bs, 0.95), organs=len(x)))
    gt = pd.read_csv(HERE.parents[1] / "output" / "organs_14_genes.csv", index_col="organ")
    for name, x, y in B.readouts(gt):
        rows.append(dict(readout="14 genes", waveform="GTEx (A395)", slope=name, value=B.slope(x, y), organs=len(x)))
    T = pd.DataFrame(rows)
    comp = {}
    for col, name in (("b_t", "b_t corrected"), ("b_t_uncorrected", "b_t uncorrected"), ("scatter_slope", "scatter slope corrected")):
        v = ok[col].dropna()
        c = T[(T.readout == "14 genes") & (T.waveform == "baboon-fitted") & (T.slope == name)].value.iloc[0]
        comp[name] = dict(companions=len(v), q05=v.quantile(0.05), median=v.median(), q95=v.quantile(0.95),
                          clock=c, pct_companions_below_clock=100 * np.mean(v <= c))
    C = pd.DataFrame(comp).T
    h = comp["b_t corrected"]
    C["H407.1_clock_below_q05"] = [bool(h["clock"] < h["q05"]) if n == "b_t corrected" else np.nan for n in C.index]
    C["median_companion_b_t_below_0.5"] = bool(h["median"] < 0.5)
    T.to_csv(OUT / "tests_clock.csv", index=False)
    C.to_csv(OUT / "tests.csv")
    print(T.round(3).to_string(index=False))
    print(C.round(3).to_string())


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true", help="3 sets, shapes and timings only")
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        power(a.jobs, 3 if a.check else N_SETS)
    else:
        score()


if __name__ == "__main__":
    main()
