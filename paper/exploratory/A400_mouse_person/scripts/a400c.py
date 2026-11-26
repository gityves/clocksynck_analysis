"""A400c (descriptive check, registered 2026-11-26): is the mouse time-of-day radius a projection artefact of reading
mouse clock genes with the human GTEx waveform, and is it beyond noise?

Waveforms: H (human GTEx, as A400), Mx (mouse, cross-fitted: `fit_waveform` on one mouse of each ZT x arm cell reads
the other), Mf (mouse, fitted on all mice; sensitivity). Per waveform and readout: R^2 of the half-host radius on the
ZT x arm cells and incremental R^2 of ZT harmonics 1-4; ellipticity (second-harmonic amplitude / mean radius) with a
mouse bootstrap; F-test of the cells; A400 leverage under V0 / V1 / V2 with the mice-permuted null (H companions: the
A400b values; Mx companions: 50 + 50 sets read with the fold waveforms). Mouse vs human gene amplitudes and phases (Mf,
200 mouse bootstraps).

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A400_mouse_person/scripts/a400c.py --jobs 4
    python a400c.py --check      # 2 sets, 20 permutations / bootstraps; nothing written
Outputs (`../output/`): c_waveforms.csv, c_radius.csv, c_leverages.csv, c_summary.csv.
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
from scipy.stats import f as fdist  # noqa: E402

import a385b as R  # noqa: E402
import a394 as M  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform, fit_waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
VARIANTS = ("V0 A400", "V1 ZT harmonics 1-4", "V2 ZT x arm cells")
TRIO = ["PER1", "PER2", "NR1D1"]
N_SETS, N_PERM, N_BOOT, N_GENE_BOOT = 50, 1000, 1000, 200
S = {}


# ---- waveforms -------------------------------------------------------------------------------------------------------
def fit_on(mice):
    """Mouse waveform fitted on the samples of `mice` (14 genes, peripheral tissues, group tissue x arm)."""
    pb = S["mouse"]
    d = pb[pb.donor.isin(mice)]
    y = np.log((d[list(R.WF.genes)].to_numpy(float) + 0.5) / d.depth.to_numpy(float)[:, None])
    return fit_waveform(y, S["zt"].reindex(d.donor).to_numpy(float), R.WF.genes,
                        group=d.compartment.astype(str) + "|" + d.stratum.astype(str), name="mouse")


def as_set(wf, genes, rd, name):
    """Waveform with `genes` at the positions of `wf` (companions: as `pseudo_waveforms`); core = core positions."""
    w = Waveform(tuple(genes), np.asarray(wf.amp), np.asarray(wf.psi), None, None, name)
    return w if rd == "14 genes" else Waveform(tuple(w.genes[k] for k in R.CORE_IDX), w.amp[R.CORE_IDX],
                                                w.psi[R.CORE_IDX], None, None, name)


# ---- states ----------------------------------------------------------------------------------------------------------
def states(w, mice):
    pb = S["mouse"]
    d = pb[pb.donor.isin(mice)]
    ds = cs.Dataset(d[M.KEYS + list(w.genes)], anchors=S["per"], waveform=w, name=w.name, min_donors=5)
    return R.states(ds, S["halvings"])


def cross_states(genes, rd, name):
    """Mx: each fold read with the waveform fitted on the other fold; fold indicator appended to the covariates."""
    parts = []
    for f in (0, 1):
        w = as_set(S["wf_fold"][1 - f], genes, rd, name)
        parts.append((f, states(w, S["folds"][f])))
    out = []
    for h in range(len(S["halvings"])):
        rows = [(p[h], f) for f, p in parts]
        j = np.concatenate([r[0] for r, _ in rows])
        a = np.concatenate([r[1] for r, _ in rows])
        b = np.concatenate([r[2] for r, _ in rows])
        cov = np.vstack([np.column_stack([r[3], np.full(len(r[0]), f)]) for r, f in rows])
        out.append((j, a, b, cov))
    return out


# ---- statistics ------------------------------------------------------------------------------------------------------
def harm(j, ks):
    w = S["zt"].reindex(j).to_numpy(float) * 2 * np.pi / 24
    return np.column_stack([g(k * w) for k in ks for g in (np.cos, np.sin)]) if ks else np.empty((len(j), 0))


def dummies(v):
    return pd.get_dummies(pd.Series(v), drop_first=True).to_numpy(float)


def leverage(st, variant, rng=None):
    L = []
    for j, a, b, cov in st:
        arm = S["arm"].reindex(j).to_numpy()
        if rng is not None:
            b = b.copy()
            for p in np.unique(arm):
                m = np.flatnonzero(arm == p)
                b[m] = b[rng.permutation(m)]
        X = np.column_stack([cov, harm(j, range(1, 5))]) if variant == VARIANTS[1] else cov
        grp = S["cell"].reindex(j).to_numpy() if variant == VARIANTS[2] else arm
        L.append(R.pooled_leverage(a, b, X, grp)[0])
    return float(np.nanmedian(L))


def ols_rss(y, X):
    return float(np.sum((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2))


def radius_stats(st):
    """Cell R^2, incremental harmonic R^2, cell F-test p (medians over halvings and halves)."""
    r2c, inc, pf = [], [], []
    for j, a, b, cov in st:
        base = np.column_stack([np.ones(len(j)), dummies(S["arm"].reindex(j).to_numpy()), cov[:, 2:]])
        cells = np.column_stack([np.ones(len(j)), dummies(S["cell"].reindex(j).to_numpy()), cov[:, 2:]])
        for z in (a, b):
            r = np.abs(z)
            tss = float(np.sum((r - r.mean()) ** 2))
            rss0, rss1 = ols_rss(r, base), ols_rss(r, cells)
            r2c.append(1 - rss1 / tss)
            rk1 = np.linalg.matrix_rank(cells)
            df1, dfd = rk1 - np.linalg.matrix_rank(base), len(r) - rk1
            pf.append(float(fdist.sf(((rss0 - rss1) / df1) / (rss1 / dfd), df1, dfd)))
            prev, row = rss0, []
            for k in range(1, 5):
                rk = ols_rss(r, np.column_stack([base, harm(j, range(1, k + 1))]))
                row.append((prev - rk) / tss)
                prev = rk
            inc.append(row)
    inc = np.median(np.array(inc), axis=0)
    return dict(r2_cells=float(np.median(r2c)), f_p_cells=float(np.median(pf)),
                **{f"inc_r2_h{k + 1}": float(inc[k]) for k in range(4)})


def ellipticity(st, idx=None):
    """Second ZT harmonic amplitude of the half-A radius / mean radius, median over halvings; idx: bootstrap mice."""
    e = []
    for j, a, _, cov in st:
        if idx is not None:
            k = np.concatenate([np.flatnonzero(j == m) for m in idx])
            j, a, cov = j[k], a[k], cov[k]
        r = np.abs(a)
        X = np.column_stack([np.ones(len(j)), dummies(S["arm"].reindex(j).to_numpy()), cov[:, 2:], harm(j, (1, 2))])
        c = np.linalg.lstsq(X, r, rcond=None)[0][-2:]
        e.append(np.hypot(*c) / r.mean())
    return float(np.median(e))


def task(job):
    warnings.filterwarnings("ignore")
    kind, i, rd, wname, genes = job
    try:
        if wname == "Mx":
            st = cross_states(genes, rd, f"{kind}{i}")
        else:
            st = states(as_set(S["wf"][wname], genes, rd, f"{kind}{i}"), S["mice"])
            st = [(j, a, b, np.column_stack([cov, np.zeros(len(j))])) for j, a, b, cov in st]  # no folds
        out = [dict(kind=kind, set=i, readout=rd, waveform=wname, variant=v, leverage=leverage(st, v)) for v in VARIANTS]
        if kind == "clock":
            rng = np.random.default_rng([400, 5])
            for o in out:
                null = np.array([leverage(st, o["variant"], rng) for _ in range(S["n_perm"])])
                o.update(p_perm=float((1 + np.sum(null >= o["leverage"])) / (1 + len(null))),
                         null_q95=float(np.quantile(null, 0.95)))
            rad = radius_stats(st)
            rb = np.random.default_rng([400, 6])
            mice = np.unique(np.concatenate([j for j, *_ in st]))
            boots = [ellipticity(st, rb.choice(mice, len(mice))) for _ in range(S["n_boot"])]
            rad.update(kind=kind, readout=rd, waveform=wname, ellipticity=ellipticity(st),
                       ell_lo90=float(np.quantile(boots, 0.05)), ell_hi90=float(np.quantile(boots, 0.95)),
                       chance_r2=(S["cell"].nunique() - 1) / (len(mice) - 1))
            return out, [rad]
        return out, []
    except Exception as e:  # counted, not dropped silently
        return [dict(kind=kind, set=i, readout=rd, waveform=wname, variant=None, error=f"{type(e).__name__}: {e}"[:200])], []


def gene_table(n_boot):
    """Mouse (Mf) vs GTEx per gene: relative amplitude difference and centred phase difference, mouse bootstrap."""
    def compare(w):
        ra = np.log(np.asarray(w.amp) / np.mean(w.amp)) - np.log(R.WF.amp / np.mean(R.WF.amp))
        dp = np.angle(np.exp(1j * (np.asarray(w.psi) - R.WF.psi)))
        dp = np.angle(np.exp(1j * (dp - np.angle(np.mean(np.exp(1j * dp))))))
        return ra, dp * 24 / (2 * np.pi)
    ra, dp = compare(S["wf"]["Mf"])
    rng = np.random.default_rng([400, 7])
    bb = [compare(fit_on(rng.choice(S["mice"], len(S["mice"])))) for _ in range(n_boot)]
    A, P = np.array([b[0] for b in bb]), np.array([b[1] for b in bb])
    return pd.DataFrame(dict(gene=R.WF.genes, log_rel_amp_mouse_minus_gtex=ra, lo90=np.quantile(A, 0.05, 0),
                             hi90=np.quantile(A, 0.95, 0), phase_diff_h=dp, phase_lo90=np.quantile(P, 0.05, 0),
                             phase_hi90=np.quantile(P, 0.95, 0), core=[g not in TRIO for g in R.WF.genes]))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    pb, per, _ = M.mouse_table()
    pb = pb[pb.compartment.isin(per)].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    lm = M.log_mean(pb)
    n = 2 if a.check else N_SETS
    sets = {k: pseudo_waveforms(lm, N_SETS, np.random.default_rng([394, t]), R.WF, kind=k)[:n]
            for k, t in (("matched", 2), ("responsive", 3))}
    genes = sorted({g for ws in sets.values() for w in ws for g in w.genes} | set(R.WF.genes))
    arm = pb.groupby("donor").stratum.first()
    zt = pd.Series(arm.index.str.extract(r"_ZT(\d+)_", expand=False).astype(float), index=arm.index)
    cell = pd.Series(arm.to_numpy() + "_" + zt.astype(int).astype(str).to_numpy(), index=arm.index)
    fold = pd.Series(0, index=arm.index)
    for ms in cell.groupby(cell).groups.values():
        fold[sorted(ms)[1::2]] = 1
    S.update(mouse=pb[M.KEYS + genes].copy(), per=per, arm=arm, zt=zt, cell=cell, mice=np.asarray(arm.index),
             folds=[np.asarray(fold.index[fold == f]) for f in (0, 1)], halvings=R.halves(per, [400, 1]),
             n_perm=20 if a.check else N_PERM, n_boot=20 if a.check else N_BOOT)
    S["wf"] = {"H": R.WF, "Mf": fit_on(S["mice"])}
    S["wf_fold"] = [fit_on(S["folds"][f]) for f in (0, 1)]
    print(f"{arm.size} mice, {cell.nunique()} cells, folds {[len(f) for f in S['folds']]} ({time.time() - t0:.0f}s)", flush=True)

    jobs = [("clock", 0, rd, wn, R.WF.genes) for wn in ("H", "Mx", "Mf") for rd in ("core", "14 genes")]
    jobs += [(k, j, rd, "Mx", w.genes) for k, ws in sets.items() for j, w in enumerate(ws) for rd in ("14 genes", "core")]
    rows, rad = [], []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, (o, r) in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows += o
            rad += r
            if (k + 1) % 10 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    C, Rd = pd.DataFrame(rows), pd.DataFrame(rad)
    h0 = C[(C.kind == "clock") & (C.waveform == "H") & (C.variant == VARIANTS[0])].set_index("readout").leverage
    print(f"A400 reproduction (H, V0): core {h0['core']:.4f} vs 0.4671; 14 genes {h0['14 genes']:.4f} vs 0.7857", flush=True)
    G = gene_table(5 if a.check else N_GENE_BOOT)
    if a.check:
        print(f"check ok; failed {int(C.variant.isna().sum())}", flush=True)
        return
    prev = pd.read_csv(OUT / "b_leverages.csv")
    prev = prev[prev.kind != "clock"].assign(waveform="H").replace(
        {"variant": {"V0 unadjusted": VARIANTS[0], "V2 ZT x arm cells": VARIANTS[2], "V1 ZT harmonics 1-4": VARIANTS[1]}})
    comp = pd.concat([C[C.kind != "clock"], prev[prev.variant.isin(VARIANTS)]])
    summ = []
    for wn in ("H", "Mx", "Mf"):
        for rd in ("core", "14 genes"):
            for v in VARIANTS:
                clk = C[(C.kind == "clock") & (C.waveform == wn) & (C.readout == rd) & (C.variant == v)].iloc[0]
                cw = comp[(comp.waveform == ("H" if wn == "H" else "Mx")) & (comp.readout == rd) & (comp.variant == v)]
                q = {k: float(np.nanquantile(cw[cw.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}
                summ.append(dict(waveform=wn, readout=rd, variant=v, clock=clk.leverage, p_perm=clk.p_perm,
                                 null_q95=clk.null_q95, matched_q95=q["matched"], responsive_q95=q["responsive"],
                                 above_both=bool(clk.leverage > max(q.values()))))
    OUT.mkdir(parents=True, exist_ok=True)
    C.to_csv(OUT / "c_leverages.csv", index=False)
    Rd.to_csv(OUT / "c_radius.csv", index=False)
    G.to_csv(OUT / "c_waveforms.csv", index=False)
    Sm = pd.DataFrame(summ)
    Sm.to_csv(OUT / "c_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(Rd.round(4).to_string(index=False))
    print(Sm.round(4).to_string(index=False))
    print(G.round(3).to_string(index=False), f"\nfailed: {int(C.variant.isna().sum())} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
