"""A316 (registered): within-tissue cell-to-cell clock dispersion of the aorta during re-entrainment after a 6 h advance.

Auerbach 2026 aorta (A201 loader). Per library x cell type (>= 100 cells): D from `clocksync.cells.dispersion`
(cross-gene moments of split gene halves, log scale, 30 splits); the library's pseudobulk clock state b_l (same
projector on the log mean expression). Per condition x cell type: D_bar, coherent power A^2 (spread of b_l around the
sex means, minus Poisson variance), F = A^2 / (A^2 + D_bar), P = A^2 + D_bar. Library jackknife within condition.

    python score_a316.py --stage feasibility   # aligned arm only: UMI, D vs null, placebo splits
    python score_a316.py --stage score         # misaligned vs aligned, Bmal1 KO descriptive
Outputs (`../output/`): feasibility.csv, placebo.csv | libraries.csv, tests.csv, celltypes.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A201_aorta_jetlag" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a201 as A201  # noqa: E402
from clocksync.cells import dispersion as DSP  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
MAIN, EXTRA = ["SMC", "Fibroblast"], ["EC", "Macrophage"]
MIN_CELLS, N_SPLITS, SEED = 100, 30, 316
BOUND, SE_MAX = 0.26, 0.105


def design(genes):
    wf = cs.gtex()
    idx = {g.capitalize(): i for i, g in enumerate(wf.genes)}
    X = np.asarray(wf.X, float)
    return np.array([X[idx[g]] for g in genes])


def library_rows(o, Xc, X, conds, rng, null=False):
    """One row per library x cell type: D (and its permuted-acrophase null mean if `null`), pseudobulk state, UMI."""
    L = DSP._projector(X, np.ones(len(X)))
    rows = []
    for (lib, ct), idx in o[o.cond.isin(conds)].groupby(["sample_id", "cell_type"], observed=True).groups.items():
        if ct not in MAIN + EXTRA or len(idx) < MIN_CELLS:
            continue
        i = np.asarray(idx)
        C, N = Xc[i], o.total_counts.to_numpy(float)[i]
        d = DSP.dispersion(C, N, X, n_splits=N_SPLITS, rng=rng)
        tot = C.sum(0)
        lp = np.log((tot + 0.5) / N.sum())
        r = o.iloc[i[0]]
        row = dict(library=lib, cell_type=ct, cond=r.cond, sex=r.sex, zt=r.zt, cells=len(i), D=d["D"],
                   umi_median=float(np.median(C.sum(1))), b_x=np.nan, b_y=np.nan, b_var=float(np.trace(L @ np.diag(1 / (tot + 0.5)) @ L.T)))
        row["_lp"] = lp
        if null:
            row["D_null"] = float(np.mean([DSP.dispersion(C, N, X[rng.permutation(len(X))], n_splits=N_SPLITS, rng=rng)["D"]
                                           for _ in range(20)]))
        rows.append(row)
    t = pd.DataFrame(rows)
    # pseudobulk states: library log means centred over the libraries of the same condition x cell type
    for (cond, ct), g in t.groupby(["cond", "cell_type"]):
        M = np.vstack(g._lp.to_numpy())
        b = (M - M.mean(0)) @ L.T
        t.loc[g.index, "b_x"], t.loc[g.index, "b_y"] = b[:, 0], b[:, 1]
    return t.drop(columns="_lp")


def stats(t):
    """F and P of one condition x cell type table (libraries)."""
    b = t[["b_x", "b_y"]].to_numpy()
    b = b - t.groupby("sex")[["b_x", "b_y"]].transform("mean").to_numpy()
    A2 = float(np.mean(np.sum(b ** 2, 1)) - t.b_var.mean())
    D = float(t.D.mean())
    return dict(A2=A2, D=D, F=A2 / (A2 + D) if A2 + D > 0 else np.nan, P=A2 + D)


def log_ratio(ta, tb, key):
    sa, sb = stats(ta), stats(tb)
    return np.log(sa[key] / sb[key]) if sa[key] > 0 and sb[key] > 0 else np.nan


def jack(ta, tb, key):
    """Pooled (mean over MAIN) log ratio a / b and its library jackknife SE (libraries of both arms left out in turn)."""
    def pooled(xa, xb):
        return float(np.nanmean([log_ratio(xa[xa.cell_type == c], xb[xb.cell_type == c], key) for c in MAIN]))
    est = pooled(ta, tb)
    reps = []
    for side, tt in (("a", ta), ("b", tb)):
        libs = tt.library.unique()
        for l in libs:
            reps.append(pooled(ta[ta.library != l] if side == "a" else ta, tb[tb.library != l] if side == "b" else tb))
    reps = np.array(reps)
    n = len(reps)
    se = float(np.sqrt((n - 1) / n * np.sum((reps - reps.mean()) ** 2)))
    return est, se


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["feasibility", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    o, Xc, genes = A201.load()
    X = design(genes)
    rng = np.random.default_rng(SEED)
    print(f"loaded {len(o)} cells, {len(genes)} genes ({time.time() - t0:.0f}s)", flush=True)
    if a.stage == "feasibility":
        t = library_rows(o, Xc, X, ["ali"], rng, null=True)
        feas = []
        for ct in MAIN:
            s = t[t.cell_type == ct]
            feas.append(dict(cell_type=ct, libraries=len(s), umi_median=float(s.umi_median.median()), D=float(s.D.mean()),
                             D_null=float(s.D_null.mean()), D_over_null=bool(s.D.mean() > s.D_null.mean())))
        libs = t[["library", "sex", "zt"]].drop_duplicates()
        pl = []
        tries = 0
        while len(pl) < 50 and tries < 5000:
            tries += 1
            h = rng.permutation(len(libs)) % 2 == 0
            A, B = libs[h], libs[~h]
            if set(A.zt) != set(libs.zt) or set(B.zt) != set(libs.zt) or A.sex.nunique() < 2 or B.sex.nunique() < 2:
                continue
            ta, tb = t[t.library.isin(A.library)], t[t.library.isin(B.library)]
            est, se = jack(ta, tb, "F")
            pl.append(dict(split=len(pl), logF=est, se=se))
        pl = pd.DataFrame(pl)
        pl.to_csv(OUT / "placebo.csv", index=False)
        f = pd.DataFrame(feas)
        res = dict(umi_ok=bool((f.umi_median >= 2).all()), D_over_null=bool(f.D_over_null.all()),
                   placebo_splits=len(pl), placebo_mean=float(pl.logF.mean()), placebo_sd=float(pl.logF.std(ddof=1)),
                   se_median=float(pl.se.median()))
        res["gate"] = bool(res["umi_ok"] and res["D_over_null"] and res["se_median"] <= SE_MAX and abs(res["placebo_mean"]) <= 0.1)
        f.to_csv(OUT / "feasibility.csv", index=False)
        pd.Series(res).to_csv(OUT / "feasibility_gate.csv", header=["value"])
        print(f.round(3).to_string(index=False))
        print(pd.Series(res).to_string(), flush=True)
        return
    t = library_rows(o, Xc, X, ["ali", "mis", "ko"], rng)
    t.to_csv(OUT / "libraries.csv", index=False)
    ali, mis, ko = (t[t.cond == c] for c in ("ali", "mis", "ko"))
    tests = []
    for key in ("F", "P"):
        est, se = jack(mis, ali, key)
        tests.append(dict(test=f"log {key} mis/ali", est=est, se=se, lo90=est - 1.645 * se, hi90=est + 1.645 * se))
    F = tests[0]
    tests.append(dict(test="H316.1 dephasing", passes=bool(F["hi90"] < 0)))
    tests.append(dict(test="H316.2 block", passes=bool(F["lo90"] > -BOUND and F["hi90"] < BOUND)))
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    rows = []
    for ct in MAIN + EXTRA:
        for cond, tt in (("ali", ali), ("mis", mis), ("ko", ko)):
            s = tt[tt.cell_type == ct]
            if len(s) >= 3:
                rows.append(dict(cell_type=ct, cond=cond, libraries=len(s), **stats(s), umi_median=float(s.umi_median.median())))
    ctab = pd.DataFrame(rows)
    ctab.to_csv(OUT / "celltypes.csv", index=False)
    print(ctab.round(4).to_string(index=False))
    print(pd.DataFrame(tests).round(4).to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
