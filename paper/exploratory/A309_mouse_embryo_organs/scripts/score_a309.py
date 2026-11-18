"""A309 scoring: between-embryo clock covariation of five cell types in the late mouse embryo (H309.1), the A307
glucocorticoid / waveform statistics where it passes (H309.2), and all-age descriptives.

Registered in paper/docs/registered_experiments.md (A309). Only readouts with an adult pass rate >= 0.8 at the full clock
(`power_passrates.csv`) are scored.

    python score_a309.py [--n-perm 10000] [--n-pseudo 200] [--jobs 4]
Outputs (../output/): l1.csv, verdicts.csv, gc_pattern.csv, gc_genes.csv, all_ages.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

import lib_a309 as L
from clocksync.config import limit_threads
from clocksync.paired import pseudo_waveforms
from common import data
from common import development as dev

sys.path.insert(0, str(L.HERE.parents[1] / "A307_fetal_gc_pattern" / "scripts"))
S = {}
_log_mean = dev.log_mean


def _log_mean_keep(counts, depth, min_count: float = 1.0):
    """`development.log_mean` that keeps the waveform genes whatever their count (the >= min_count rule then applies
    to the replacement pool only); plumbing fix: CIART is below 1 count in a late-window readout."""
    full = _log_mean(counts, depth, 0.0)
    keep = (counts.mean(0) >= min_count) | counts.columns.isin(L.GENES)
    return full[keep[full.index].to_numpy()]


dev.log_mean = _log_mean_keep


def task(r):
    j = list(L.READOUTS).index(r)
    o = S["obs"][S["obs"].readout == r]
    return r, dev.l1_test(S["counts"], o, L.cov, L.WF, S["n_perm"], S["n_pseudo"], [L.SEED, 1, j], "block")


def gc_pattern(r, n_perm, n_pseudo):
    """A307 L and W on readout r (within-batch nulls, matched companions)."""
    import score_a307 as G
    j = list(L.READOUTS).index(r)
    o = S["obs"][S["obs"].readout == r]
    counts = S["counts"].loc[o.index]
    y, groups = G.adjusted(counts, o, L.cov, L.GENES + L.GC)
    Y, g = y[:, :len(L.GENES)], G.gc_score(y[:, len(L.GENES):])
    si = [L.GENES.index(x) for x in L.SI]
    co = [L.GENES.index(x) for x in L.CORE]
    phi = dict(zip(L.GENES, L.WF.psi))
    cosm = np.cos(np.subtract.outer([phi[x] for x in L.SI], [phi[x] for x in L.CORE]))
    rng = np.random.default_rng([L.SEED, 3, j, 0])
    nL = np.array([G.L_stat(Y[:, si], G.perm_vec(g, groups, rng)) for _ in range(n_perm)])
    rng = np.random.default_rng([L.SEED, 3, j, 2])
    nW = np.array([G.W_stat(dev._permuted(Y, groups, rng), si, co, cosm) for _ in range(n_perm)])
    lm = dev.log_mean(counts, o.depth)
    ws = pseudo_waveforms(lm, n_pseudo + 100, np.random.default_rng([L.SEED, 3, j, 1]), L.WF.drop(L.CORE))
    cL = [G.L_stat(G.adjusted(counts, o, L.cov, list(w.genes))[0], g) for w in ws if not set(w.genes) & set(L.GC)][:n_pseudo]
    cW = [G.W_stat(G.adjusted(counts, o, L.cov, list(w.genes))[0], si, co, cosm)
          for w in pseudo_waveforms(lm, n_pseudo, np.random.default_rng([L.SEED, 3, j, 3]), L.WF)]
    Lo, Wo = G.L_stat(Y[:, si], g), G.W_stat(Y, si, co, cosm)
    rows = [dict(readout=r, statistic=s, value=v, p=float((1 + np.sum(n >= v)) / (1 + len(n))),
                 matched_median=float(np.median(c)), matched_q95=float(np.quantile(c, .95)),
                 matched_pct=float(100 * np.mean(np.array(c) < v)), n_sets=len(c))
            for s, v, n, c in (("L", Lo, nL, cL), ("W", Wo, nW, cW))]
    genes = [dict(readout=r, gene=x, r_gc=float(v)) for x, v in zip(L.GENES, G.corr_with(Y, g))]
    return rows, genes


def all_ages(n_perm):
    pb, don = data.load("mouse_qiu")
    pb = pb.join(don[["day"]], on="donor")
    out = []
    for k, (r, ct) in enumerate(L.READOUTS.items()):
        x = pb[(pb.compartment == ct) & (pb.n_cells >= L.MIN_CELLS) & (pb.day < 19.5)].reset_index(drop=True)
        tab = dev.adjust(x, L.GENES, lambda g: dev.design(g, "day", 2))[0]
        c, p, q = dev.coherence(tab, L.WF, n_perm, np.random.default_rng([L.SEED, 4, k]))
        out.append(dict(readout=r, embryos=len(x), day_min=float(x.day.min()), day_max=float(x.day.max()), C=c, p=p,
                        null_q95=q))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    rows = L.late_rows()
    counts = L.genome(rows)
    limit_threads(1)
    rates = pd.read_csv(L.OUT / "power_passrates.csv").pivot(index="readout", columns="setting", values="pass_rate")
    run = [r for r in L.READOUTS if rates.loc[r, "a1"] >= 0.8]
    S.update(obs=rows, counts=counts, n_perm=a.n_perm, n_pseudo=a.n_pseudo)
    res = {}
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for r, out in pool.imap_unordered(task, run):
            res[r] = out
            print(r, "done", f"{time.time() - t0:.0f}s", flush=True)
    l1 = pd.DataFrame([dict(readout=r, cell_type=L.READOUTS[r], adult_a1=rates.loc[r, "a1"],
                            adult_a075=rates.loc[r, "a0.75"], adult_a05=rates.loc[r, "a0.5"],
                            adult_shuffled=rates.loc[r, "shuffled"], **res.get(r, {})) for r in L.READOUTS])
    sc = l1.readout.isin(run)
    l1.loc[sc, "p_holm"] = dev.holm(l1.p[sc].to_numpy())
    l1["verdict"] = np.where(~sc, "not estimable (adult power < 0.8)",
                             np.where(l1.p_holm >= 0.05, "fail",
                                      np.where(l1.clock_specific.fillna(False).astype(bool), "PASS", "not clock-specific")))
    l1.to_csv(L.OUT / "l1.csv", index=False)
    pd.set_option("display.width", 250)
    print(l1.drop(columns=["cell_type"]).round(4).to_string(index=False), flush=True)
    gp, gg = [], []
    for r in l1.readout[l1.verdict == "PASS"]:
        t, g = gc_pattern(r, a.n_perm, a.n_pseudo)
        gp += t
        gg += g
    if gp:
        gp = pd.DataFrame(gp)
        for s in ("L", "W"):
            m = gp.statistic == s
            gp.loc[m, "p_holm"] = dev.holm(gp.p[m].to_numpy())
        gp["pass"] = (gp.p_holm < 0.05) & (gp.value > gp.matched_q95)
        gp.to_csv(L.OUT / "gc_pattern.csv", index=False)
        pd.DataFrame(gg).to_csv(L.OUT / "gc_genes.csv", index=False)
        print(gp.round(4).to_string(index=False), pd.DataFrame(gg).pivot(index="gene", columns="readout", values="r_gc")
              .round(3).to_string(), sep="\n", flush=True)
    aa = all_ages(a.n_perm)
    aa.to_csv(L.OUT / "all_ages.csv", index=False)
    print(aa.round(4).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
