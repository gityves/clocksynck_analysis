"""A267 (registered in `paper/docs/registered_experiments.md`; follow-up of A265): is the pairing of COFE's tumour and normal
labels a clock signal or a patient-shared non-clock effect?

Data and pairing statistic as A265 (`lib_a265`, `score_a265.pooled_tests`: R_c = |mean e^{i (phi_T - phi_N)}| per cancer,
pooled with weights n pairs, 10,000 within-cancer permutations of the normals for the observed sets).
- `clock`: COFE (`run_cofe`, as published: standardised features, 3 x 5-fold CV, 5 restarts) on the 14 GTEx clock genes
  of the tumours, log2(cpm + 1) of the protein-coding counts, no expression filter; normals projected (COFE's rule).
  Companions: COFE on 100 expression-matched and 100 responsive 14-gene sets per cancer (`pseudo_waveforms`, matching
  variable = mean log fraction over the cancer's tumours, genes with mean count >= 1; set j of every cancer pooled);
  each set's pooled R_c. H267.1: the clock-gene pooled R_c has p < 0.05 and lies above the q95 of both families.
- `axes`: generic-axis companion of the genome-wide COFE pairing of A265 (H267.2): on the A265 COFE matrices (standardised
  over the tumours), the top 10 principal components of the tumours, normals projected; for each of the 45 pairs of
  components the angle atan2(PC_b, PC_a) of tumours and normals and its pooled R_c. H267.2: the genome-wide COFE pooled
  R_c of A265 (0.302) lies above the q95 of the 45 pairs.
- descriptive: genome-wide COFE without the 365 genes of the A265 simulation program (clock and GTEx-rhythmic genes);
  circular correlation of clock-gene COFE tumour phases with the clocksync tumour state phase (`cohort_states`).

    python score_a267.py [--jobs 6] [--n-sets 100]
Outputs (`../output/`): a267_sets.csv (one row per gene set: kind, pooled R_c), a267_axes.csv, a267_tests.csv,
a267_per_cancer.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import lib_a265 as L  # noqa: E402
import score_a265 as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import cohort_states, pseudo_waveforms  # noqa: E402

OUT = A.OUT
SEED, W = 267, A.W
CLOCK = A.CLOCK
D = {}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def cofe_pairs(c, genes, seed):
    """Tumour and normal COFE phases (radians) of the pairs of cancer c on `genes`."""
    d = D[c]
    lt, ln = D["log_t"][c][list(genes)], D["log_n"][c][list(genes)]
    lt = lt.loc[:, lt.std() > 0]
    r = A.run_cofe(lt, seed=seed, ncores=1, test=ln[lt.columns].loc[d["pairs"]])
    return r["phase"].loc[d["pairs"]].to_numpy() * W, r["test_phase"].loc[d["pairs"]].to_numpy() * W


def pooled_R(per):
    w = np.array([len(v[0]) for v in per.values()], float)
    return float(np.average([abs(np.exp(1j * (pt - pn)).mean()) for pt, pn in per.values()], weights=w))


def one_set(args):
    kind, j, sets = args
    per = {c: cofe_pairs(c, g, SEED + 1000 * j + i) for i, (c, g) in enumerate(sets.items())}
    return dict(kind=kind, set=j, R=pooled_R(per))


def perm_p(per, n_perm, rng):
    obs = pooled_R(per)
    w = np.array([len(v[0]) for v in per.values()], float)
    null = np.zeros(n_perm)
    for (pt, pn), wi in zip(per.values(), w):
        perms = np.argsort(rng.random((n_perm, len(pt))), axis=1)
        null += wi * np.abs(np.exp(1j * (pt[None, :] - pn[perms])).mean(1))
    null /= w.sum()
    return obs, float(null.mean()), float((1 + np.sum(null >= obs)) / (1 + n_perm))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--n-sets", type=int, default=100)
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng([SEED, 0])
    D["log_t"], D["log_n"] = {}, {}
    for c in L.CANCERS:
        D[c] = L.cancer(c)
        D["log_t"][c], D["log_n"][c] = L.log2cpm(D[c]["counts_t"]), L.log2cpm(D[c]["counts_n"])
        log("loaded", c)
    clock_cols = {c: L.B.T.clock_columns(D[c]["counts_t"].columns) for c in L.CANCERS}
    # gene sets: clock and companions (set j: one draw per cancer)
    sets = [("clock", 0, {c: clock_cols[c] for c in L.CANCERS})]
    for k, kind in enumerate(("matched", "responsive"), 1):
        per_c = {}
        for i, c in enumerate(L.CANCERS):
            ct = D[c]["counts_t"]
            lm = np.log(ct.div(ct.sum(1), axis=0).mean() + 1e-12)[ct.mean() >= 1]
            lm = lm.rename(dict(zip(clock_cols[c], CLOCK)))
            per_c[c] = [list(w.genes) for w in pseudo_waveforms(lm, a.n_sets, np.random.default_rng([SEED, k, i]), A.q.WF, kind)]
        sets += [(kind, j, {c: per_c[c][j] for c in L.CANCERS}) for j in range(a.n_sets)]
    # observed clock-gene COFE with permutation p and the clocksync agreement
    per = {c: cofe_pairs(c, clock_cols[c], SEED + i) for i, c in enumerate(L.CANCERS)}
    obs, null_mean, p = perm_p(per, 10_000, rng)
    rows = []
    for c, (pt, pn) in per.items():
        d = D[c]
        cols = clock_cols[c]
        case = d["counts_t"].loc[d["pairs"], cols].set_axis(CLOCK, axis=1).assign(depth=d["depth_t"].loc[d["pairs"]].to_numpy())
        ref = d["counts_n"].loc[d["pairs"], cols].set_axis(CLOCK, axis=1).assign(depth=d["depth_n"].loc[d["pairs"]].to_numpy())
        bt, _, _, _ = cohort_states(case, ref, True, A.q.WF)
        ang = np.arctan2(bt[:, 1] - bt[:, 1].mean(), bt[:, 0] - bt[:, 0].mean())
        R_same = abs(np.exp(1j * (pt - ang)).mean())
        R_flip = abs(np.exp(1j * (-pt - ang)).mean())
        rows.append(dict(cancer=c, pairs=len(pt), R_c=abs(np.exp(1j * (pt - pn)).mean()), cofe_vs_clocksync_R=max(R_same, R_flip)))
    pd.DataFrame(rows).to_csv(OUT / "a267_per_cancer.csv", index=False)
    log("clock-gene COFE pooled R", round(obs, 4), "p", p)
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        res = []
        for k, r in enumerate(pool.imap_unordered(one_set, sets[1:])):
            res.append(r)
            if (k + 1) % 20 == 0:
                log(f"{k + 1}/{len(sets) - 1} sets")
    S = pd.DataFrame([dict(kind="clock", set=0, R=obs), *res])
    S.to_csv(OUT / "a267_sets.csv", index=False)
    # generic axes on the A265 genome-wide COFE matrices
    axr = []
    pcs = {}
    for c in L.CANCERS:
        lt, ln = L.cofe_matrices(D[c]["counts_t"], D[c]["counts_n"])
        mu, sd = lt.mean(), lt.std().replace(0, 1)
        X = ((lt - mu) / sd).to_numpy()
        Xn = ((ln.loc[D[c]["pairs"]] - mu) / sd).to_numpy()
        _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
        pcs[c] = (pd.DataFrame(X @ Vt[:10].T, index=lt.index).loc[D[c]["pairs"]].to_numpy(), Xn @ Vt[:10].T)
    for i in range(10):
        for j in range(i + 1, 10):
            per_ax = {c: (np.arctan2(t[:, j], t[:, i]), np.arctan2(n[:, j], n[:, i])) for c, (t, n) in pcs.items()}
            axr.append(dict(pc_a=i + 1, pc_b=j + 1, R=pooled_R(per_ax)))
    AX = pd.DataFrame(axr)
    AX.to_csv(OUT / "a267_axes.csv", index=False)
    genome = float(pd.read_csv(OUT / "real_tests.csv").cofe_R.iloc[0])
    # descriptive: genome-wide COFE without the clock and GTEx-rhythmic genes of the A265 program
    prog = set(A.program().index)
    per_np = {}
    for i, c in enumerate(L.CANCERS):
        lt, ln = L.cofe_matrices(D[c]["counts_t"], D[c]["counts_n"])
        keep = [g for g in lt.columns if g not in prog]
        r = A.run_cofe(lt[keep], seed=SEED + 500 + i, ncores=a.jobs, test=ln[keep].loc[D[c]["pairs"]])
        per_np[c] = (r["phase"].loc[D[c]["pairs"]].to_numpy() * W, r["test_phase"].loc[D[c]["pairs"]].to_numpy() * W)
        log("genome without program", c)
    np_R, np_null, np_p = perm_p(per_np, 10_000, rng)
    t = dict(clock_R=obs, clock_null_mean=null_mean, clock_p=p)
    spec = p < 0.05
    for kind in ("matched", "responsive"):
        x = S[S.kind == kind].R
        t |= {f"{kind}_n": int(x.notna().sum()), f"{kind}_q95": float(np.nanquantile(x, 0.95)),
              f"{kind}_median": float(np.nanmedian(x)), f"{kind}_pct_below_clock": float(100 * np.mean(x < obs))}
        spec &= obs > t[f"{kind}_q95"]
    t["H267_1_pass"] = bool(spec)
    t |= dict(genome_R=genome, axes_q95=float(AX.R.quantile(0.95)), axes_median=float(AX.R.median()),
              axes_max=float(AX.R.max()), genome_pct_above_axes=float(100 * np.mean(AX.R < genome)))
    t["H267_2_pass"] = bool(genome > t["axes_q95"])
    t |= dict(genome_without_program_R=np_R, genome_without_program_null_mean=np_null, genome_without_program_p=np_p)
    pd.DataFrame([t]).to_csv(OUT / "a267_tests.csv", index=False)
    print(pd.Series(t).to_string(), flush=True)


if __name__ == "__main__":
    main()
