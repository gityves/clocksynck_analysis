"""A307: is the fetal clock-gene covariation the glucocorticoid-induction pattern or the waveform's phase relations?

Registered in paper/docs/registered_experiments.md (A307; NOT BLIND). Rows and adjustment as A264 (adrenal cortex) and
A269 (lung, liver): log fractions residualised within organ on cohort indicators and common quadratic age.
  L  GC link: mean over SI = (PER1, PER2, NR1D1) of r^2 with the GC score (mean of the standardised adjusted FKBP5,
     TSC22D3, ZBTB16, KLF9). Null: GC score permuted within cohort. Companion: 200 expression-matched non-clock 3-gene
     sets in place of SI (sets with a GC target dropped).
  W  waveform cross-agreement: mean over (SI, CORE) pairs of r_gh cos(phi_g - phi_h). Null: every gene permuted on its
     own within cohort. Companion: 200 expression-matched sets on the full waveform, same positions.
Holm over the three organs per statistic; pass = Holm p < 0.05 and above the matched q95.

    python score_a307.py --stage power|score [--n-perm 10000] [--n-pseudo 200]
Outputs (../output/): power.csv (null only), tests.csv, genes.csv (signed r of every clock gene with the GC score),
companions.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
for p in ("A269_fetal_organs_drive", "A266_fetal_adrenal_drive", "A264_fetal_adrenal_covariates", "A244_fetal_adrenal_classes"):
    sys.path.insert(0, str(EXP / p / "scripts"))
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import design_a264 as d264  # noqa: E402
import design_a269 as d269  # noqa: E402
import prep_a244 as PR  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import development as dev  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 307
WF = d269.WF
GENES = list(WF.genes)
SI = list(d269.SIGNAL_INDUCED)
CORE = [g for g in GENES if g not in SI]
GC = list(d269.GC)
ORGANS = ("cortex", "lung", "liver")
PHI = dict(zip(GENES, WF.psi))


def load():
    """{organ: (counts of the rows (all genes), obs rows, covariate function)}; read before the thread limit."""
    out = {}
    obs = d264.d244.sheet()
    rows = d264.d244.pooled(obs, obs.compartment == "cortex")
    pb = pd.read_parquet(PR.PB)
    out["cortex"] = (pb.loc[obs[rows].index], obs[rows], d264.cov())
    del pb
    o = d269.d243.sheet()
    lp = pd.read_parquet(d269.d243.SRC / "pseudobulk.parquet")
    for org in ("lung", "liver"):
        r = d269.d243.pooled(o, o.organ == org)
        out[org] = (lp.loc[o[r].index], o[r], d269.cov())
    return out


def adjusted(counts, obs, cov, genes) -> tuple[np.ndarray, list]:
    """Centred adjusted log fractions (rows x genes) and the cohort row groups."""
    tab, _ = dev.adjust(dev.table(counts, obs, genes), list(genes), cov)
    y = np.log(tab[list(genes)].to_numpy(float) / tab.depth.to_numpy(float)[:, None])
    return y - y.mean(0), list(tab.groupby("block", sort=False).indices.values())


def gc_score(ygc) -> np.ndarray:
    return ((ygc - ygc.mean(0)) / ygc.std(0)).mean(1)


def corr_with(Y, g) -> np.ndarray:
    Yc, gc = Y - Y.mean(0), g - g.mean()
    return (Yc * gc[:, None]).sum(0) / np.sqrt((Yc ** 2).sum(0) * (gc ** 2).sum())


def L_stat(Ysi, g) -> float:
    return float(np.mean(corr_with(Ysi, g) ** 2))


def W_stat(Y, si_idx, core_idx, cosm) -> float:
    R = np.corrcoef(Y, rowvar=False)
    return float(np.mean(R[np.ix_(si_idx, core_idx)] * cosm))


def perm_vec(g, groups, rng):
    out = g.copy()
    for i in groups:
        out[i] = g[rng.permutation(i)]
    return out


def summary(null) -> dict:
    x = np.asarray(null)[np.isfinite(null)]
    return dict(null_mean=float(x.mean()), null_sd=float(x.std(ddof=1)), null_q95=float(np.quantile(x, .95)),
                mde=float(np.quantile(x, .95) + 0.84 * x.std(ddof=1)))


def nulls(y, groups, si_idx, core_idx, cosm, g, n_perm, oi):
    rng = np.random.default_rng([SEED, oi, 0])
    nL = np.array([L_stat(y[:, si_idx], perm_vec(g, groups, rng)) for _ in range(n_perm)])
    rng = np.random.default_rng([SEED, oi, 2])
    nW = np.array([W_stat(dev._permuted(y, groups, rng), si_idx, core_idx, cosm) for _ in range(n_perm)])
    return nL, nW


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    a = ap.parse_args()
    t0 = time.time()
    data = load()
    from clocksync.config import limit_threads
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    si_idx = [GENES.index(x) for x in SI]
    core_idx = [GENES.index(x) for x in CORE]
    cosm = np.cos(np.subtract.outer([PHI[x] for x in SI], [PHI[x] for x in CORE]))
    prow, trow, grow, crow = [], [], [], []
    for oi, org in enumerate(ORGANS):
        counts, obs, cov = data[org]
        y, groups = adjusted(counts, obs, cov, GENES + GC)
        Y, g = y[:, :len(GENES)], gc_score(y[:, len(GENES):])
        nL, nW = nulls(Y, groups, si_idx, core_idx, cosm, g, a.n_perm, oi)
        prow += [dict(organ=org, statistic=s, n=len(obs), cohorts=len(groups), **summary(n)) for s, n in (("L", nL), ("W", nW))]
        print(org, "nulls done", f"{time.time() - t0:.0f}s", flush=True)
        if a.stage == "power":
            continue
        L, W = L_stat(Y[:, si_idx], g), W_stat(Y, si_idx, core_idx, cosm)
        r = corr_with(Y, g)
        grow += [dict(organ=org, gene=x, set="SI" if x in SI else "CORE", r_gc=float(v)) for x, v in zip(GENES, r)]
        lm = dev.log_mean(counts, obs.depth)
        # companions of L: matched 3-gene sets in place of SI (sets containing a GC target dropped)
        ws = pseudo_waveforms(lm, a.n_pseudo + 100, np.random.default_rng([SEED, oi, 1]), WF.drop(CORE))
        ws = [w for w in ws if not set(w.genes) & set(GC)][:a.n_pseudo]
        cL = []
        for w in ws:
            yw, _ = adjusted(counts, obs, cov, list(w.genes))
            cL.append(L_stat(yw, g))
        # companions of W: matched sets on the full waveform, same positions
        cW = []
        for w in pseudo_waveforms(lm, a.n_pseudo, np.random.default_rng([SEED, oi, 3]), WF):
            yw, _ = adjusted(counts, obs, cov, list(w.genes))
            cW.append(W_stat(yw, si_idx, core_idx, cosm))
        crow += [dict(organ=org, statistic="L", set=j, value=v) for j, v in enumerate(cL)]
        crow += [dict(organ=org, statistic="W", set=j, value=v) for j, v in enumerate(cW)]
        for s, o, n, c in (("L", L, nL, cL), ("W", W, nW, cW)):
            trow.append(dict(organ=org, statistic=s, value=o, p=float((1 + np.sum(n >= o)) / (1 + len(n))),
                             null_q95=float(np.quantile(n, .95)), matched_median=float(np.median(c)),
                             matched_q95=float(np.quantile(c, .95)), matched_pct=float(100 * np.mean(np.array(c) < o)),
                             n_sets=len(c)))
        print(org, "scored", f"{time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(prow).to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(prow).round(4).to_string(index=False), flush=True)
    if a.stage == "power":
        return
    T = pd.DataFrame(trow)
    for s in ("L", "W"):
        m = T.statistic == s
        T.loc[m, "p_holm"] = dev.holm(T.p[m].to_numpy())
    T["clock_specific"] = T.value > T.matched_q95
    T["pass"] = (T.p_holm < 0.05) & T.clock_specific
    T.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(grow).to_csv(OUT / "genes.csv", index=False)
    pd.DataFrame(crow).to_csv(OUT / "companions.csv", index=False)
    print(T.round(4).to_string(index=False), pd.DataFrame(grow).pivot(index="gene", columns="organ", values="r_gc").round(3)
          .to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
