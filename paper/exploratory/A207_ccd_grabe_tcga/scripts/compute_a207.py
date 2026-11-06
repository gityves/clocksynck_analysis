"""A207: ΔCCD (Shilts 2018) and the Grabe 2024 order-based procedure on the TCGA tumour-normal pairs.

Registered in paper/docs/registered_experiments.md (A207). 10 cancers of the A138 table (`load_tcga`, `pair_frames`).
  * ΔCCD per cancer = CCD(tumours) - CCD(normals), `common.benchmark.ccd.ccd` (12 genes, Spearman, Shilts reference).
  * Grabe: PCA of the normals' centred log expression (14 clock genes), phase = atan2(PC2, PC1); the paired tumour
    takes the phase of its normal; per gene cosinor R^2 on that phase in normals and in tumours; per cancer
    dR2 = median over genes (tumour - normal).
  * Null: tumour/normal labels swapped within patient (1,000 swaps, seed 207; the same swap masks for ΔCCD, ΔR^2 and
    the pseudo sets); pooled = mean over cancers weighted by n_pairs, the same swaps jointly.
  * Companion (H207.2): 200 expression-matched pseudo-clock sets per cancer (A52 rule, `paired.pseudo_waveforms`,
    genome counts via `_shared.tcga_full_counts`, as A138 `pseudo_pooled`), the same dR2 pooled over cancers by set
    index. H207.2 counts only if the clock pooled dR2 is below the 5th percentile of the pseudo-set pooled dR2.

Run: RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A207_ccd_grabe_tcga/scripts/compute_a207.py
Test: ... compute_a207.py --cancers BLCA UCEC --swaps 20 --pseudo 5 --out <dir>

Choices not in the registration:
  * Log expression = log2(count / library depth * 1e6 + 1); depth = the table's `depth` (clock genes) or the genome
    column sum (pseudo sets, as A138). CCD and Grabe use the same log expression.
  * PCA = SVD of the column-centred (not scaled) log matrix; PC scores = U*S.
  * Swap p one-sided, (1 + #{null >= obs}) / (B + 1) for ΔCCD and (1 + #{null <= obs}) / (B + 1) for ΔR^2.
  * Verdict: H207.k passes if its Holm-adjusted (over the two) p < 0.05; H207.2 additionally needs the companion.
  * Pseudo sets are drawn per cancer (rng [207, 44, i]); pooled by set index across cancers (A138 rule). Pseudo
    swap nulls use the same masks; per-set one-sided swap p is reported in pseudo.csv (descriptive).
  * Mean log fraction for pseudo matching from normals with mean count >= 1 (A138 rule).
  * Spearman across cancers (A138 corr vs dCCD / dR2) descriptive, rows 'desc_*' of tests.csv.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve()
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "paper"))
sys.path.insert(0, str(REPO / "paper/results/figure_5/scripts/cancer"))
import _shared as sh  # noqa: E402
from common.benchmark.ccd import ccd, reference  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

SEED = 207
GENES = sh.GENES
PER_UNIT = REPO / "paper/exploratory/A138_cancer_estimability/output/score/per_unit.csv"
OUT = HERE.parents[1] / "output"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def lcpm(counts: pd.DataFrame, depth) -> pd.DataFrame:
    return np.log2(counts.div(np.asarray(depth, float), axis=0) * 1e6 + 1)


def delta_r2(Xt: np.ndarray, Xn: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Median-over-genes (tumour - normal) cosinor R^2 for each swap mask (B, n); row of True = labels swapped.
    Xt, Xn: (n, G) log expression of the tumours and normals. Phase = PCA of the (swapped) normals."""
    M = mask[:, :, None]
    A = np.where(M, Xt, Xn)  # reference (normal) rows
    T = np.where(M, Xn, Xt)
    U, S, _ = np.linalg.svd(A - A.mean(1, keepdims=True), full_matrices=False)
    ph = np.arctan2(U[:, :, 1] * S[:, None, 1], U[:, :, 0] * S[:, None, 0])
    X = np.stack([np.ones_like(ph), np.cos(ph), np.sin(ph)], axis=2)  # (B, n, 3)
    H = X @ np.linalg.pinv(X)  # (B, n, n)

    def r2(Y):
        Yc = Y - Y.mean(1, keepdims=True)
        res = Y - H @ Y
        return 1 - (res ** 2).sum(1) / (Yc ** 2).sum(1)  # (B, G)

    return np.median(r2(T) - r2(A), axis=1)


def delta_ccd(Xt: pd.DataFrame, Xn: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    ref_genes = [g for g in reference().index if g in Xt.columns]
    t, n = Xt[ref_genes], Xn[ref_genes]
    out = np.empty(len(mask))
    for b, m in enumerate(mask):
        a = pd.DataFrame(np.where(m[:, None], t.to_numpy(), n.to_numpy()), columns=ref_genes)  # normals'
        c = pd.DataFrame(np.where(m[:, None], n.to_numpy(), t.to_numpy()), columns=ref_genes)  # tumours'
        out[b] = ccd(c)[0] - ccd(a)[0]
    return out


def pool(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted mean over axis 0 (cancers)."""
    return (x * w.reshape(-1, *[1] * (x.ndim - 1))).sum(0) / w.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancers", nargs="*", default=sh.CANCERS)
    ap.add_argument("--swaps", type=int, default=1000)
    ap.add_argument("--pseudo", type=int, default=200)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    B, P = a.swaps, a.pseudo

    samples = sh.load_tcga()
    a138 = pd.read_csv(PER_UNIT).query("hypothesis == 'H138.4'").set_index("unit")
    rows, n_pairs = [], []
    dccd, dr2, pdr2 = [], [], []  # (cancer, B+1), (cancer, B+1), (cancer, P, B+1)
    for i, c in enumerate(a.cancers):
        case, ref = sh.pair_frames(samples, c)
        n = len(case)
        mask = np.vstack([np.zeros(n, bool), np.random.default_rng([SEED, i]).random((B, n)) < 0.5])
        Xt = lcpm(case[GENES], case.depth)
        Xn = lcpm(ref[GENES], ref.depth)
        dd = delta_ccd(Xt, Xn, mask)
        dr = delta_r2(Xt.to_numpy(), Xn.to_numpy(), mask)
        ct, cn = sh.tcga_full_counts(c, case.index)
        dt, dn = ct.sum(axis=1), cn.sum(axis=1)
        lm = np.log(cn.div(dn, axis=0).mean() + 1e-12)[cn.mean() >= 1]
        pr = []
        for w in pseudo_waveforms(lm, P, np.random.default_rng([SEED, 44, i]), sh.WF):
            g = list(w.genes)
            pr.append(delta_r2(lcpm(ct[g], dt).to_numpy(), lcpm(cn[g], dn).to_numpy(), mask))
        dccd.append(dd); dr2.append(dr); pdr2.append(np.array(pr)); n_pairs.append(n)
        cc, cr = ccd(Xt[[g for g in reference().index]])[0], ccd(Xn[[g for g in reference().index]])[0]
        rows.append(dict(cancer=c, n_pairs=n, ccd_normal=cr, ccd_tumour=cc, delta_ccd=dd[0], dR2=dr[0],
                         a138_corr=a138.loc[c, "statistic"], a138_p=a138.loc[c, "p"]))
        log(c, n, f"dCCD {dd[0]:+.4f} dR2 {dr[0]:+.4f}")
    w = np.array(n_pairs, float)
    per = pd.DataFrame(rows)
    per.to_csv(a.out / "per_cancer.csv", index=False)

    pc, pr_, pp = pool(np.array(dccd), w), pool(np.array(dr2), w), pool(np.array(pdr2), w)  # (B+1), (B+1), (P, B+1)
    p1 = (1 + (pc[1:] >= pc[0]).sum()) / (B + 1)
    p2 = (1 + (pr_[1:] <= pr_[0]).sum()) / (B + 1)
    pseudo_obs = pp[:, 0]
    pct = float((pseudo_obs <= pr_[0]).mean())  # share of pseudo sets with dR2 <= clock dR2
    q05 = float(np.quantile(pseudo_obs, 0.05))
    comp = bool(pr_[0] < q05)
    order = np.argsort([p1, p2])
    holm = np.empty(2); raw = np.array([p1, p2])
    holm[order] = np.minimum(1, np.maximum.accumulate((2 - np.arange(2)) * raw[order]))
    t = pd.DataFrame([
        dict(hypothesis="H207.1", statistic="pooled dCCD", observed=pc[0], null_mean=pc[1:].mean(), p_raw=p1, p_holm=holm[0],
             companion_percentile=np.nan, companion_q05=np.nan, verdict="pass" if holm[0] < 0.05 else "fail"),
        dict(hypothesis="H207.2", statistic="pooled dR2", observed=pr_[0], null_mean=pr_[1:].mean(), p_raw=p2, p_holm=holm[1],
             companion_percentile=pct, companion_q05=q05,
             verdict="pass" if (holm[1] < 0.05 and comp) else ("fail (companion)" if holm[1] < 0.05 else "fail")),
    ])
    for nm, x, y in (("desc_spearman_dCCD_vs_a138corr", per.delta_ccd, per.a138_corr),
                     ("desc_spearman_dR2_vs_a138corr", per.dR2, per.a138_corr)):
        r = spearmanr(x, y)
        t = pd.concat([t, pd.DataFrame([dict(hypothesis=nm, statistic="Spearman rho (n=%d cancers)" % len(per),
                                             observed=r.statistic, p_raw=r.pvalue, verdict="descriptive")])], ignore_index=True)
    t.to_csv(a.out / "tests.csv", index=False)
    pd.DataFrame(dict(set=np.arange(P), pooled_dR2=pseudo_obs,
                      p_swap_lower=(1 + (pp[:, 1:] <= pp[:, :1]).sum(1)) / (B + 1))).to_csv(a.out / "pseudo.csv", index=False)
    log("done", t.to_string())


if __name__ == "__main__":
    main()
