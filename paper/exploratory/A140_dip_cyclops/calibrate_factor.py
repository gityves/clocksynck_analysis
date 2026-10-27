"""A140 Q1 revision: calibration of the non-clock patient factor from the genome-wide tumour immune pseudobulks.

Registered (paper/docs/registered_experiments.md, A140 "Q1 revision"): per compartment (tumour T_NK, tumour myeloid) and
cohort, log CPM centred within cohort; the projection on the 14 clock genes' host plane removed; the expected Poisson
sampling variance of log CPM subtracted per gene; PCA; median over cohorts = factor variance at scale 1.

Concrete choices (README, "Choices not in the registration"):
- log CPM = ln((count + 0.5) / depth), depth = the pseudobulk's genome-wide UMI (natural log: the unit of the simulator's
  log-fold factor).
- clock genes are found under their old symbols too (ARNTL = BMAL1, CIART = C1orf51); cohorts with < 5 samples in a
  compartment are skipped (SRP273482: 7 T_NK and no myeloid; SRP335254: none in the main table).
- genes: mean count >= 5 per pseudobulk in both compartments of the cohort (the sampling variance is then well defined).
- host plane of a compartment = the span of the first two sample-space principal components of its 14 clock genes
  (centred ln-CPM); each gene's vector over the cohort's samples is made orthogonal to it (and to the constant).
- Poisson variance of a gene in a pseudobulk of depth N: Var ln((X + 0.5) / N), X ~ Poisson(N p_g), p_g = the compartment's
  pooled fraction in the cohort (Monte Carlo, 200 draws). After centring and the 2-D projection a compartment of n samples
  keeps n - 3 degrees of freedom of iid noise, so the diagonal subtracted from Y'Y is sum_c (n_c - 3) mean_j v_gjc.
- PCA on the two compartments' rows stacked (each compartment centred and projected separately): eigenvalues of
  (Y'Y - D) / (n_rows - 1) by Lanczos (genes x genes, never formed); "before" = the same without D.

Output (git-ignored): output/q1/tumour_immune_genome/calibration_cohorts.csv (top 5 eigenvalues before/after per cohort),
calibration_table.csv (median over cohorts) and calibration_lambda.csv (lam used by the simulator at scale 1).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
from build_tumour_immune_genome import COMPS, OUT, targets  # noqa: E402

CLOCK = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"]
MIN_MEAN_COUNT, N_MC, NEIG = 5.0, 200, 5


def poisson_var(mu: np.ndarray, N: np.ndarray, rng) -> np.ndarray:
    """Var of ln((X+0.5)/N), X ~ Poisson(mu), per row (depth N_j) and gene; mu: rows x genes. Monte Carlo."""
    s1 = np.zeros_like(mu)
    s2 = np.zeros_like(mu)
    for _ in range(N_MC):
        y = np.log((rng.poisson(mu) + 0.5) / N[:, None])
        s1 += y
        s2 += y * y
    return np.maximum(s2 / N_MC - (s1 / N_MC) ** 2, 0) * N_MC / (N_MC - 1)


def cohort(study: str, rng):
    X = sp.load_npz(OUT / f"{study}.npz").tocsr()
    rows = pd.read_csv(OUT / f"{study}_rows.csv")
    genes = np.array((OUT / f"{study}_genes.txt").read_text().split("\n"))
    blocks, dvec, n_rows = [], None, 0
    parts = {}
    for c in COMPS:
        i = np.flatnonzero((rows.compartment == c).to_numpy())
        if len(i) < 5:
            return None
        parts[c] = (X[i].toarray(), rows.depth.to_numpy()[i])
    keep = np.all([parts[c][0].mean(0) >= MIN_MEAN_COUNT for c in COMPS], axis=0)
    alias = {"ARNTL": ("ARNTL", "BMAL1"), "CIART": ("CIART", "C1orf51")}
    ci = [int(np.flatnonzero(np.isin(genes, alias.get(g, (g,))))[0]) for g in CLOCK if np.isin(genes, alias.get(g, (g,))).any()]
    if len(ci) < 14:
        return None
    for c in COMPS:
        C, N = parts[c]
        n = len(N)
        Y = np.log((C[:, keep] + 0.5) / N[:, None])
        Yc = np.log((C[:, ci] + 0.5) / N[:, None])
        Y -= Y.mean(0)
        Yc -= Yc.mean(0)
        u = np.linalg.svd(Yc, full_matrices=False)[0][:, :2]  # host plane (sample space)
        Y -= u @ (u.T @ Y)
        p = C[:, keep].sum(0) / N.sum()
        v = poisson_var(N[:, None] * p[None], N, rng).mean(0)
        blocks.append(Y)
        dvec = (n - 3) * v if dvec is None else dvec + (n - 3) * v
        n_rows += n
    Y = np.vstack(blocks)
    return Y, dvec, n_rows, int(keep.sum()), [len(parts[c][1]) for c in COMPS]


def eig(Y, d, n_rows):
    k = min(NEIG, Y.shape[0] - 6)
    G = Y.shape[1]
    op = LinearOperator((G, G), matvec=lambda v: (Y.T @ (Y @ v) - d * v) / (n_rows - 1), dtype=float)
    return np.sort(eigsh(op, k=k, which="LA", return_eigenvectors=False))[::-1]


def main():
    rng = np.random.default_rng(0)
    out = []
    for s in targets():
        r = cohort(s, rng)
        if r is None:
            print(s, "skipped", flush=True)
            continue
        Y, d, nr, ng, ns = r
        before = eig(Y, np.zeros_like(d), nr)
        after = eig(Y, d, nr)
        row = dict(study=s, n_genes=ng, n_tnk=ns[0], n_myeloid=ns[1], total_var_before=float((Y ** 2).sum() / (nr - 1)),
                   poisson_total=float(d.sum() / (nr - 1)))
        for k in range(NEIG):
            row[f"before_{k + 1}"] = before[k] if k < len(before) else np.nan
            row[f"after_{k + 1}"] = after[k] if k < len(after) else np.nan
        out.append(row)
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    t = pd.DataFrame(out)
    t.to_csv(OUT / "calibration_cohorts.csv", index=False)
    med = t.drop(columns="study").median()
    med.rename("median_over_cohorts").to_frame().to_csv(OUT / "calibration_table.csv")
    lam = pd.DataFrame({"k": range(1, 4), "lam": [med[f"after_{k}"] for k in range(1, 4)],
                        "lam_before_poisson": [med[f"before_{k}"] for k in range(1, 4)]})
    lam.to_csv(OUT / "calibration_lambda.csv", index=False)
    print(f"{len(t)} cohorts\n", lam.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
