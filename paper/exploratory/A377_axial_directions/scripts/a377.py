"""A377 (registered 2026-11-25): is the brain-death loss of the organ arrangement a displacement along tissue-specific
axes rather than a change of the relations between organ clocks?

Axial factor: rank-1 Takagi fit of the pseudo-covariance P = E[z_a z_b] on the observed off-diagonal pairs (diagonal
imputed with the off-diagonal pairs, as registered: tissue-own non-circular noise does not enter). Axis alignment
A = |sum_a w_a^2| / sum_a |w_a|^2 (1: one common direction); axial size: pooled pseudo-covariance ratio (A109
`circ_ratio`).

Stages:
  power  V377.1 gate on the A133 simulator (ventilated design, fast-death offsets, n 79, 200 replicates per cell, f 0.2
         and 0.44, common vs tissue-specific directions); SD of D_A = A(post-mortem) - A(organ donors) and of the
         circularity-ratio difference under the A317 relabelling scheme (2000 relabellings of Hardy 0 donors into two
         groups of 79). No group-split observed value is computed.
  score  (added after the power note)

    python a377.py --stage power
Outputs (`../output/`): gate.csv (cells), gate_reps.csv, power.csv. Aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A317_braindeath_confirm" / "scripts"))
sys.path.insert(0, str(EXP / "A133_axial_decomposition" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import compute_v133_sims as S  # noqa: E402
import score_a317 as A  # noqa: E402
from axial import takagi_rank1  # noqa: E402
from score_a317 import B, H, gtex, raw  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
N_REP, N_PERM, SEED = 200, 2000, 377


def fit_offdiag(P, observed, n_iter=500, tol=1e-10):
    """Rank-1 Takagi factor w (P ~ w w^T) fitted on the observed off-diagonal entries; every other entry, the diagonal
    included, imputed from the current fit (A133 `fit_axial` with the diagonal treated as unobserved)."""
    obs = observed | observed.T
    np.fill_diagonal(obs, False)
    P = (P + P.T) / 2
    w = takagi_rank1(np.where(obs, P, 0))
    for _ in range(n_iter):
        w_new = takagi_rank1(np.where(obs, P, np.outer(w, w)))
        done = np.linalg.norm(np.outer(w_new, w_new) - np.outer(w, w)) <= tol * max(np.linalg.norm(np.outer(w, w)), 1e-300)
        w = w_new
        if done:
            break
    return w


def alignment(w):
    return float(abs(np.sum(w ** 2)) / max(np.sum(abs(w) ** 2), 1e-300))


def circ_ratio(C, P, observed):
    iu = np.triu_indices(len(C), 1)
    o = observed[iu]
    return float(np.mean(abs(P[iu][o])) / np.mean(abs(C[iu][o])))


def readout(Z, M, observed):
    C, P, _ = H.cross_moments(Z, M)
    return alignment(fit_offdiag(P, observed)), circ_ratio(C, P, observed)


def centred(Z, M, idx):
    m = M[idx]
    mu = (Z[idx] * m).sum(0) / np.maximum(m.sum(0), 1)
    return (Z[idx] - mu) * m, m


def gate():
    a, V, phi = S.design()
    rng = np.random.default_rng(SEED)
    rows = []
    for f in (0.2, 0.44):
        for u in ("common", "pertissue"):
            for r in range(N_REP):
                Z, M = S.simulate(79, f, 0.0, u, a, V, phi, rng)
                _, _, N = H.cross_moments(Z, M)
                al, cr = readout(Z, M, N >= S.MIN_N)
                rows.append(dict(f=f, directions=u, rep=r, A=al, circ_ratio=cr))
    reps = pd.DataFrame(rows)
    cells = reps.groupby(["f", "directions"]).agg(A_median=("A", "median"), A_q05=("A", lambda x: x.quantile(0.05)),
                                                 A_q95=("A", lambda x: x.quantile(0.95)),
                                                 circ_ratio_median=("circ_ratio", "median")).reset_index()
    auc = {f: roc_auc_score((reps[reps.f == f].directions == "common").astype(int), reps[reps.f == f].A)
           for f in (0.2, 0.44)}
    cells["auc"] = cells.f.map(auc)
    c44 = cells[cells.f == 0.44].set_index("directions")
    passes = bool(c44.A_median["common"] >= 0.8 and c44.A_median["pertissue"] <= 0.5 and auc[0.44] >= 0.9)
    cells["V377_1_passes"] = passes
    return reps, cells


def power():
    limit_threads(4)
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(A.SEED)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    Z, M = raw(T, "clock", 0, T.R)
    mask = B.pair_mask(M, fast, pm, *od_draws[:1]) & B.pair_mask(M, od)  # the A317 pairs
    obs = np.zeros((len(H.BODY), len(H.BODY)), bool)
    for (a, b), m in zip(H.PAIRS, mask):
        obs[a, b] = obs[b, a] = m
    rng = np.random.default_rng(SEED)
    d_a, d_c = [], []
    for _ in range(N_PERM):  # relabelled Hardy 0 donors only: no group-split value
        p = rng.permutation(h0)
        ra, rc = readout(*centred(Z, M, p[:len(pm)]), obs)
        oa, oc = readout(*centred(Z, M, p[len(pm):2 * len(pm)]), obs)
        d_a.append(ra - oa)
        d_c.append(rc - oc)
    d_a, d_c = np.array(d_a), np.array(d_c)
    return pd.Series(dict(pairs=int(mask.sum()), n_group=len(pm), sd_D_A=float(np.std(d_a, ddof=1)),
                          mde_D_A=float(2.49 * np.std(d_a, ddof=1)), q95_D_A=float(np.quantile(d_a, 0.95)),
                          sd_D_circ=float(np.std(d_c, ddof=1)), mde_D_circ=float(2.49 * np.std(d_c, ddof=1))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power"], required=True)
    ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    reps, cells = gate()
    reps.to_csv(OUT / "gate_reps.csv", index=False)
    cells.to_csv(OUT / "gate.csv", index=False)
    print(cells.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)", flush=True)
    pw = power()
    pw.to_csv(OUT / "power.csv", header=["value"])
    print(pw.round(4).to_string(), f"\n({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
