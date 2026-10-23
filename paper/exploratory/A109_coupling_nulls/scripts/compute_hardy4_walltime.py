"""A109c (secondary, own alpha 0.05, outside the confirmatory family): is the Hardy 4 result inherited from wall time?

Registered in paper/docs/registered_experiments.md, section A109c (H109c.1, H109c.2; D109c.3/4 are elsewhere).
  H109c.1  R, the A92/A93 wall-time alignment (per tissue |mean exp(i (phase of the clock state - w t_death))|, states
           centred within Hardy 4, median over the 12 tissues, tissues with >= 20 donors; equal n 103 = all Hardy 4
           donors), against 1,000 permutations of the time of death among Hardy 4 donors; one-sided add-one p.
           MDE: (i) null q95 of R; (ii) R of the simulated alignment (angle = w t + von Mises noise of concentration
           kappa, observed donors and tissues) that gives p < 0.05 in 80 % of simulations.
  H109c.2  H109.2 in Hardy 4 after regressing every tissue's state on (cos, sin) of time of death within Hardy 4;
           fast-death pair phases (unresidualised) as reference; Fisher-Lee over 66 pairs, 10,000 joint tissue-label
           permutations; p_closer by donor-label permutation between fast and Hardy 4 (residualised); pass: rho > 0,
           p < 0.05 and p_closer > 0.05. Rows: H4_raw (A109 secondary), H4_residual, and the residual variant without the
           Hardy 4 ventilator cases (DTHVNT = 1) and in hospital inpatients (DTHPLCE).
Controlled access: time of death and phenotypes in memory only; aggregates in the outputs.

    python compute_hardy4_walltime.py [--subsample N] [--seed S]
Outputs (paper/exploratory/A109_coupling_nulls/output/): a109c_walltime.csv, a109c_offsets_residual.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_coupling_nulls as N  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

W = 2 * np.pi / 24
BODY = H.BODY


def tissue_arrays(T, idx):
    """Per tissue: (centred complex states, time of death) of the donors idx present in it."""
    out = []
    for t in BODY:
        r = T.R[t][idx]
        m = r >= 0
        if m.sum() < H.MIN_N:
            continue
        z = T.Z["clock"][t][0][r[m]]
        out.append((z - z.mean(), T.tod[idx][m]))
    return out


def align_R(arr, tod_perm=None):
    r = []
    for z, tod in arr:
        r.append(abs(np.exp(1j * (np.angle(z) - W * tod)).mean()))
    return float(np.median(r))


def walltime(T, rng, sz):
    idx = T.group_donors("H4")
    arr = tissue_arrays(T, idx)
    R = align_R(arr)
    tod_d = T.tod[idx]
    pos = [np.flatnonzero(T.R[t][idx] >= 0) for t in BODY if (T.R[t][idx] >= 0).sum() >= H.MIN_N]
    ang = [np.angle(z) for z, _ in arr]

    def R_of(todv):  # todv: tod of the donors idx (permuted among donors)
        return float(np.median([abs(np.exp(1j * (a - W * todv[p])).mean()) for a, p in zip(ang, pos)]))

    null = np.array([R_of(rng.permutation(tod_d)) for _ in range(sz["n_perm"])])
    q95 = float(np.quantile(null, 0.95))
    rows = []
    for kappa in (0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5):
        rs = []
        for _ in range(sz["n_sim"]):
            vm = rng.vonmises(0, kappa, len(idx))
            a_sim = [W * tod_d[p] + vm[p] for p in pos]
            rs.append(float(np.median([abs(np.exp(1j * (a - W * tod_d[p])).mean()) for a, p in zip(a_sim, pos)])))
        rs = np.array(rs)
        rows.append(dict(kappa=kappa, R_median=np.median(rs), power=float(np.mean(rs > q95))))
    sim = pd.DataFrame(rows)
    ok = sim[sim.power >= 0.8]
    mde = float(np.interp(0.8, sim.power, sim.R_median)) if len(ok) else np.nan
    out = pd.DataFrame([dict(group="H4", R=R, null_median=np.median(null), null_sd=null.std(ddof=1), null_q95=q95,
                             p_one_sided=float((1 + np.sum(null >= R)) / (1 + len(null))), n_donors=len(idx),
                             n_tissues=len(arr), n_perm=len(null), mde_null_q95=q95, mde_80pct_power_R=mde)])
    return out, sim


def residualise(Zraw, M, idx, tod):
    """Zraw rows idx regressed, per tissue, on (1, cos, sin) of the time of death (donors present in the tissue)."""
    Z = Zraw.copy()
    w = W * tod
    for j in range(Zraw.shape[1]):
        m = M[idx, j] > 0
        if m.sum() < 4:
            continue
        X = np.c_[np.ones(m.sum()), np.cos(w[m]), np.sin(w[m])]
        Y = np.c_[Zraw[idx[m], j].real, Zraw[idx[m], j].imag]
        e = Y - X @ np.linalg.lstsq(X, Y, rcond=None)[0]
        Z[idx[m], j] = e[:, 0] + 1j * e[:, 1]
    return Z


def offsets_row(Zraw, M, fast, gd, rng, sz, label):
    th_f = H.pair_phase(N._H_of(Zraw, M, fast))
    Hg = N._H_of(Zraw, M, gd)
    th_g = H.pair_phase(Hg)
    r = ST.fisher_lee(th_f, th_g)
    null = np.array([ST.fisher_lee(th_f, H.permute_labels(Hg, rng.permutation(len(BODY)))) for _ in range(sz["n_fl"])])
    d_f, d_g = H.folded_distance(th_f), H.folded_distance(th_g)
    pool, nf = np.concatenate([fast, gd]), len(fast)
    nl = []
    for _ in range(sz["n_closer"]):
        p = rng.permutation(pool)
        nl.append(H.folded_distance(H.pair_phase(N._H_of(Zraw, M, p[:nf])))
                  - H.folded_distance(H.pair_phase(N._H_of(Zraw, M, p[nf:]))))
    p_rho, p_cl = N.pval(r, null), N.pval(d_f - d_g, nl)
    return dict(variant=label, rho_fisher_lee=r, p=p_rho, null_median=np.median(null),
                median_folded_phase_h=H.folded_median(th_g), folded_distance_h=d_g, folded_distance_fast_h=d_f,
                p_closer=p_cl, pass_rule=bool(r > 0 and p_rho < 0.05 and p_cl > 0.05), n_pairs=len(th_g),
                n_donors=len(gd), n_fast=nf, n_perm=len(null))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    q = run["quick"]
    sz = dict(n_perm=100 if q else 1000, n_sim=30 if q else 200, n_fl=200 if q else 10_000,
              n_closer=100 if q else 1000)
    out = exploratory_dir(N.FOLDER, q)
    rng = np.random.default_rng(1209 + args.seed)
    T = H.Tissues(1, args.seed, args.subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    wt, sim = walltime(T, rng, sz)
    wt.to_csv(out / "a109c_walltime.csv", index=False)
    sim.to_csv(out / "a109c_walltime_power.csv", index=False)
    print(wt.round(4).to_string(index=False), "\n", sim.round(3).to_string(index=False), flush=True)

    ph = gtex.phenotypes(["DTHVNT", "DTHPLCE"]).reindex(T.donors)
    Zraw, M, _ = N._raw_matrix(T)
    fast, h4 = T.group_donors("fast"), T.group_donors("H4")
    vent = (ph.DTHVNT.to_numpy() == 1)
    inpt = (ph.DTHPLCE.to_numpy() == "Hospital inpatient")
    rows = [offsets_row(Zraw, M, fast, h4, rng, sz, "H4_raw")]
    for label, sel in (("H4_residual", h4), ("H4_residual_no_ventilator", h4[~vent[h4]]),
                       ("H4_residual_inpatients", h4[inpt[h4]])):
        Zr = residualise(Zraw, M, sel, T.tod[sel])
        rows.append(offsets_row(Zr, M, fast, sel, rng, sz, label))
    d = pd.DataFrame(rows)
    d.to_csv(out / "a109c_offsets_residual.csv", index=False)
    print(d.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
