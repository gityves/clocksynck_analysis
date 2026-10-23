"""V14 (amendment A19, part 3a; paper-only, no clocksync call): the ORPHEUS phase-dispersion estimator when subject phases carry error.

ORPHEUS (Hollis et al. 2026, bioRxiv 10.64898/2026.03.23.713759) reads the phase dispersion between cells, sigma_c,
from the rhythms of the mean and of the variance between cells over subject phases (CYCLOPS phases in human brain).
Model per gene: cell mean mu = M + A cos(t + phi - psi), phi ~ N(0, sigma_c^2); counts ~ NB(mu, theta)
(var = mu + theta mu^2). Over cells at subject phase t (u = t - psi):
    E[count]   = M + A e^{-s/2} cos u,                                   s = sigma_c^2
    Var[count] = DC + A e^{-s/2} (1 + 2 theta M) cos u
                 + [(1 + theta) A^2 / 2 (e^{-2s} - e^{-s}) + theta A^2 e^{-s} / 2] cos 2u.
Moment solution (with A_m = A e^{-s/2}, B = 24-h variance coefficient, C = 12-h coefficient):
    theta = (B / A_m - 1) / (2 M),    e^{-s} = 1 + (2 C / A_m^2 - theta) / (1 + theta).
If the subject phase used for t has an error of SD sigma_s, the 24-h coefficients shrink by e^{-sigma_s^2/2} and the
12-h coefficient by e^{-2 sigma_s^2}. theta is unaffected (a ratio of two 24-h terms), but 2C/A_m^2 shrinks by
e^{-sigma_s^2}. When the negative-binomial term dominates the 12-h rhythm (2C/A_m^2 > 0), the estimated dispersion
then rises with sigma_s: subject-phase error, or a weaker coupling of subjects to the ordering, reads as
desynchrony between cells.

Simulation: per setting, `--subjects` subjects at phases uniform on the day, `--cells` cells each, `--genes` genes
(A/M in [0.15, 0.75], M in [0.5, 5], theta in [0.05, 1.5], as the ORPHEUS liver proxy), cell-specific additive
noise N(0, 0.3), assigned phase = true + N(0, sigma_s^2). Per gene: cosinor of the subject means on the assigned
phase; the subject variances regressed on cos/sin u and cos/sin 2u (u relative to the fitted mean acrophase); the
moment solution above; genes with a non-real solution are dropped (as ORPHEUS keeps only solvable genes). Reported:
median estimated sigma_c (h) over genes, per (sigma_c, sigma_s); and the analytic prediction.

    python compute_orpheus_bias.py [--subsample N] [--subjects 400 --cells 300 --genes 60 --reps 5]
`--subsample N`: check mode (N subjects, 40 cells, 10 genes, 1 replicate).
Outputs: paper/results/extended_data/output/calibration/orpheus_bias_{simulation,summary,analytic}.csv
Legacy: legacy/scripts/v14_orpheus_bias.py (outputs/v14_orpheus_bias/). Same generator and seed: identical numbers.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "calibration"
H = 24 / (2 * np.pi)


def simulate(sigma_c_h, sigma_s_h, n_sub, n_cells, n_genes, rng):
    sc, ss = sigma_c_h / H, sigma_s_h / H
    t = rng.uniform(0, 2 * np.pi, n_sub)
    t_hat = t + rng.normal(0, ss, n_sub) if ss > 0 else t
    rows = []
    for g in range(n_genes):
        M = rng.uniform(0.5, 5)
        A = rng.uniform(0.15, 0.75) * M
        th = rng.uniform(0.05, 1.5)
        psi = rng.uniform(0, 2 * np.pi)
        phi = rng.normal(0, sc, (n_sub, n_cells)) if sc > 0 else 0.0
        mu = M + A * np.cos(t[:, None] + phi - psi) + rng.normal(0, 0.3, (n_sub, n_cells))
        mu = np.maximum(mu, 1e-6)
        x = rng.negative_binomial(1 / th, 1 / (1 + th * mu))  # mean mu, var mu + th mu^2
        m, v = x.mean(1), x.var(1, ddof=1)
        rows.append(estimate(m, v, t_hat) | dict(gene=g, M=M, A=A, theta=th))
    return pd.DataFrame(rows)


def estimate(m, v, t):
    """ORPHEUS moment solution from subject means m and variances v at phases t; sigma_c in hours (nan if none)."""
    X = np.column_stack([np.ones_like(t), np.cos(t), np.sin(t)])
    b = np.linalg.lstsq(X, m, rcond=None)[0]
    M, Am, acro = b[0], np.hypot(b[1], b[2]), np.arctan2(b[2], b[1])
    u = t - acro
    Z = np.column_stack([np.ones_like(u), np.cos(u), np.sin(u), np.cos(2 * u), np.sin(2 * u)])
    c = np.linalg.lstsq(Z, v, rcond=None)[0]
    B, C = c[1], c[3]
    theta = (B / Am - 1) / (2 * M)
    e = 1 + (2 * C / Am ** 2 - theta) / (1 + theta)
    s = -np.log(e) if 0 < e <= 1 else np.nan
    return dict(M_hat=M, Am=Am, theta_hat=theta, ratio=2 * C / Am ** 2, sigma_c_hat_h=np.sqrt(s) * H if s >= 0 else np.nan)


def analytic(sigma_c_h, sigma_s_h, theta):
    """Estimated sigma_c (h) from the moment solution when the true coefficients are attenuated by subject error."""
    s, s2 = (sigma_c_h / H) ** 2, (sigma_s_h / H) ** 2
    r = theta + (1 + theta) * (np.exp(-s) - 1)  # 2 C / A_m^2 without subject error
    e = 1 + (r * np.exp(-s2) - theta) / (1 + theta)
    return np.sqrt(-np.log(e)) * H if 0 < e <= 1 else np.nan


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--subjects", type=int, default=400)
    ap.add_argument("--cells", type=int, default=300)
    ap.add_argument("--genes", type=int, default=60)
    ap.add_argument("--reps", type=int, default=5)
    a = ap.parse_args()
    opt = cli.setup(a)
    if opt["quick"]:
        a.subjects, a.cells, a.genes, a.reps = a.subsample, 40, 10, 1
    out = out_dir(FOLDER, opt["quick"])
    rng = np.random.default_rng(a.seed)
    sim = []
    for sc in (1.0, 2.0, 3.0):
        for ss in (0.0, 1.0, 2.0, 3.0):
            for rep in range(a.reps):
                d = simulate(sc, ss, a.subjects, a.cells, a.genes, rng)
                sim.append(d.assign(sigma_c_h=sc, sigma_s_h=ss, rep=rep))
            print(sc, ss, flush=True)
    sim = pd.concat(sim, ignore_index=True)
    sim.to_csv(out / "orpheus_bias_simulation.csv", index=False)
    summ = (sim.groupby(["sigma_c_h", "sigma_s_h"])
            .agg(median_sigma_hat_h=("sigma_c_hat_h", "median"), solvable=("sigma_c_hat_h", lambda x: x.notna().mean()),
                 median_theta_hat=("theta_hat", "median"), median_theta=("theta", "median"))
            .reset_index())
    summ["analytic_median_theta"] = [analytic(r.sigma_c_h, r.sigma_s_h, r.median_theta) for r in summ.itertuples()]
    summ.to_csv(out / "orpheus_bias_summary.csv", index=False)
    grid = pd.DataFrame([dict(sigma_c_h=sc, sigma_s_h=ss, theta=th, sigma_hat_h=analytic(sc, ss, th))
                         for sc in np.arange(1, 4.01, 0.25) for ss in np.arange(0, 4.01, 0.25)
                         for th in (0.1, 0.3, 0.5, 1.0)])
    grid.to_csv(out / "orpheus_bias_analytic.csv", index=False)
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
