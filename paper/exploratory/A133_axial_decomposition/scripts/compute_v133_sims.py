"""A133 V133.1: simulation gate of the circular / axial decomposition (no GTEx donor-level data is used).

Generative model (A109 `simulate_group` extended with a host phase law and per-tissue directions): for donor d, tissue a
    z_da = a_a [ sqrt(c f) s_d e^{i g_a} + sqrt(c (1 - f)) e^{i (theta_d + phi_a)} + sqrt(1 - c) xi_da ] + noise_da,
s_d ~ N(0, 1) real displacement, theta_d host phase (uniform, or von Mises with Rayleigh R), xi circular Gaussian,
noise Poisson-like (circular Gaussian with variance V_a), states centred within tissue over the observed donors.
  a_a^2, V_a: ventilated (Hardy 0) corrected clock-plane variance and Poisson variance of A109
  (`clock_plane_variance.csv`, group vent, full); c = 0.45 (A109 circularity_sims C_target 0.451).
  phi_a: tissue phases fitted (circular least squares, leading eigenvector) to the A109 fast-death pair offsets
  (`pair_relations.csv`, aggregates; spread SD of pair offsets 1.3 h).
  g_a: 'common' (all 17.7 h, the glucocorticoid direction, as A109) or 'pertissue' (g_a ~ U[0, pi) per tissue and rep).
Estimator: `axial.circular_part` on C, P from A109 `cross_moments` (pairs with >= 20 shared donors observed in P).
Readouts: Fisher-Lee rho of the 66 pair offsets (raw C, C_circ) with truth, and of C_circ with raw C; median folded
distance to 0 / 12 h (h).

    python compute_v133_sims.py [--subsample] [--seed S]
Outputs: output/v133_grid.csv, output/v133_gate.csv (output/subsample/ for --subsample).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "A109_coupling_nulls" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import brentq  # noqa: E402
from scipy.special import i0e, i1e  # noqa: E402

import _helpers as H  # noqa: E402
from axial import circular_part  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

FOLDER = "A133_axial_decomposition"
A109 = HERE.parents[1] / "A109_coupling_nulls" / "output"
BODY, PAIRS, MIN_N = H.BODY, H.PAIRS, H.MIN_N
C_COUP, GC_RAD = 0.45, 2 * np.pi * H.GC_H / 24
FS, NS, RS, US = (0, 0.2, 0.44, 0.6, 0.8), (79, 298), (0.0, 0.15, 0.3), ("common", "pertissue")


def design():
    """(a, V, phi): amplitudes, Poisson variances (ventilated, A109) and tissue phases (rad) of the fast-death pairs."""
    v = pd.read_csv(A109 / "clock_plane_variance.csv")
    v = v[(v.group == "vent") & (v.donors == "full")].set_index("tissue").loc[BODY]
    rel = pd.read_csv(A109 / "pair_relations.csv")
    rel = rel[rel.group == "fast"]
    idx = {t: i for i, t in enumerate(BODY)}
    M = np.eye(len(BODY), dtype=complex)
    for _, r in rel.iterrows():
        a, b = idx[r.tissue_a], idx[r.tissue_b]
        M[a, b] = np.exp(1j * r.offset_h / H.HOURS)
        M[b, a] = np.conj(M[a, b])
    w, U = np.linalg.eigh(M)
    phi = np.angle(U[:, -1] * np.conj(U[0, -1]))
    return np.sqrt(v.clock_plane_variance.to_numpy()), v.V.to_numpy(), phi


def kappa(R: float) -> float:
    return 0.0 if R == 0 else brentq(lambda k: i1e(k) / i0e(k) - R, 1e-6, 50)


def simulate(n, f, R, u, a, V, phi, rng, presence=1.0):
    T = len(a)
    theta = rng.vonmises(0, kappa(R), n)[:, None] if R > 0 else rng.uniform(0, 2 * np.pi, n)[:, None]
    g = np.full(T, GC_RAD) if u == "common" else rng.uniform(0, np.pi, T)
    s = rng.standard_normal(n)[:, None]
    xi = (rng.standard_normal((n, T)) + 1j * rng.standard_normal((n, T))) / np.sqrt(2)
    z = a * (np.sqrt(C_COUP * f) * s * np.exp(1j * g) + np.sqrt(C_COUP * (1 - f)) * np.exp(1j * (theta + phi))
             + np.sqrt(1 - C_COUP) * xi)
    z = z + np.sqrt(V / 2) * (rng.standard_normal((n, T)) + 1j * rng.standard_normal((n, T)))
    M = (rng.uniform(size=(n, T)) < presence).astype(float)
    z = z * M
    mu = z.sum(0) / np.maximum(M.sum(0), 1)
    return (z - mu * M) * M, M


def pair_phases(Cm):
    return np.array([np.angle(Cm[a, b]) for a, b in PAIRS])


def folded(theta):
    return H.folded_distance(theta)


def one_rep(n, f, R, u, a, V, phi, rng, presence=1.0):
    Z, M = simulate(n, f, R, u, a, V, phi, rng, presence)
    C, P, N = H.cross_moments(Z, M)
    Cc, _ = circular_part(C, P, N >= MIN_N)
    raw, circ = pair_phases(C), pair_phases(Cc)
    truth = np.array([phi[x] - phi[y] for x, y in PAIRS])
    fl = H.fisher_lee_fast
    return dict(rho_raw=fl(raw, truth), rho_circ=fl(circ, truth), rho_circ_raw=fl(circ, raw),
                fold_raw=folded(raw), fold_circ=folded(circ), fold_truth=folded(truth))


def grid(n_rep, seed, cells=None, presence=1.0):
    a, V, phi = design()
    rows = []
    for k, (f, n, R, u) in enumerate(cells or itertools.product(FS, NS, RS, US)):
        rng = np.random.default_rng(np.random.SeedSequence([seed, 133, int(presence * 100), k]))
        r = pd.DataFrame([one_rep(n, f, R, u, a, V, phi, rng, presence) for _ in range(n_rep)])
        rows.append(dict(f=f, n=n, R=R, u=u, presence=presence, reps=n_rep, **r.mean().to_dict(),
                         **{f"{c}_sd": r[c].std() for c in r.columns if c != "fold_truth"},
                         rho_circ_gt_raw=float((r.rho_circ > r.rho_raw).mean())))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.set_defaults(subsample=None)
    args = ap.parse_args(argv)
    cli.setup(args)
    quick = args.subsample is not None
    out = exploratory_dir(FOLDER, quick)
    n_rep = 10 if quick else 100
    g = grid(n_rep, args.seed)
    # extra registered cells: f = 0.5, n = 79, R = 0 (criterion i), 12 tissues
    g5 = grid(n_rep, args.seed, [(0.5, 79, 0.0, u) for u in US])
    sup = grid(n_rep, args.seed, [(f, n, 0.0, u) for f in (0, 0.44) for n in NS for u in US], presence=0.8)
    full = pd.concat([g, g5, sup], ignore_index=True)
    full.to_csv(out / "v133_grid.csv", index=False)
    rows = []

    def add(name, cells, col, thr, extra=""):
        val = cells[col].min()
        rows.append(dict(criterion=name, value=round(float(val), 4), threshold=thr, passed=bool(val >= thr),
                         detail=extra + " cells: " + "; ".join(
                             f"n{int(r.n)} {r.u} R{r.R} f{r.f}: {getattr(r, col):.3f}" for r in cells.itertuples())))

    p1 = full[full.presence == 1]
    c5 = p1[(p1.f == 0.5) & (p1.n == 79)]
    add("i. truth recovery, f=0.5, n=79, rho(C_circ, truth) (worst u)", c5, "rho_circ", 0.8)
    ok = bool((c5.rho_circ > c5.rho_raw).all())
    rows.append(dict(criterion="i. raw offsets do worse (rho_raw < rho_circ in every u)", value=round(float((c5.rho_circ - c5.rho_raw).min()), 4),
                     threshold="> 0", passed=ok,
                     detail="; ".join(f"{r.u}: raw {r.rho_raw:.3f} circ {r.rho_circ:.3f}" for r in c5.itertuples())))
    c0 = p1[(p1.f == 0) & (p1.R == 0) & (p1.n == 79)]
    add("ii. no displacement, f=0, n=79, rho(C_circ, raw) (worst u)", c0, "rho_circ_raw", 0.95)
    c3 = p1[(p1.f == 0) & (p1.R > 0) & (p1.n == 79)]
    add("iii. concentrated host, f=0, n=79, R 0.15-0.3, rho(C_circ, truth) (worst cell)", c3, "rho_circ", 0.9)
    rows.append(dict(criterion="iii. companion: raw rho with truth, same cells (min)", value=round(float(c3.rho_raw.min()), 4),
                     threshold="descriptive", passed=np.nan, detail=""))
    pd.DataFrame(rows).to_csv(out / "v133_gate.csv", index=False)
    print(pd.DataFrame(rows)[["criterion", "value", "threshold", "passed"]].to_string(index=False))
    cols = ["f", "n", "R", "u", "presence", "rho_raw", "rho_circ", "rho_circ_raw", "fold_raw", "fold_circ", "fold_truth"]
    print(full[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
