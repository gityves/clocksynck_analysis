"""V11: calibration of the von Mises moment model of per-cell clock amplitude and synchrony, by simulation at human depths.

Post hoc (2026-10-08, legacy V11; no amendment). Estimand: per-cell amplitude lam, synchrony R = I1(kappa) / I0(kappa), within-group
dispersion D = lam^2 (1 - R^2) and rho = lam R recovered from simulated single-cell counts (`clocksync.cells.vonmises.simulate`,
`group_stats`, `solve`). Gene baselines: mean clock-gene fractions of the tumour T/NK pseudobulks of the single-cell cancer cohorts
(`cancer_sc`); waveform: global GTEx (c_g = A_g e^{-i psi_g}); group phases uniform over 08-18 h. Grid: lam in {1.5, 3}, concentration
kappa in {0.5, 2, 8} (R = 0.24, 0.70, 0.93), clock UMI per cell in {2.4 (T/NK), 8 (fibroblast)}, metacells of m in {1, 10} random cells
(the moment rescaled by m), 4 replicates, 100 groups of 600 cells at depth 3000. Reported: bias and SD of lam, R, D, rho.
Note: the stored legacy summary (outputs/v11_vm_calibration/summary.csv) was made with an earlier estimator (per-group phase mu_s, `fit`); the
rotation-invariant estimator used here is the current one in both `clocksync` and legacy `sccc` (identical code, checked) and is unstable at
these depths (lam_hat up to 1e5, R_hat at 0 or 1 in single replicates): see the README.

    python compute_vonmises_calibration.py [--subsample N]     (N: groups per simulation; 1 replicate, 100 cells)
`--summarise-only`: recompute the summaries (also `vonmises_D_recovery.csv`: bias, SD and median ratio of D_hat, share of degenerate (lam, R) fits, by clock UMI and
metacell size) from the stored runs.
Outputs: paper/results/extended_data/output/calibration/vonmises_{runs,summary,D_recovery}.csv
Legacy: legacy/scripts/v11_vm_calibration.py (outputs/v11_vm_calibration/summary.csv)
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import special  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.cells import dispersion, vonmises  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "calibration"
N_GROUPS, N_CELLS, DEPTH, REPS = 100, 600, 3000.0, 4


def baseline(subsample) -> np.ndarray:
    """Mean clock-gene fractions of the tumour T/NK pseudobulks (all cohorts of `cancer_sc` with that compartment)."""
    pb, _ = data._builder("cancer_sc").from_legacy()
    t = pb[pb.compartment == "tumour:T_NK"]
    f = (t[list(cs.gtex().genes)].sum() / t.depth.sum()).to_numpy()
    return f / f.sum()


def split_projectors(c, n_splits, rng):
    """local: not in clocksync (`cells.vonmises._split_projectors` is private). Gene-half projectors of the rotation-invariant estimator."""
    X = np.column_stack([c.real, -c.imag])
    P = lambda Xh: dispersion._projector(Xh, np.ones(len(Xh)))
    return [(P(X[A]), P(X[B]), A, B) for A, B in dispersion.splits(X, n_splits, rng)]


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--summarise-only", action="store_true", help="recompute the summaries from the stored runs (no simulation)")
    args = ap.parse_args()
    opt = cli.setup(args)
    quick = opt["quick"]
    n_groups, n_cells, reps = (args.subsample, 100, 1) if quick else (N_GROUPS, N_CELLS, REPS)
    out = out_dir(FOLDER, quick)
    if args.summarise_only:
        return summarise(pd.read_csv(out / "vonmises_runs.csv"), out)
    wf = cs.gtex()
    c = wf.amp * np.exp(-1j * wf.psi)
    shape = baseline(args.subsample)
    rng = np.random.default_rng(0)
    LA = split_projectors(c, 30, rng)
    rows = []
    for lam, kappa, umi, m, rep in itertools.product((1.5, 3.0), (0.5, 2.0, 8.0), (2.4, 8.0), (1, 10), range(reps)):
        p = shape * umi / DEPTH
        mus = rng.uniform(8, 18, n_groups) * 2 * np.pi / 24
        sims = vonmises.simulate(n_groups, n_cells, c, p, lam, kappa, DEPTH, mus, rng)
        st = [vonmises.group_stats(C, N, LA, m=m, rng=rng) for C, N in sims]
        f = vonmises.solve(st, c, LA)
        R = float(special.ive(1, kappa) / special.ive(0, kappa))
        rows.append(dict(lam=lam, kappa=kappa, R=R, umi=umi, m=m, rep=rep, lam_hat=f["lam"], R_hat=f["R"], D=lam ** 2 * (1 - R ** 2),
                         D_hat=f["D"], ratio_hat=f["ratio"], rho=lam * R, rho_hat=f["rho"]))
        print(rows[-1], flush=True)
    summarise(pd.DataFrame(rows), out)


def summarise(d: pd.DataFrame, out) -> None:
    """Per-cell summary (legacy columns) and the D re-derivation table: what the RESULTS statements on V11 rest on."""
    d.to_csv(out / "vonmises_runs.csv", index=False)
    s = d.groupby(["lam", "R", "umi", "m"]).agg(lam_hat=("lam_hat", "mean"), lam_sd=("lam_hat", "std"), R_hat=("R_hat", "mean"),
                                                 R_sd=("R_hat", "std"), D=("D", "first"), D_hat=("D_hat", "mean"), D_sd=("D_hat", "std"),
                                                 rho=("rho", "first"), rho_hat=("rho_hat", "mean"), ratio_hat=("ratio_hat", "mean")).reset_index()
    s.to_csv(out / "vonmises_summary.csv", index=False)
    d = d.assign(err=d.D_hat - d.D, degenerate=(d.R_hat < 1e-3) | (d.R_hat > 0.999) | (d.lam_hat > 20))
    r = d.groupby(["umi", "m"]).agg(runs=("D", "size"), D_bias=("err", "mean"), D_sd_runs=("err", "std"), D_ratio_median=("D_hat", lambda x: float(np.median(x / d.loc[x.index, "D"]))),
                                    share_degenerate_lam_R=("degenerate", "mean")).reset_index()
    sd_within = d.groupby(["lam", "R", "umi", "m"]).D_hat.std().groupby(["umi", "m"]).mean().rename("D_sd_within_cell").reset_index()
    r = r.merge(sd_within, on=["umi", "m"])
    r.to_csv(out / "vonmises_D_recovery.csv", index=False)
    pd.set_option("display.width", 220)
    print(s.round(3).to_string(index=False))
    print(r.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
