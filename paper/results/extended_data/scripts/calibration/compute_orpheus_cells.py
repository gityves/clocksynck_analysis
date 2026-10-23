"""A19 part 3b: does the within-sample dispersion of excitatory neurons rise in Alzheimer's disease? (the ORPHEUS discrepancy, data arm)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A19 part 3, written before running; exploratory there). ORPHEUS (Hollis et al. 2026) reads a lower
synchrony between cells of ROSMAP excitatory neurons in AD (phase dispersion 3.76 vs 3.23 h); A1 found no AD change of the between-donor
coupling (MSSM, RADC = ROSMAP, SEA-AD). 3b: within-sample in-plane dispersion D and synchrony R of excitatory neurons, AD minus control, on `AD_mssm`,
`AD_radc` and `AD_seaad_dlpfc` with the A18 per-cell estimator (`backbone/compute_per_cell.py`: cell groups of >= 100 cells,
D minus its acrophase-permuted null, rho^2 = s / (1 - R_W^2), donor bootstrap within groups; eligibility rule of A18). ORPHEUS predicts higher
D (lower R) in AD. This script runs those three targets (cells: `$SCCC_DATA/cells/{psychad_mssm,psychad_radc,seaad_dlpfc}`), keeps the excitatory
compartment and pools the relative change of D, d_D / D_control, by inverse variance across the three cohorts (SE from the bootstrap interval
width / 3.92; fixed effect; SEA-AD DLPFC is one cohort). Reading (plan): D unchanged and the 3a bias unable to produce the shift -> compatible
with no change in synchrony between cells; D higher in AD -> the two methods agree on desynchrony within samples. Also reported: the between-donor
side (A17 anchored coupling contrast of the excitatory neurons) is `resilience/compute_contrasts_a17.py`, not recomputed here.
Not ported: the host-uncertainty proxy of 3b (mean posterior covariance tr S_u of the host, AD against control; `Fit.host` has posterior means only).

    python compute_orpheus_cells.py [--subsample N] [--jobs 3] [--n-boot B]
`--subsample N`: N random donors per target, 5 bootstraps.
Outputs: paper/results/extended_data/output/calibration/orpheus_cells_{split,contrasts,pooled}.csv
Legacy: legacy/desync/common/cell_split.py (outputs: common_output/a18_cell_split_contrasts.csv, AD rows)
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("backbone")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

import compute_per_cell as PC  # noqa: E402

FOLDER = "calibration"
TARGETS = ("AD_mssm", "AD_radc", "AD_seaad_dlpfc")
EXCITATORY = "brain:excitatory"


def pool_relative(c: pd.DataFrame) -> pd.DataFrame:
    """Inverse-variance (fixed-effect) pool of the relative change d_D / D_control over the eligible cohorts."""
    c = c[c.eligible].copy()
    c["rel"] = c.d_D / c.D_control
    c["se"] = (c.d_D_hi - c.d_D_lo) / 3.92 / c.D_control.abs()
    w = 1 / c.se ** 2
    m = float((w * c.rel).sum() / w.sum())
    se = float(1 / np.sqrt(w.sum()))
    return pd.DataFrame([dict(n_cohorts=len(c), rel_change=m, lo=m - 1.96 * se, hi=m + 1.96 * se, se=se)])


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_boot = args.n_boot if args.n_boot is not None else (5 if run["quick"] else 100)
    out = out_dir(FOLDER, run["quick"])
    ns = SimpleNamespace(subsample=args.subsample, seed=args.seed)
    todo = [t for t in TARGETS if (DATA / "cells" / PC.DATASET[t] / "cells.parquet").exists()]
    res = Parallel(n_jobs=min(args.jobs, len(todo)))(delayed(PC.run_target)(t, n_boot, args.seed, ns) for t in todo)
    rows = pd.DataFrame([r for x, _ in res for r in x]).drop(columns="boot")
    con = pd.DataFrame([r for _, y in res for r in y])
    ex = con[con.compartment == EXCITATORY]
    rows.to_csv(out / "orpheus_cells_split.csv", index=False)
    con.to_csv(out / "orpheus_cells_contrasts.csv", index=False)
    pool = pool_relative(ex) if len(ex) else pd.DataFrame()
    pool.to_csv(out / "orpheus_cells_pooled.csv", index=False)
    pd.set_option("display.width", 220)
    print(ex[["target", "D_control", "D_case", "d_D", "d_D_lo", "d_D_hi", "R_control", "R_case", "eligible"]].round(3).to_string(index=False))
    print(pool.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
