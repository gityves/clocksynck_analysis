"""PI Q3 follow-up (2026-10-27): the A21 R21.4 epithelial phase shift with lymphatic EC removed from the anchors.

Re-estimation, not a new hypothesis. A22 V17.2 (positive control) failed: lymphatic EC, known to decouple after DSS, was not
dropped by the confirmation rule (coherence 0.448 -> 0.048), so the R-DC host includes a decoupled anchor. Here R21.4 is
recomputed exactly as in `figure_1/scripts/species/compute_colon_injury.py` (untimed: each arm's own anchored host; timed: arg of
the known-time loadings; both relative to the circular mean of the anchors), with the anchors fibroblast, smooth_muscle, pericyte
and endothelial only. Lymphatic EC becomes a target, so its own shift is reported too. 500 animal bootstraps (as M3).

    python a21_shift_without_lymphatic.py [--subsample N] [--jobs J]
Output: paper/exploratory/WP0_integrity/output/a21_shift_without_lymphatic.csv
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
SPECIES = PAPER / "results/figure_1/scripts/species"
sys.path.insert(0, str(SPECIES))

from common import cli, data  # noqa: E402

_spec = importlib.util.spec_from_file_location("colon_injury", SPECIES / "compute_colon_injury.py")
A21 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A21)

OUT = PAPER / "exploratory/WP0_integrity/output"


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    A21.STROMAL_VASCULAR = tuple(k for k in A21.STROMAL_VASCULAR if k != "lymphatic_EC")  # read by arm_dataset
    pb, donors = data.load("colon_chrono")
    n_boot = 20 if run["quick"] else 500
    ph = A21.r214(pb, donors, donors.zt, n_boot, args.jobs)
    ph.insert(0, "anchors", ",".join(A21.STROMAL_VASCULAR))
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / ("a21_shift_without_lymphatic_quick.csv" if run["quick"] else "a21_shift_without_lymphatic.csv")
    ph.to_csv(f, index=False)
    print(ph.round(2).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
