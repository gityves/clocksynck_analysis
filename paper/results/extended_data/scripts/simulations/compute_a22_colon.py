"""A22 V17.2: the group-confirmation rule on the mouse colon (positive control, written before the human contrasts).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A22, V17 part 2): fixed stromal-vascular anchors (fibroblast, smooth_muscle,
pericyte, endothelial, lymphatic_EC present), C-DC (control) against R-DC (case), 8 mice per arm, classes present in both arms. Rule
(`host_identity/compute_group_anchors.confirm`, the coherence rule of A22: pair coherence among the anchors only, 200 shuffled
pairings, 200 donor bootstraps per group; delta_k = e_case - e_control, SE = hypot; drop if delta < -1.96 SE and e_case <= 0, or the
mirror image): it must drop lymphatic EC. The A21 contrasts are then recomputed on the confirmed host (reported: R21.2 Spearman and the
attenuation of the untimed against the timed d*; descriptive, the lymphatic result is known).
Run note of A22 (2026-10-13): V17.2 fails (lymphatic EC not dropped); this script reports it as run. When an anchor is dropped the
recomputation of the A21 contrasts is printed as not made (`species/compute_colon_injury.py` has no anchor-drop option).

    python compute_a22_colon.py [--subsample N] [--seed S]   (seed S + 3000, as the legacy)
`--subsample N`: 20 shuffles and N bootstraps (check).
Outputs: paper/results/extended_data/output/simulations/a22_colon_confirm.csv
Legacy: legacy/desync/common/a22_group_anchors.py::colon (common_output/a22_colon_confirm.csv)
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("host_identity")))
sys.path.insert(0, str(group_dir("atlas_fits")))
from compute_group_anchors import confirm  # noqa: E402

FOLDER = "simulations"
STROMAL_VASCULAR = ("fibroblast", "smooth_muscle", "pericyte", "endothelial", "lymphatic_EC")


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    warnings.filterwarnings("ignore")
    pb, donors = data.load("colon_chrono")
    p = pb[pb.donor.isin(donors.index[donors.arm.isin(["C-DC", "R-DC"])])]
    both = p.merge(donors[["arm"]], left_on="donor", right_index=True).groupby("compartment").arm.nunique() == 2
    p = p[p.compartment.isin(both.index[both])].reset_index(drop=True)
    anchors = [k for k in sorted(p.compartment.unique()) if k in STROMAL_VASCULAR]
    lab = donors.arm.loc[p.donor.unique()].map({"C-DC": "control", "R-DC": "case"})
    n_perm, n_boot = (20, args.subsample) if run["quick"] else (200, 200)
    r = confirm(SimpleNamespace(pb=p), anchors, lab, 3000 + args.seed, n_perm, n_boot)
    r.to_csv(out / "a22_colon_confirm.csv", index=False)
    print(r.round(3).to_string(index=False))
    drop = list(r.anchor[r["drop"]])
    print("V17.2 lymphatic_EC dropped:", "lymphatic_EC" in drop, "| anchors dropped:", drop)
    if drop:
        print("recomputation of the A21 contrasts on the confirmed host: not made (see docstring)")


if __name__ == "__main__":
    main()
