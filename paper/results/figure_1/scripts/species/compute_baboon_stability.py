"""A36: is the baboon host identified? Stability bootstrap, descriptive gate and leave-one-animal-out similarity.

Pre-registered as A36 (2026-10-15; docs/PLAN_internal_desync.md (archive)). The primary (ordering), R2 and R3 are in `compute_baboon.py`; this script ports the
remaining legacy readouts of `score.py`: (i) host stability, the share of 200 donor-bootstrap refits of the host (63 tissues as anchors) whose state
correlates > 0.7 with the full-data host (`Dataset.stability`), with the median |complex correlation|; (ii) the descriptive gate of the 63 tissues (200
permutations; a failed gate would be uninformative); (iii) similarity of the host fitted without each animal to the full host (|complex correlation|, over the 11
common animals), next to the ordering readouts without that animal (`a36_stability.csv`). Estimands: stability, median similarity, clock share and null.

    python compute_baboon_stability.py [--subsample N] [--n-boot B] [--n-perm P]
Outputs: paper/results/figure_1/output/species/a36_host_stability.csv (stability, gate and LOAO similarity rows)
Legacy: legacy/desync/validation/baboon/scripts/score.py -> output/score.log (stability, gate), stability.csv (sim_to_full)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "species"


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = args.n_boot or (10 if run["quick"] else 200)
    n_perm = args.n_perm or (50 if run["quick"] else 200)
    pb, donors = data.load("baboon")
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, donors=donors, name="baboon")
    rows = []
    st = ds.stability(n_boot=n_boot, seed=args.seed)
    rows.append(dict(readout="host stability (bootstrap)", **st))
    g = ds.gate(n_perm=n_perm, seed=args.seed)
    rows.append(dict(readout="gate (63 tissues)", clock_share=g.clock_share, share_null_mean=g.share_null_mean, p_pairing=g.p_pairing,
                     p_acrophase=g.p_acrophase, n_pairs=g.n_pairs, n=g.n_donors, passes=g.passes))
    full = ds.host()[["x", "y"]]
    sims = []
    for a in (donors.index[:3] if run["quick"] else donors.index):
        h = ds.subset(donors=[x for x in donors.index if x != a]).host()[["x", "y"]]
        common = h.index.intersection(full.index)
        sims.append(abs(cs.stats.ccorr(h.loc[common].to_numpy(), full.loc[common].to_numpy())))
    rows.append(dict(readout="leave-one-animal-out host similarity to the full host", n=len(sims), min=float(np.min(sims)),
                     median=float(np.median(sims))))
    R = pd.DataFrame(rows)
    R.to_csv(out / "a36_host_stability.csv", index=False)
    print(R.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
