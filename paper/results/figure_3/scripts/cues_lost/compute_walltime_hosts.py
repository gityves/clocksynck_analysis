"""A27 part 2 (D4/D5 sanity), ventilated donors: does the internal clock state read the recorded time of death?

Post hoc (A27, disclosed there; the readout is reused as context in A46, A92, A93). The host is fitted without any time on
the clock genes of the 15 reference tissues (`Dataset`, every tissue an anchor), separately for fast deaths (Hardy 1-2),
ventilator cases (Hardy 0) and, for context, the intermediate and slow deaths (Hardy 3, Hardy 4). `Dataset.timed(tod)`:
R = |mean exp(i (host phase - w t_death))| (1: the host phase equals clock time up to one offset), its permutation p (times
shuffled over donors, 10,000), the common offset, the median absolute error after it and the circular correlations.
In ventilated donors the host is expected to be decoupled from wall time (legacy R 0.14) while it stays shared between the
tissues (the gate passes).

    python compute_walltime_hosts.py [--subsample N]
Outputs: paper/results/figure_2/output/organs/walltime_hosts.csv (one row per cohort; aggregates only)
Legacy: legacy/desync/validation/glucocorticoid/scripts/diag5_int_ext.py (sanity lines of `main`), a46_fast_context.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "cues_lost"
GROUPS = {"fast": (1, 2), "ventilator": (0,), "hardy3": (3,), "hardy4": (4,)}


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 500 if run["quick"] else 10_000
    rows = []
    for name, hardy in GROUPS.items():
        pb, donors = gtex.load(tissues=list(S.V2_TISSUES), hardy=hardy, genes="clock", depth="library_size")
        if args.subsample:
            keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:max(args.subsample, 20)]
            pb = pb[pb.donor.isin(keep)]
        ds = cs.Dataset(pb, anchors=sorted(pb.compartment.unique()), name=name)
        g = ds.gate(n_perm=100 if run["quick"] else 200)
        r = ds.timed(donors.tod, n_perm=n_perm)
        rows.append(dict(cohort=name, donors=pb.donor.nunique(), p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                         clock_share=g.clock_share, **r.to_dict()))
        print(pd.DataFrame(rows[-1:]).round(3).to_string(index=False), flush=True)
    pd.DataFrame(rows).to_csv(out / "walltime_hosts.csv", index=False)


if __name__ == "__main__":
    main()
