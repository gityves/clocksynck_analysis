"""A108 (descriptive) and A108c (V2): when and in which lineage does the mouse embryo carry a shared clock state?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A108 emergence; A108c V2, descriptive).
Estimand: the model-free gate (shared clock-plane covariance of the cell types of an embryo against shuffled
pairing and permuted acrophases; pass = p_pairing and p_acrophase < 0.05; 1000 permutations).
  A108:   all cell types (>= 50 cells, present in >= 75 % of the window's embryos) in four developmental windows
          E8.0-E10.0, E10.25-E12.75, E13.0-E14.75, E15.0-E18.75 (batch is the stratum). Earlier windows were staged
          by somites or limb morphology at unrecorded clock times, so a failed gate is uninformative.
  A108c:  per lineage (trajectory) with >= 3 cell types meeting the window rule, the gate over its cell types in each
          window; "first appearance" = earliest window that passes. No correction (descriptive).

    python compute_embryo_emergence.py [--subsample N] [--jobs 4] [--seed S]
Outputs: paper/results/figure_2/output/embryo/a108_gates.csv, a108c_lineage_gates.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a108_embryo.py (gates), a108c_lineages.py (lineage_gate).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "embryo"


def window_gate(pb, donors, name, lo, hi, seed, n_perm, traj=None, lineage=None):
    limit_threads(1)
    t, don = S.window(pb, donors, lo, hi)
    key = dict(window=name) if lineage is None else dict(lineage=lineage, window=name)
    if lineage is not None:
        t = t[t.compartment.map(traj) == lineage].reset_index(drop=True)
    n = t.compartment.nunique()
    if n < 3 or t.donor.nunique() < 5:
        return dict(**key, embryos=len(don), cell_types=n, status="too few")
    g = S.dataset(t).gate(n_perm=n_perm, seed=seed)
    return dict(**key, embryos=t.donor.nunique(), cell_types=n, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                clock_share=g.clock_share, share_null_mean=g.share_null_mean, passed=g.passes)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    pb, donors = data.load("mouse_qiu", subsample=args.subsample, seed=args.seed)
    rows = [window_gate(pb, donors, w, lo, hi, j, run["n_perm"]) for j, (w, (lo, hi)) in enumerate(S.WINDOWS.items())]
    g = pd.DataFrame(rows)
    g.to_csv(out / "a108_gates.csv", index=False)
    print(g.round(4).to_string(index=False), flush=True)
    traj = S.trajectory_of()
    jobs = [(lin, w, lo, hi) for lin in sorted(traj.unique()) for w, (lo, hi) in S.WINDOWS.items()]
    res = Parallel(n_jobs=args.jobs)(delayed(window_gate)(pb, donors, w, lo, hi, 100 + j, run["n_perm"], traj, lin)
                                     for j, (lin, w, lo, hi) in enumerate(jobs))
    v2 = pd.DataFrame(res)
    v2.to_csv(out / "a108c_lineage_gates.csv", index=False)
    print(v2.dropna(subset=["p_acrophase"]).round(3).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
