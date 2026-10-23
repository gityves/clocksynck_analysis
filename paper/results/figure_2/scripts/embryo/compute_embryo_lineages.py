"""A108 (secondary) and A108c (V1): which mouse lineages already follow the time of day before birth?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A108 secondary; A108c V1, descriptive). Per trajectory (lineage)
pseudobulk of the timed window (clock time t = 24 x frac(E-day), one embryo per time point), the complex clock-plane
state of the lineage is regressed on a rotating wave:
  A108 secondary:  beta = mu + b E-day + a exp(i w t), E15.0-E18.75, trajectories with >= 50 cells in >= 12 embryos;
                   |a| tested by permuting t over embryos (E-day kept), 10,000; Benjamini-Hochberg over trajectories.
  A108c V1:        beta = mu_batch + a exp(i w t) in the early (E15.0-E16.75) and the late half (E17.0-E18.75), 8
                   embryos each; t permuted within batch, 10,000; "earlier" lineages: p < 0.05 in the early half.
Estimands: amplitude |a| (clock-plane units) and peak phase arg a (h); no correction in V1 (descriptive).

    python compute_embryo_lineages.py [--subsample N] [--seed S]
Outputs: paper/results/figure_2/output/embryo/a108_trajectories.csv, a108c_timed_halves.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a108_embryo.py (trajectory_rhythms), a108c_lineages.py
         (timed_half). local: the complex cosinor with permutation (_shared.timed_wave); the states are `Dataset.data`.
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
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "embryo"
HALVES = {"E15.0-E16.75": (15.0, 16.75), "E17.0-E18.75": (17.0, 18.75)}


def rhythms(pb, donors, don, min_cells, min_n, n_perm, seed, covariate=None, by_batch=False):
    """Per trajectory: |a|, phase and permutation p of the timed rotating wave."""
    t = S.eligible(pb, min_cells=min_cells, donors=don)
    n = t.groupby("compartment").donor.nunique()
    t = t[t.compartment.isin(n.index[n >= min_n])].reset_index(drop=True)
    ds = S.dataset(t, donors)
    d = ds.data
    rng = np.random.default_rng(seed)
    rows = []
    for k in sorted(set(t.compartment)):
        i = np.flatnonzero(d.compartment == k)
        if len(i) < min_n:
            continue
        r = d.resection[i]
        y = d.beta[i, 0] + 1j * d.beta[i, 1]
        tt = donors.t_clock.reindex(r).to_numpy()
        kw = dict(covariates=donors.day.reindex(r).to_numpy()) if covariate else {}
        if by_batch:
            bb = donors.batch.reindex(r).to_numpy()
            kw = dict(groups=bb, blocks=bb)
        rows.append(dict(trajectory=k, n=len(i), **S.timed_wave(y, tt, n_perm, rng, **kw)))
    return pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 500 if run["quick"] else 10_000
    pb, donors = data.load("mouse_qiu", level="trajectory", subsample=args.subsample, seed=args.seed)
    don = S.embryos(donors, *S.LATE)
    tr = rhythms(pb, donors, don, 50, 12 if not run["quick"] else 5, n_perm, 3, covariate=True)
    tr = tr.sort_values("p").reset_index(drop=True)
    tr["q_bh"] = cs.stats.bh(tr.p.to_numpy())
    tr.to_csv(out / "a108_trajectories.csv", index=False)
    print(tr.round(4).to_string(index=False), flush=True)
    halves = []
    for j, (h, (lo, hi)) in enumerate(HALVES.items()):
        try:
            halves.append(rhythms(pb, donors, S.embryos(donors, lo, hi), 50, 6, n_perm, 10 + j, by_batch=True)
                          .assign(half=h))
        except ValueError as e:  # too few embryos in a subsample
            print(h, "skipped:", e)
    if halves:
        v1 = pd.concat(halves).pivot(index="trajectory", columns="half", values=["amp", "phase_h", "p"])
        v1.columns = [f"{a}_{b}" for a, b in v1.columns]
        v1 = v1.sort_values(v1.columns[v1.columns.str.startswith("p_")][0])
        v1.to_csv(out / "a108c_timed_halves.csv")
        print(v1.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
