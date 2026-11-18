"""A314 (registered in `paper/docs/registered_experiments.md`; descriptive, post hoc split): does the body's clock read the
hour of death in ventilated donors with a living brain and in brain-dead organ donors?

GTEx Hardy 0 split by COHORT (Postmortem; Organ Donor (OPO)). Per group, as A27 `compute_walltime_hosts.py`: the 15
reference tissues, clock genes, library-size depth; host fitted without time, every tissue an anchor; `Dataset.timed`
(R = |mean exp(i (host phase - w t_death))|, permutation p, 10,000). Added: the null median of R (1,000 of the same kind
of permutations), and organ donors at the post-mortem n (20 random draws, seed 314). Aggregates only (time of death
and COHORT are controlled access, held in memory).

    python score_a314.py
Outputs (`../output/`): groups.csv (one row per group), draws.csv (organ donors at equal n, one row per draw).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER.parent))
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.timed import W, agreement  # noqa: E402
from common import cli, gtex  # noqa: E402
from results.figure_2.scripts.organs._shared import V2_TISSUES  # noqa: E402

OUT = HERE.parents[1] / "output"
GROUPS = {"ventilated_postmortem": "Postmortem", "organ_donor": "Organ Donor (OPO)"}
N_PERM, N_NULL, N_DRAW, SEED = 10_000, 1_000, 20, 314


def readout(pb, tod):
    ds = cs.Dataset(pb, anchors=sorted(pb.compartment.unique()))
    h = A.host(ds.data, ds.anchors)
    phi = pd.Series(np.angle(h.f[:, 0] + 1j * h.f[:, 1]), index=h.resections)
    t = tod.reindex(phi.index).astype(float)
    ok = t.notna().to_numpy()
    a, b = phi.to_numpy()[ok], t.to_numpy()[ok] * W
    r = agreement(a, b, n_perm=N_PERM, seed=SEED)
    rng = np.random.default_rng(SEED)
    null = [abs(np.exp(1j * (a - rng.permutation(b))).mean()) for _ in range(N_NULL)]
    return dict(donors=int(ok.sum()), R=r.R, p_R=r.p_R, null_median_R=float(np.median(null)),
                null_q95_R=float(np.quantile(null, 0.95)), offset_h=r.offset_h)


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    cli.setup(args)
    OUT.mkdir(parents=True, exist_ok=True)
    pb, donors = gtex.load(tissues=list(V2_TISSUES), hardy=(0,), genes="clock", depth="library_size")
    cohort = gtex.phenotypes(["COHORT"]).COHORT
    pb = pb.assign(cohort=pb.donor.astype(str).map(cohort))
    tod = donors.tod.rename(index=str)
    rows, groups = [], {}
    for name, c in GROUPS.items():
        g = pb[pb.cohort == c].drop(columns="cohort")
        groups[name] = g
        rows.append(dict(group=name, **readout(g, tod)))
        print(pd.DataFrame(rows[-1:]).round(3).to_string(index=False), flush=True)
    n_pm = groups["ventilated_postmortem"].donor.nunique()
    od = groups["organ_donor"]
    ids = np.array(sorted(od.donor.unique()))
    rng = np.random.default_rng(SEED)
    draws = []
    for k in range(N_DRAW):
        keep = rng.choice(ids, n_pm, replace=False)
        draws.append(dict(draw=k, **readout(od[od.donor.isin(keep)], tod)))
    d = pd.DataFrame(draws)
    d.to_csv(OUT / "draws.csv", index=False)
    rows.append(dict(group="organ_donor_equal_n", donors=n_pm, R=d.R.median(), p_R=d.p_R.median(),
                     null_median_R=d.null_median_R.median(), null_q95_R=d.null_q95_R.median(),
                     offset_h=np.nan, R_min=d.R.min(), R_max=d.R.max()))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "groups.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
