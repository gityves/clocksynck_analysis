"""A132b (POST HOC, PI request): offset preservation, ventilated post-mortem cases minus brain-dead organ donors.

Registered in paper/docs/registered_experiments.md, section A132b. Fisher-Lee rho of the pair offsets vs fast deaths
(A109 H109.2 form) in each Hardy 0 group at equal n (the A132 organ-donor draw, default_rng(132 + seed)), over the
pairs with >= 20 shared donors in fast deaths and both groups (fixed mask). Statistic: rho(post-mortem) - rho(organ
donors); null: donor labels permuted between the two groups (1,000), two-sided p.
Controlled access (dbGaP phs000424.v8): phenotypes and time of death in memory only; aggregates in the output.

    python compute_posthoc_rho.py [--subsample N] [--seed S]
Output (paper/exploratory/A132_brain_death/output/): a132b_rho_difference.csv
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import compute_brain_death as B
from compute_brain_death import H, N, ST, cli, exploratory_dir, gtex


def rho(Zraw, M, fast, g, mask) -> float:
    return ST.fisher_lee(B.phase(Zraw, M, fast, mask), B.phase(Zraw, M, g, mask))


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_perm = 100 if run["quick"] else 1000
    out = exploratory_dir(B.FOLDER, run["quick"])
    rng = np.random.default_rng(13200 + args.seed)
    T = H.Tissues(1, args.seed, args.subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    Zraw, M, _ = N._raw_matrix(T)
    h0, fast = T.group_donors("vent"), T.group_donors("fast")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    ns = min(len(od), len(pm))
    rd = np.random.default_rng(132 + args.seed)  # first draw of A132's score(): the same organ donors
    od_ = rd.choice(od, ns, replace=False) if len(od) > ns else od
    pm_ = rd.choice(pm, ns, replace=False) if len(pm) > ns else pm
    mask = B.pair_mask(M, fast, od_, pm_)
    r_od, r_pm = rho(Zraw, M, fast, od_, mask), rho(Zraw, M, fast, pm_, mask)
    obs = r_pm - r_od
    pool, null = np.concatenate([pm_, od_]), []
    for _ in range(n_perm):
        p = rng.permutation(pool)
        null.append(rho(Zraw, M, fast, p[:ns], mask) - rho(Zraw, M, fast, p[ns:], mask))
    null = np.asarray(null)
    res = pd.DataFrame([dict(contrast="A132b_rho_postmortem_minus_organ_donor", rho_postmortem=r_pm,
                             rho_organ_donor=r_od, difference=obs, p_two_sided=N.pval(abs(obs), np.abs(null)),
                             null_sd=float(np.std(null, ddof=1)), n_each=ns, n_pairs=int(mask.sum()),
                             n_fast=len(fast), n_perm=len(null))])
    res.to_csv(out / "a132b_rho_difference.csv", index=False)
    print(res.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
