"""A132b descriptive (POST HOC, PI request): pulled toward 0 / 12 h with the order kept, or reordered?

A shared non-temporal displacement added to running clocks pulls every pair phase toward its nearest pole (0 or 12 h)
without changing which pole or the order of the deviations from it; lost clock relations reorder them. Per group,
over the pairs with >= 20 shared donors in fast deaths and the group: the share of pairs on the same pole as in fast
deaths; Spearman of the signed deviations from the pole (fast vs group), joint tissue-label permutation p (1,000,
one-sided > 0); the median ratio of deviations (group / fast) over same-pole pairs with |deviation in fast| >= 1 h
(1 = no pull, 0 = full pull). Also the Rayleigh R of the times of death per group (uniform coverage of the cycle;
aggregate only). Groups: Hardy 0 ventilated post-mortem and organ donors (the A132 equal-n draw), all organ donors,
Hardy 4.
Controlled access (dbGaP phs000424.v8): phenotypes and time of death in memory only; aggregates in the output.

    python compute_posthoc_shrinkage.py [--subsample N] [--seed S]
Output (paper/exploratory/A132_brain_death/output/): a132b_shrinkage.csv
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import compute_brain_death as B
from compute_brain_death import H, N, cli, exploratory_dir, gtex


def deviation(theta: np.ndarray):
    """Nearest pole (0 or 12 h) and signed deviation from it (h) of pair phases (rad)."""
    h = theta * H.HOURS
    pole = np.abs(h) > 6
    return pole, np.where(pole, h - np.sign(h) * 12, h)


def readouts(th_f, th_g):
    pf, df = deviation(th_f)
    pg, dg = deviation(th_g)
    same = pf == pg
    big = same & (np.abs(df) >= 1)
    return dict(same_pole=float(same.mean()), spearman=float(spearmanr(df, dg)[0]),
                ratio_median=float(np.median(dg[big] / df[big])) if big.any() else np.nan, n_ratio=int(big.sum()))


def group_row(Zraw, M, fast, g, tod, label, rng, n_perm) -> dict:
    mask = B.pair_mask(M, fast, g)
    th_f = B.phase(Zraw, M, fast, mask)
    Hg = N._H_of(Zraw, M, g)
    r = readouts(th_f, H.pair_phase(Hg)[mask])
    null = [spearmanr(deviation(th_f)[1], deviation(H.permute_labels(Hg, rng.permutation(len(B.BODY)))[mask])[1])[0]
            for _ in range(n_perm)]
    w = 2 * np.pi * tod[g] / 24
    return dict(group=label, n=len(g), n_pairs=int(mask.sum()), **r, p_spearman=N.pval(r["spearman"], null),
                rayleigh_R_tod=float(np.abs(np.exp(1j * w).mean())), n_perm=n_perm)


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_perm = 100 if run["quick"] else 1000
    out = exploratory_dir(B.FOLDER, run["quick"])
    rng = np.random.default_rng(13210 + args.seed)
    T = H.Tissues(1, args.seed, args.subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    Zraw, M, _ = N._raw_matrix(T)
    h0, fast = T.group_donors("vent"), T.group_donors("fast")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    ns = min(len(od), len(pm))
    rd = np.random.default_rng(132 + args.seed)  # first draw of A132's score(): the same organ donors
    od_ = rd.choice(od, ns, replace=False) if len(od) > ns else od
    groups = (("fast", fast), ("ventilated_postmortem", pm), ("organ_donor", od_), ("organ_donor_full", od),
              ("hardy4", T.group_donors("H4")))
    res = pd.DataFrame([group_row(Zraw, M, fast, g, T.tod, lab, rng, n_perm) for lab, g in groups])
    res.to_csv(out / "a132b_shrinkage.csv", index=False)
    print(res.round(3).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
