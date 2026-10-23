"""A86: is a clock observable within malignant cells? Single-cell cancer cohorts with copy-number calls.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A86). Per compartment, pooled over the cohorts (donor = study|patient,
stratum = study|stratum): log((n + 0.5) / depth) of the 14 clock genes centred within study x stratum; clock-plane share
of the covariance with the Poisson diagonal (mean 1 / (n + 0.5)) subtracted; null from 1,000 permuted-acrophase
planes; compartments with >= 100 units (malignant per cancer type with >= 30, descriptive).
H86: malignant minus adjacent-normal epithelium share in the patients with both; two-sided p from within-patient label
swaps (1,000).
Reading (as A83): a lower share is not by itself a broken clock; a share at the null is "no observable clock".

    python compute_malignant_gate.py [--subsample N] [--seed S]
Outputs (paper/results/figure_5/output/cancer/): a86_compartments.csv, a86_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a86_malignant_gate.py.
The pooled, Poisson-corrected gate has no public clocksync call: local helpers (_shared.centred_units / pooled_share).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
SEED = 86


def gate(d: pd.DataFrame, rng, n_perm: int, P0: np.ndarray) -> dict:
    """local: pooled clock-plane share (Poisson diagonal subtracted) against permuted-acrophase planes."""
    y, pv, k = sh.centred_units(d, "stratum")
    obs = sh.pooled_share(y, pv, k, P0)
    null = np.array([sh.pooled_share(y, pv, k, sh.plane(sh.WF.permuted(rng))) for _ in range(n_perm)])
    return dict(units=len(y), share=obs, share_null=float(null.mean()), p=float((1 + np.sum(null >= obs)) / (1 + n_perm)))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 50 if run["quick"] else 1000
    min_units, min_cancer = (10, 8) if run["quick"] else (100, 30)
    rng = np.random.default_rng(SEED + args.seed)
    P0 = sh.plane()
    pb, don = sh.load_sc(args.subsample, args.seed)
    rows = []
    for comp, x in pb.groupby("compartment"):
        if len(x) >= min_units:
            rows.append(dict(compartment=comp, studies=x.study.nunique(), **gate(x, rng, n_perm, P0)))
            print(rows[-1], flush=True)
    mal = pb[pb.compartment == sh.MAL].assign(cancer=lambda x: x.donor.map(don.cancer).str.replace(" cancer", ""))
    for c, x in mal.groupby("cancer"):
        if len(x) >= min_cancer:
            rows.append(dict(compartment=f"{sh.MAL}|{c}", studies=x.study.nunique(), **gate(x, rng, n_perm, P0)))
    g = pd.DataFrame(rows)

    # H86: the same patients, within-patient label swaps
    m, nrm = pb[pb.compartment == sh.MAL], pb[pb.compartment == sh.NORM]
    both = sorted(set(m.donor) & set(nrm.donor))
    m = m.set_index("donor").loc[both].groupby(level=0).head(1).reset_index()
    nrm = nrm.set_index("donor").loc[both].groupby(level=0).head(1).reset_index()

    def diff(flip):
        a, b = pd.concat([nrm[flip], m[~flip]]), pd.concat([m[flip], nrm[~flip]])
        return sh.pooled_share(*sh.centred_units(a, "stratum"), P0) - sh.pooled_share(*sh.centred_units(b, "stratum"), P0)

    obs = diff(np.zeros(len(m), bool))
    null = np.array([diff(rng.random(len(m)) < 0.5) for _ in range(n_perm)])
    two = float(min(1.0, 2 * min(cs.stats.perm_p(null, obs, "less", tie=0.0), cs.stats.perm_p(-null, -obs, "less", tie=0.0))))
    v = dict(patients=len(m), studies=m.study.nunique(), malignant_minus_normal=obs, p_two=two, n_swap=n_perm)
    g.to_csv(out / "a86_compartments.csv", index=False)
    pd.Series(v).to_csv(out / "a86_verdict.csv", header=False)
    print(g.round(3).to_string(), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
