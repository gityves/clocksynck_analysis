"""A105 (H105.2 reference, descriptive): how often does a gate over 14 donors detect a sharing that is present?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A105: GTEx reference). The four tissues of the Delorey organs (Lung,
Heart - Left Ventricle, Liver, Kidney - Cortex) in GTEx, fast deaths (Hardy 1-2) and ventilator cases (Hardy 0) apart:
among donors with >= 3 of the four tissues, 50 random subsets of 14 donors, each gated (200 permutations);
the readout is the share of subsets that pass (p_pairing and p_acrophase < 0.05). It calibrates the power of the
n = 14 gate of H105.2. Controlled access: aggregates only (pass rates and per-subset p values, no donor identity).

    python compute_gtex_reference.py [--subsample N] [--jobs 4] [--seed S]     (--subsample: number of subsets)
Outputs: paper/results/figure_4/output/natural_experiments/a105_gtex_reference.csv (one row per subset),
         a105_gtex_reference_summary.csv (pass rate per cohort)
Legacy:  legacy/desync/validation/natural_experiments/scripts/a105_covid.py (gtex_reference).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "natural_experiments"
TISSUES = ["Lung", "Heart - Left Ventricle", "Liver", "Kidney - Cortex"]
COHORTS = {"fast": gtex.FAST, "ventilator": gtex.VENTILATOR}


def one(cohort, j, t, subset, n_perm):
    limit_threads(1)
    u = t[t.donor.isin(subset)].reset_index(drop=True)
    g = S.dataset(u).gate(n_perm=n_perm, seed=j)
    return dict(cohort=cohort, j=j, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, passed=g.passes)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_sub = args.subsample or 50
    rows = []
    for cohort, hardy in COHORTS.items():
        t, _ = gtex.load(tissues=TISSUES, hardy=hardy)
        n = t.groupby("donor").compartment.nunique()
        elig = n.index[n >= 3].to_numpy()
        rng = np.random.default_rng(11)
        subs = [rng.choice(elig, 14, replace=False) for _ in range(n_sub)]
        rows += Parallel(n_jobs=args.jobs)(delayed(one)(cohort, j, t, s, 200) for j, s in enumerate(subs))
        print(cohort, "eligible donors", len(elig), flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(out / "a105_gtex_reference.csv", index=False)
    summ = res.groupby("cohort").passed.agg(subsets="size", pass_rate="mean")
    summ.to_csv(out / "a105_gtex_reference_summary.csv")
    print(summ, flush=True)


if __name__ == "__main__":
    main()
