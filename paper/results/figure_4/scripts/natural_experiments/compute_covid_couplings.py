"""A105 (H105.2 coupling, D1): couplings of the organs and of the lung classes in fatal COVID-19.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A105: H105.2 per-organ coupling, descriptive; D1 per-class coupling
in COVID-19 vs surgical controls, descriptive).
  H105.2  Delorey whole-organ nuclei (>= 50 nuclei): coupling of each organ to the host of the other organs
          (`Dataset.fit`, every organ an anchor, leave-one-out), c = z2 V / (z2 V + 2 psi), donor-bootstrap interval
          (host refitted) and a coupling p (organ's donors permuted against a fixed host).
  D1      Melms lung classes (>= 20 cells), COVID-19 vs the 7 living surgical controls (procurement-confounded):
          `Dataset.contrast` on the endothelial + fibroblast host fitted on the pooled donors, so both arms are read on
          one host: c per arm, d, d* at a common spread, and the spread-free log ratios of z2 and psi.
          (legacy `per_arm` fitted each arm on its own host: `Dataset.fit` needs >= 3 anchors, so the pooled host of
          `contrast` is used; see the README.)
No gate is run here (gates are in compute_covid_gates.py); coupling is reported whatever its p, never as a verdict.

    python compute_covid_couplings.py [--subsample N] [--seed S]
Outputs: paper/results/figure_4/output/natural_experiments/a105_delorey_organ_coupling.csv, a105_arms.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a105_covid.py (target_coupling in natexp.py, per_arm),
         rerun_loo.py (leave-one-out fix).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "natural_experiments"
ANCH = ["Endothelial cells", "Fibroblasts"]


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    organ, cl = data.load("covid_delorey", level="organ", subsample=args.subsample, seed=args.seed)
    t = S.eligible(organ, 50)
    ds = S.dataset(t, cl)  # every organ an anchor: each is read on the host of the other three
    fit = ds.fit(n_boot=run["n_boot"], n_perm=run["n_perm"], seed=5, gate=False)
    tab = fit.table.copy()
    tab.to_csv(out / "a105_delorey_organ_coupling.csv")
    print(tab[["n", "coupling", "coupling_lo", "coupling_hi", "z2", "offset_h", "psi", "boundary",
               "p_coupling"]].round(3).to_string(), flush=True)

    melms, ms = data.load("covid_melms", subsample=args.subsample, seed=args.seed)
    lung = S.eligible(melms, 20)
    ds = S.dataset(lung, ms, anchors=ANCH)
    c = ds.contrast("group", levels=("Control", "COVID-19"), blocks="none", n_perm=run["n_perm"],
                    n_boot=run["n_boot"], seed=1)
    c.table.to_csv(out / "a105_arms.csv")
    print(c.table[["anchor", "donors_Control", "donors_COVID-19", "c_Control", "c_COVID-19", "d", "p", "d_star",
                   "log_z2_ratio", "log_psi_ratio", "status"]].round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
