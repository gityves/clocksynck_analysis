"""A107 (secondary S1, S2, descriptive early/late): fibroblast/stromal gate, per-organ coupling, gate by gestational age.
Suo 2022, 25 fetuses.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A107: S1, S2, descriptive). The primary gates (H107.1, H107.2) are in
embryo/compute_fetal_atlases.py. Compartments with >= 8 fetuses, stratum = the fetus's majority assay.
  S1   fibroblast/stromal cells, one compartment per organ (>= 50 cells): the gate (1000 permutations).
  S2   per organ (whole organs, >= 1000 cells), coupling to the host of the other organs (`Dataset.fit`, every organ an
       anchor, leave-one-out): donor-bootstrap interval (host refitted) and coupling p (the organ's donors permuted
       against the fixed host).
  D    the H107.1 gate in early (<= median gestational age) and late fetuses (fewer than 8 each: descriptive only).

    python compute_fetal_secondary.py [--subsample N] [--seed S]
Outputs: paper/results/figure_4/output/natural_experiments/a107_secondary_gates.csv, a107_organ_coupling.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a107_fetal.py (S1, S2, early/late), rerun_loo.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402
from results.figure_2.scripts.embryo.compute_fetal_atlases import ORGANS, gate_or_status  # noqa: E402

FOLDER = "natural_experiments"


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    min_don = 5 if run["quick"] else 8
    organ, donors = data.load("fetal_suo", level="organ_all", subsample=args.subsample, seed=args.seed)
    cls, _ = data.load("fetal_suo", level="class", subsample=args.subsample, seed=args.seed)
    whole = S.eligible(organ, 1000, min_don, compartments=ORGANS)
    fib = S.eligible(cls, 50, min_don, compartments=[f"{o}:fibroblast_stromal" for o in ORGANS])
    rows = [gate_or_status(fib, "S1 fibroblast/stromal", 3, run["n_perm"])]
    ga = donors.ga_pcw.reindex(whole.donor.unique())
    med = ga.median()
    for half, dn in (("early", ga.index[ga <= med]), ("late", ga.index[ga > med])):
        rows.append(gate_or_status(whole[whole.donor.isin(dn)].reset_index(drop=True),
                                   f"D early/late: {half} (median {med} pcw)", 4, run["n_perm"]))
    g = pd.DataFrame(rows)
    g.to_csv(out / "a107_secondary_gates.csv", index=False)
    print(g.drop(columns=["compartments"], errors="ignore").round(4).to_string(index=False), flush=True)

    fit = S.dataset(whole, donors).fit(n_boot=run["n_boot"], n_perm=run["n_perm"], seed=8, gate=False)
    tab = fit.table.copy()
    tab.to_csv(out / "a107_organ_coupling.csv")
    print(tab[["n", "coupling", "coupling_lo", "coupling_hi", "z2", "offset_h", "psi", "boundary",
               "p_coupling"]].round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
