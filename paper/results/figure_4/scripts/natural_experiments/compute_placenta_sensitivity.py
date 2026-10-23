"""A106 (secondary S1, S2): per fetal compartment cross gates with their ambient checks, the fetal-only gate and the
couplings of the fetal compartments on the maternal host. Arutyunyan 2023, first trimester.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A106: S1, S2; S3 is reported only if >= 8 pregnancies qualify, and
only 4 do, so it is not run). The primary readouts (H106.1, H106.2 and its ambient check) are in
embryo/compute_placenta.py.
  S1   per fetal compartment (EVT, Hofbauer, fibroblast) the cross gate against the maternal compartments (1000
       permutations) and its ambient check (synthetic pseudobulks, a = 10 x XIST soup share; EVT and Hofbauer carry
       2-4 % soup); coupling of every fetal compartment (VCT, SCT, EVT, Hofbauer, fibroblast) on the host of the
       maternal anchors (dS, endothelial, perivascular): `Dataset.fit` with a donor-bootstrap interval (host
       refitted) and a coupling p (fetal donors permuted against the fixed host).
  S2   gate over the fetal compartments only: does the fetal side carry a shared clock state of its own?
Units, eligibility and strata as in compute_placenta.py.

    python compute_placenta_sensitivity.py [--subsample N] [--seed S]
Outputs: paper/results/figure_4/output/natural_experiments/a106_sensitivity_gates.csv, a106_fetal_coupling.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a106_placenta.py, natexp.py (target_coupling).
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


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    pb, donors = data.load("placenta_arutyunyan", subsample=args.subsample, seed=args.seed)
    t = S.placenta_table(pb, min_donors=5 if run["quick"] else 8)
    comps = sorted(t.compartment.unique())
    mat = [c for c in comps if c in S.MAT]
    fet = [c for c in comps if c in S.FET]
    prim = [c for c in S.PRIMARY_FET if c in comps]
    share = S.soup_share(pb, donors)
    ds = S.dataset(t, donors)
    rows = []
    for k in fet:
        if k in prim:
            continue
        g = ds.gate(n_perm=run["n_perm"], seed=3, targets=[k], hosts=mat)
        rows.append(S.gate_row(g, test=f"S1 cross gate {k}", donors=g.n_donors))
        syn = S.dataset(S.synthetic_ambient(t, [k], share), donors).gate(n_perm=run["n_perm"], seed=3, targets=[k],
                                                                         hosts=mat)
        rows.append(S.gate_row(syn, test=f"S1 ambient check {k}", donors=syn.n_donors))
    if len(fet) >= 2:
        rows.append(S.run_gate(t[t.compartment.isin(fet)].reset_index(drop=True), "S2 fetal-only gate", 4,
                               run["n_perm"]))
    gates = pd.DataFrame(rows)
    gates.to_csv(out / "a106_sensitivity_gates.csv", index=False)
    print(gates.drop(columns=["compartments", "strata"], errors="ignore").round(4).to_string(index=False), flush=True)

    anchors = [a for a in S.ANCH_PLACENTA if a in comps]
    fit = S.dataset(t, donors, anchors=anchors).fit(n_boot=run["n_boot"], n_perm=run["n_perm"], seed=5, gate=False,
                                                    targets=fet)
    tab = fit.table.copy()
    tab.to_csv(out / "a106_fetal_coupling.csv")
    print(tab[["n", "coupling", "coupling_lo", "coupling_hi", "z2", "offset_h", "psi", "boundary",
               "p_coupling"]].round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
