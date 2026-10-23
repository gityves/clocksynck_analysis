"""A105 / A105b: is the clock state of fatal COVID-19 lungs and organs coherent? (Melms 2021, Delorey 2021)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A105: H105.1, H105.2, D1; A105b: H105b.1, S1).
Estimand: the model-free gate (shared clock-plane covariance of the classes of a donor against shuffled pairing and
permuted acrophases; 1000 permutations; pass = p_pairing and p_acrophase < 0.05). A failed gate is uninformative.
  H105.1   Melms COVID-19 arm (19 donors): main classes with >= 20 cells, classes with >= 10 donors (sensitivity:
           >= 50 cells). D1: the 7 surgical controls (>= 5 donors per class; procurement-confounded, descriptive).
  H105.2   Delorey whole-organ nuclei (>= 50 nuclei), organs as compartments: the gate across organs.
  H105b.1  within each organ (heart, liver, kidney, lung) of the Delorey decedents: classes with >= 20 nuclei and >= 10
           donors; Holm across the four organs on max(p_pairing, p_acrophase); lung also without the donors that are
           in Melms (D25-D32).
  S1       endothelial cells as one compartment per organ (heart, liver, lung), if >= 8 donors have >= 2 of them.

    python compute_covid_gates.py [--subsample N] [--seed S]
Outputs: paper/results/figure_4/output/natural_experiments/a105_gates.csv (H105.1, D1, H105.2), a105b_gates.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a105_covid.py (gates), a105b_organs.py.
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
ORGANS = ("heart", "liver", "kidney", "lung")
OVERLAP = {f"D{i}" for i in range(25, 33)}  # Delorey donors that are also Melms patients


def gate(pb, test, seed, n_perm, min_cells=20, min_donors=10, donors=None, **extra):
    t = S.eligible(pb, min_cells, min_donors, donors=donors)
    try:
        return S.run_gate(t, test, seed, n_perm, **extra)
    except ValueError as e:  # too few donors or classes (subsamples)
        return dict(test=test, status=str(e)[:100], **extra)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n, md = run["n_perm"], (4 if run["quick"] else 10)  # a subsample has fewer donors per class
    melms, ms = data.load("covid_melms", subsample=args.subsample, seed=args.seed)
    covid, ctrl = ms.index[ms.group != "Control"], ms.index[ms.group == "Control"]
    a105 = [gate(melms, "H105.1", 0, n, 20, md, covid, data="Melms COVID"),
            gate(melms, "D1", 0, n, 20, md // 2, ctrl, data="Melms controls (surgical)"),
            gate(melms, "H105.1 sens >=50 cells", 0, n, 50, md, covid, data="Melms COVID")]
    organ, _ = data.load("covid_delorey", level="organ", subsample=args.subsample, seed=args.seed)
    a105.append(gate(organ, "H105.2", 0, n, 50, md, data="Delorey organs"))
    a105 = pd.DataFrame(a105)
    a105.to_csv(out / "a105_gates.csv", index=False)
    print(a105.drop(columns=["compartments"], errors="ignore").round(4).to_string(index=False), flush=True)

    cls, _ = data.load("covid_delorey", level="class", subsample=args.subsample, seed=args.seed)
    rows = []
    for j, o in enumerate(ORGANS):
        t = S.eligible(cls, 20, md)
        t = t[t.compartment.str.startswith(f"{o}:")]
        rows.append(gate(t, f"H105b.1 {o}", 10 + j, n, 0, 0))
        if o == "lung":
            rows.append(gate(t[~t.donor.isin(OVERLAP)], "H105b.1 lung without Melms overlap", 20, n, 0, 0))
    b = pd.DataFrame(rows)
    ok = b.p_pairing.notna().to_numpy() & b.test.isin([f"H105b.1 {o}" for o in ORGANS]).to_numpy()
    b = S.holm_max(b, mask=ok) if ok.any() else b
    endo = S.eligible(cls, 20, 0, compartments=[f"{o}:endothelial" for o in ("heart", "liver", "lung")])
    per = endo.groupby("donor").compartment.nunique()
    if (per >= 2).sum() >= 8 or run["quick"]:
        b = pd.concat([b, pd.DataFrame([gate(endo, "S1 endothelial across organs", 30, n, 0, 0)])], ignore_index=True)
    else:
        print("S1 not run: donors with >= 2 endothelial organs:", int((per >= 2).sum()))
    b.to_csv(out / "a105b_gates.csv", index=False)
    print(b.drop(columns=["compartments"], errors="ignore").round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
