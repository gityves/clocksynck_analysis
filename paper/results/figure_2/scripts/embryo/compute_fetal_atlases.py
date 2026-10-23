"""A107 (primary) and A107b (primary): do the organs of a human fetus share a clock state?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A107: H107.1, H107.2; A107b: H107b.1-4).
Estimand: the model-free gate (1000 permutations; pass = p_pairing and p_acrophase < 0.05), compartments with >= 8
fetuses (the gate is reliable from 8 donors), stratum = the fetus's majority assay.
  A107  (Suo 2022, 25 fetuses, 9 organs): H107.1 whole organs (>= 1000 cells); H107.2 macrophages, one compartment per
        organ (>= 50 cells), composition-free.
  A107b (Cao 2020, 28 fetuses, sci-RNA-seq3): H107b.1 whole organs (Cerebrum, Heart, Kidney, Lung; >= 1000 cells);
        H107b.2 stromal cells per organ (Heart, Kidney, Lung, Muscle; >= 200 cells); H107b.3 within the fetal lung
        (endothelial, epithelial, stromal); H107b.4 within the fetal heart (endothelial, parenchyma, stromal).
        Holm over the four gates on max(p_pairing, p_acrophase).
A failed gate is uninformative (shallow cells, immature clocks, few fetuses), never "decoupled".

    python compute_fetal_atlases.py [--subsample N] [--seed S]
Outputs: paper/results/figure_2/output/embryo/a107_gates.csv, a107b_gates.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a107_fetal.py (H107.1-2), a107b_cao.py, a105b_organs.py
         (gate wrapper and Holm). S1, S2 and the early/late split of A107 are in natural_experiments.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "embryo"
ORGANS = ("skin", "spleen", "yolk_sac", "liver", "thymus", "gut", "bone_marrow", "kidney", "mesenteric_LN")
MIN_DON = 8
CAO = {  # test: (level, compartments, min cells)
    "H107b.1 whole organs": ("organ", ["Cerebrum", "Heart", "Kidney", "Lung"], 1000),
    "H107b.2 stromal per organ": ("class", ["Heart:stromal", "Kidney:stromal", "Lung:stromal", "Muscle:stromal"], 200),
    "H107b.3 fetal lung": ("class", ["Lung:endothelial", "Lung:epithelial", "Lung:stromal"], 200),
    "H107b.4 fetal heart": ("class", ["Heart:endothelial", "Heart:parenchyma", "Heart:stromal"], 200),
}


def gate_or_status(t, test, seed, n_perm):
    try:
        return S.run_gate(t, test, seed, n_perm)
    except ValueError as e:  # too few donors or compartments (subsamples)
        return dict(test=test, status=str(e)[:100])


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    min_don = 5 if run["quick"] else MIN_DON  # a subsample has fewer fetuses
    organ, _ = data.load("fetal_suo", level="organ_all", subsample=args.subsample, seed=args.seed)
    cls, _ = data.load("fetal_suo", level="class", subsample=args.subsample, seed=args.seed)
    whole = S.eligible(organ, 1000, min_don, compartments=ORGANS)
    mac = S.eligible(cls, 50, min_don, compartments=[f"{o}:macrophage" for o in ORGANS])
    a107 = pd.DataFrame([gate_or_status(whole, "H107.1 whole organs", 1, run["n_perm"]),
                         gate_or_status(mac, "H107.2 macrophages", 2, run["n_perm"])])
    a107.to_csv(out / "a107_gates.csv", index=False)
    print(a107.drop(columns=["compartments"], errors="ignore").round(4).to_string(index=False), flush=True)
    tabs = {lev: data.load("fetal_cao", level=lev, subsample=args.subsample, seed=args.seed)[0]
            for lev in ("organ", "class")}
    rows = []
    for j, (test, (lev, comps, mc)) in enumerate(CAO.items()):
        rows.append(gate_or_status(S.eligible(tabs[lev], mc, min_don, compartments=comps), test, 40 + j,
                                   run["n_perm"]))
    b = pd.DataFrame(rows)
    ok = b.get("p_pairing", pd.Series(dtype=float)).notna().to_numpy() if "p_pairing" in b else []
    if len(ok) and ok.any():
        b = S.holm_max(b, mask=ok)
    b.to_csv(out / "a107b_gates.csv", index=False)
    print(b.drop(columns=["compartments"], errors="ignore").round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
