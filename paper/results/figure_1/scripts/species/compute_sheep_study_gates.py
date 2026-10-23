"""A34 per-study clock gates in SheepGTEx.

Pre-registered as A34 (2026-10-16; docs/PLAN_internal_desync.md (archive)): "the clock gate per species", with the run note "the clock gate passes in sheep: clock
share 0.26 vs null 0.05; 43/45 tissues; 9/10 large studies". The all-tissue gate and the per-tissue gates (`--tissue-gates`) are in `compute_sheep.py`;
this script adds the per-study gates. Each BioProject with >= 10 animals and >= 2 tissues is gated on its own rows (a study is one stratum; the model-free gate
with the human waveform, 200 permutations). Estimand: clock share against the permuted-acrophase null, p_pairing and p_acrophase; pass = both < 0.05. A failed
gate is uninformative, not decoupled.

    python compute_sheep_study_gates.py [--subsample N] [--n-perm P]
Outputs: paper/results/figure_1/output/species/a34_gate_study.csv
Legacy: legacy/desync/validation/farm_gtex/scripts/fit_sheep.py (`study_gates`) -> output/sheep_gate_study.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "species"


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = args.n_perm or (100 if run["quick"] else 200)
    pb, donors = data.load("sheep", subsample=args.subsample, seed=args.seed)
    rows = []
    for s, sub in pb.groupby("stratum"):
        if sub.donor.nunique() < (5 if run["quick"] else 10) or sub.compartment.nunique() < 2:
            continue
        try:
            g = cs.Dataset(sub.reset_index(drop=True), anchors=None, name=f"sheep_{s}", min_donors=5).gate(n_perm=n_perm, seed=args.seed)
            rows.append(dict(stratum=s, status="ok", clock_share=g.clock_share, share_null_mean=g.share_null_mean, p_pairing=g.p_pairing,
                             p_acrophase=g.p_acrophase, n_pairs=g.n_pairs, n_animals=g.n_donors, n_tissues=g.n_compartments, passes=g.passes))
        except Exception as e:
            rows.append(dict(stratum=s, status=f"no gate: {e}"))
    R = pd.DataFrame(rows)
    R.to_csv(out / "a34_gate_study.csv", index=False)
    print(R.round(4).to_string(index=False))
    if "passes" in R:
        print(f"{int(R.passes.sum())} of {int(R.passes.notna().sum())} studies pass")


if __name__ == "__main__":
    main()
