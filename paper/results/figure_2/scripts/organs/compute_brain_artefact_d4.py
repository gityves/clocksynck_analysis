"""A103 D4: the donor-level brain-minus-periphery slope on the ischaemic delay, also in Hardy 0 (ventilated) donors.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A103, D4: "the A103.1 donor-level slope also in Hardy 0 (ventilated), where
brain collection follows a different schedule"). Descriptive. The legacy ran D4 only if >= 4 brain tissues have >= 40 ventilated
donors (tissues with fewer are dropped, as in 01_gtex_organs); this script writes the feasibility table (brain tissues and
ventilated donors) and, if the condition holds, the slope of (brain score - periphery score) on (brain - periphery mean
ischaemic time, h) with a donor bootstrap and the r* of the Hardy 0 departures (`compute_brain_artefact` functions, timed fit
on the recorded time of death as for fast deaths). Estimand, null and interval as H103.1 (bootstrap, 200).

    python compute_brain_artefact_d4.py [--subsample N] [--n-boot 200]
Outputs (paper/results/figure_2/output/organs/): a103_d4_feasibility.csv (+ a103_d4_delay.csv, a103_d4_rstar.csv when feasible)
Legacy: legacy/desync/validation/brain_artefact/scripts/a103_gtex.py (the `tv` block).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402
import compute_brain_artefact as BA  # noqa: E402

FOLDER = "organs"
MIN_BRAIN, MIN_DONORS = 4, 40


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = args.n_boot if args.n_boot is not None else (10 if run["quick"] else 200)
    from common import gtex
    pb_all, _ = gtex.load(hardy=gtex.VENTILATOR, genes="clock", depth="sum")
    n = pb_all.groupby("compartment").donor.nunique()
    feas = pd.DataFrame(dict(tissue=n.index, brain=n.index.str.startswith("Brain"), ventilated_donors=n.values,
                             eligible=n.values >= MIN_DONORS)).sort_values(["brain", "ventilated_donors"], ascending=False)
    feas.to_csv(out / "a103_d4_feasibility.csv", index=False)
    nb = int((feas.brain & feas.eligible).sum())
    print(f"brain tissues with >= {MIN_DONORS} ventilated donors: {nb} (need {MIN_BRAIN}); "
          f"max over brain tissues {int(feas[feas.brain].ventilated_donors.max())}", flush=True)
    if nb < MIN_BRAIN:
        print("D4 not feasible (as legacy): no further output")
        return
    rng = np.random.default_rng(0)
    pb, donors = S.cohort("ventilator", MIN_DONORS, genes=list(S.CLOCK), subsample=args.subsample, seed=args.seed)
    d = donors.reindex(pb.donor)
    t = pd.DataFrame({"donor": pb.donor.to_numpy(), "compartment": pb.compartment.to_numpy(), "tod": d.tod.to_numpy(),
                      "hardy": d.hardy.to_numpy(), "isch": S.sample_ischaemic_h(pb)})
    ph, _ = BA.phase_departures(pb, t, list(S.CLOCK))
    Z = S.zwide(ph, t.donor.to_numpy(), t.compartment.to_numpy())
    pd.DataFrame([BA.delay_slope(Z, t, "D4: Hardy 0, 14 genes", rng, n_boot)]).to_csv(out / "a103_d4_delay.csv", index=False)
    pd.DataFrame([BA.block(Z, "D4: Hardy 0 (descriptive)", rng, n_boot)]).to_csv(out / "a103_d4_rstar.csv", index=False)


if __name__ == "__main__":
    main()
