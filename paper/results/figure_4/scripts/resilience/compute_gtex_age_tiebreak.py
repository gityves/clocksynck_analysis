"""A94 H94.1, post hoc: sensitivity of the observed statistic to the random tie-break of the age tertiles.

Ages are integers, so tertile membership at the cut points is set by the random tie-break of A94 (`tertiles`). The observed
H94.1 (oldest minus youngest tertile, log mean radius, median over tissues at equal n) is recomputed for `--draws`
independent tie-breaks (fast deaths, ventilated, and their mean). Estimand: the spread of the observed statistic over
tie-breaks (no null, no p). Outputs: `a94_tiebreak_draws.csv` (one row per draw), `a94_tiebreak_summary.csv`. Aggregates only.
Reason: the full A94 run here (H94.1 0.122) and the legacy run (0.090) used different tie-break draws. Legacy script:
none (post hoc, 2026-10-22).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import compute_gtex_age as A  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402


def main(argv=None):
    p = cli.parser(__doc__)
    p.add_argument("--draws", type=int, default=100)
    args = p.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(A.FOLDER, run["quick"])
    s = A.load_states(A.BODY, args.subsample, args.seed).dropna(subset=["tod"])
    C = {"fast": A.Cohort(s[s.hardy.isin([1, 2])], True), "vent": A.Cohort(s[s.hardy == 0], False)}
    rows = []
    for d in range(args.draws):
        rng = np.random.default_rng(1000 + d)
        r = {"draw": d}
        for k, c in C.items():
            dn = c.donors.dropna(subset=["age"])
            r[k] = A.contrast(c, A.tertiles(dn.age, rng), 2, 0, A.radius, A.BODY, 94)
        r["pooled"] = (r["fast"] + r["vent"]) / 2
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(out / "a94_tiebreak_draws.csv", index=False)
    q = df[["fast", "vent", "pooled"]].describe(percentiles=[.05, .95]).T
    q.to_csv(out / "a94_tiebreak_summary.csv")
    print(q.round(3))


if __name__ == "__main__":
    main()
