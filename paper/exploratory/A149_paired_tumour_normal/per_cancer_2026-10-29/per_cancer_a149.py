"""A149 per cancer type, bounded statistic (DESCRIPTIVE, 2026-10-29, roadmap 7.4). No pass/fail wording.

Reuses the units, seeds, nulls and pseudo-clock sets of `bounded_2026-10-29/rescore_a149_bounded.py` (same streams, so
the pooled-all row reproduces H149.1 bounded). Per cancer and stratum (S1 malignant vs own normal, S2 beyond the stroma,
S3 malignant vs tumour stroma, S3n normal vs stroma): cohorts, patients, flagged units, the patient-weighted pooled
bounded statistic over unflagged cohorts, its permutation null (mean, q95, p), the MDE and, for S1, the pseudo-clock
companion (median, q95 over surviving sets). Not estimable: < 2 unflagged cohorts or < 15 patients in them.
    python per_cancer_a149.py [--jobs 1] [--n-perm 10000] [--n-pseudo 200]
Outputs: output/per_cancer.csv, output/per_cancer_pseudo.csv, output/strata_check.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bounded_2026-10-29"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import rescore_a149_bounded as R  # noqa: E402  (first: puts the A149 scripts on sys.path)
import a149_lib as L  # noqa: E402
import bounded_lib as B  # noqa: E402
import compute_a149 as C  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE / "output"
MIN_COH, MIN_PAT, MIN_SETS = 2, 15, 100


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    res = C.load_units()
    nl = R.nulls(res, args.n_perm)
    rows, pse_rows = [], []
    pseudo = {}
    s1_est = [x for x in nl["s1"] if x[4]]
    ps = Parallel(n_jobs=args.jobs)(delayed(R.pseudo_one)(next(r for r in res if r["study"] == x[0]), args.n_pseudo) for x in s1_est)
    pseudo = {s: a for s, a, _, _ in ps}
    print("pseudo done", flush=True)
    for k in C.KINDS:
        cancers = sorted({x[1] for x in nl[k]})
        for scope in ["ALL"] + cancers:
            allx = [x for x in nl[k] if scope == "ALL" or x[1] == scope]
            xs = [x for x in allx if x[4]]
            n_pat = sum(x[2]["n"] for x in allx)
            row = dict(kind=k, cancer=scope, n_cohorts=len(allx), n_patients=n_pat, n_unflagged_cohorts=len(xs),
                       n_unflagged_patients=sum(x[2]["n"] for x in xs), units_flagged=int(sum(B.stat(x[2]["za"], x[2]["vA"], x[2]["zb"], x[2]["vB"])["flag"] for x in allx)))
            est = scope == "ALL" or (row["n_unflagged_cohorts"] >= MIN_COH and row["n_unflagged_patients"] >= MIN_PAT)
            row["estimable"] = bool(est and xs)
            if row["estimable"]:
                w = [x[2]["n"] for x in xs]
                nul = L.pool(np.array([x[3] for x in xs]), w)
                obs = float(L.pool(np.array([[B.new_stat(x[2])] for x in xs]), w)[0])
                row.update(statistic=obs, p=B.perm_p(nul, obs), **L.summarise(nul), excess=obs - float(np.nanmean(nul)))
                if k == "s1":
                    pse = L.pool(np.array([pseudo[x[0]] for x in xs]), w)
                    ok = np.isfinite(pse)
                    row.update(pseudo_sets=int(ok.sum()))
                    if ok.sum() >= MIN_SETS:
                        row.update(pseudo_median=float(np.median(pse[ok])), pseudo_q95=float(np.quantile(pse[ok], .95)),
                                   p_pseudo=float((1 + np.sum(pse[ok] >= obs)) / (1 + ok.sum())))
            row["cohorts"] = ";".join(x[0] for x in xs)
            rows.append(row)
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "per_cancer.csv", index=False)
    print(D.drop(columns="cohorts").round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
