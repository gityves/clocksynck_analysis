"""A83: is a clock observable within TCGA tumours? Gate of tumours against the same patients' normals.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A83). Data: TCGA recount3.
(a) Ten paired cancers: clock-plane share of the covariance of the 14 centred log fractions (`paired.reference_gate`,
    200 permuted-acrophase planes) in tumours and in the normals of the same patients. Estimand: share(tumour) -
    share(normal) per cancer (H83, no direction pre-specified) and its mean over the cancers. Null: patient-level
    tumour / normal label swaps (1,000 per cancer); the mean uses the same swaps jointly; two-sided.
(b) All primary tumours of each project (11, incl. OV; descriptive): gate share and p, number of projects with p < 0.05.
Caveat as pre-registered: a lower share is not proof of a broken clock (copy-number and proliferation variance lower it).

    python compute_tumour_gate.py [--subsample N] [--seed S]
Outputs (paper/results/figure_5/output/cancer/): a83_paired.csv, a83_all_tumours.csv, a83_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a83_tumour_gate.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import log_expr, reference_gate  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
SEED = 83


def logy(x: pd.DataFrame) -> np.ndarray:
    return log_expr(x[sh.GENES].astype(float), x.depth.to_numpy(float)).to_numpy()


def share(y: np.ndarray, P: np.ndarray) -> float:
    """local: clock-plane share of a sample set (the statistic `reference_gate` reports, without its null)."""
    y = y - y.mean(axis=1, keepdims=True)
    return sh.share_cov(np.cov(y.T), P)


def two_sided(null, obs) -> float:
    return float(min(1.0, 2 * min(cs.stats.perm_p(null, obs, "less", tie=0.0),
                                  cs.stats.perm_p(-np.asarray(null), -obs, "less", tie=0.0))))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_swap, n_gate = (50, 50) if run["quick"] else (1000, 200)
    rng = np.random.default_rng(SEED + args.seed)
    P = sh.plane()
    samples = sh.load_tcga(args.subsample, args.seed)
    rows, swaps = [], []
    for c in sh.CANCERS:
        case, ref = sh.pair_frames(samples, c)
        yt, yn = logy(case), logy(ref)
        gt, gn = reference_gate(yt, sh.WF, n_gate, rng), reference_gate(yn, sh.WF, n_gate, rng)
        d = gt["share"] - gn["share"]
        flips = rng.random((n_swap, len(yt))) < 0.5
        null = np.array([share(np.where(f[:, None], yn, yt), P) - share(np.where(f[:, None], yt, yn), P) for f in flips])
        swaps.append(null)
        rows.append(dict(cancer=c, n=len(yt), share_tumour=gt["share"], p_tumour=gt["p"], share_normal=gn["share"],
                         p_normal=gn["p"], diff=d, p_two=two_sided(null, d)))
        print(rows[-1], flush=True)
    a = pd.DataFrame(rows)
    obs, null = a["diff"].mean(), np.mean(swaps, axis=0)
    v = dict(mean_diff=obs, p_two=two_sided(null, obs), tumour_gated=int((a.p_tumour < 0.05).sum()),
             normal_gated=int((a.p_normal < 0.05).sum()), cancers=len(a), n_swap=n_swap, diff_lower_in=int((a["diff"] < 0).sum()))

    rows = []
    for c in sh.PROJECTS:
        tum, _ = sh.load_tumours(c, args.subsample, args.seed)
        g = reference_gate(logy(tum), sh.WF, n_gate, rng)
        rows.append(dict(project=c, tumours=len(tum), share=g["share"], share_null=g["share_null_mean"], p=g["p"]))
        print(rows[-1], flush=True)
    b = pd.DataFrame(rows)
    v["all_tumour_projects_gated"] = int((b.p < 0.05).sum())
    a.to_csv(out / "a83_paired.csv", index=False)
    b.to_csv(out / "a83_all_tumours.csv", index=False)
    pd.Series(v).to_csv(out / "a83_verdict.csv", header=False)
    print(pd.Series(v).to_string())


if __name__ == "__main__":
    main()
