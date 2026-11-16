"""Q36 option (c), null only: is the q95 rule calibrated when the class statistic is the mean |r|^2 over pairs
(E|r|^2 = 1 / df whatever the shape of the states) instead of the median |r|? Same atlases, pairs, companions and
placebo as A276 (500 draws, seeds [276, 3, d]); no statistic on the real pairing. Not a registered test: input to Q36.

    python option_c_null.py [--draws 500] [--jobs 6]
Output (../output/option_c/): tests.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import calibrate_a276 as C  # noqa: E402
import corrected_a277 as A277  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402


def class_means_sq(absr, meta):
    return {key: np.nanmean(absr[:, idx] ** 2, 1) for key, idx in meta.groupby(["atlas", "klass"]).indices.items()}


def draw(seed, data, groups):
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    cm = {}
    for study, atl in data.items():
        rows, metas = [], []
        for n, (key, S, meta) in atl.items():
            rows.append(C.pair_abs(S, C.A272.pair_index(C.A272.placebo(key, rng))))
            metas.append(meta)
        cm[study] = class_means_sq(np.concatenate(rows, 1), pd.concat(metas, ignore_index=True))
    return C.statistics(cm, groups)


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--draws", type=int, default=500)
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(C.AMENDMENT) / "option_c"
    out.mkdir(parents=True, exist_ok=True)
    C.A150.ATLASES.update(C.A274.ALL)
    data, groups = A277.load_all(200)
    res = Parallel(n_jobs=args.jobs)(delayed(draw)([C.SEED, 3, d], data, groups) for d in range(args.draws))
    rows = []
    for t in res[0]:
        v = np.array([r[t] for r in res])
        for kind, sl in (("matched", slice(1, 201)), ("responsive", slice(201, 401))):
            rej = [bool(x[0] > np.nanquantile(x[sl], 0.95)) for x in v]
            pct = [float(np.mean(x[sl][np.isfinite(x[sl])] < x[0])) for x in v]
            rows.append(dict(test=t, family=kind, draws=len(v), rejection_rate=float(np.mean(rej)),
                             mean_percentile=float(np.mean(pct))))
    T = pd.DataFrame(rows)
    T.to_csv(out / "tests.csv", index=False)
    print(T.round(3).to_string(), "\nmean rate", round(T.rejection_rate.mean(), 4))


if __name__ == "__main__":
    main()
