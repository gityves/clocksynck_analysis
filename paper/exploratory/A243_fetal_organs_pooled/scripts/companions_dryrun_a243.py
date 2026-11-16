"""A243 design addendum: dry run of the companion gene lists (no clock statistic). Per readout of `score_a243.tests`,
the companion pool as scored (`development.log_mean`: genes measured, i.e. not missing, in every cohort of the
readout, with mean count >= 1 over its fetuses; `zero_in_a_cohort`: those of them with a zero total in one cohort) and the first
matched and responsive set drawn with the scoring seeds (`development.l1_test` streams [243, i, 1] and [243, i, 2]).

    python companions_dryrun_a243.py
Output (`../output/`): companions_dryrun.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a243 as S  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import development as dev  # noqa: E402

counts, obs = S.load()
l1, _ = S.tests(obs)
rows = []
for name, (r, *_) in l1.items():
    o = obs[r]
    x = counts.loc[o.index]
    lm = dev.log_mean(x, o.depth)
    zero = (x[lm.index].fillna(0).groupby(o.cohort).sum() == 0).any()
    row = dict(test=name, n=len(o), cohorts=o.cohort.nunique(), pool=len(lm), zero_in_a_cohort=int(zero.sum()))
    for j, kind in enumerate(("matched", "responsive"), 1):
        w = pseudo_waveforms(lm, 1, np.random.default_rng([S.SEED, j]), S.WF, kind)[0]
        row[f"{kind}_first"] = ";".join(w.genes)
    rows.append(row)
t = pd.DataFrame(rows)
t.to_csv(S.OUT / "companions_dryrun.csv", index=False)
print(t[["test", "n", "cohorts", "pool", "zero_in_a_cohort"]].to_string(index=False))
