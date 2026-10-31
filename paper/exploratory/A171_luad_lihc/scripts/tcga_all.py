"""A171 follow-up (POST HOC, descriptive; TCGA is already seen in A97): the H171 statistic (within purity-proxy tertiles, high minus
low proliferation tertile) in every TCGA project of A97, with its permutation null, to choose further replication cohorts.
    python tcga_all.py [--nperm 2000] [--jobs 4]   -> ../output/tcga_all.csv"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--jobs", type=int, default=4)
a, _ = ap.parse_known_args()
sys.path.insert(0, str(HERE.parents[2] / "results/figure_5/scripts/cancer"))
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
import _shared as S  # noqa: E402
import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402
from scipy.stats import norm  # noqa: E402
import zlib  # noqa: E402

K, SEED = norm.ppf(0.95) + norm.ppf(0.80), zlib.crc32(b"A171 tcga_all") % (2 ** 31)
tertile = L.A166.tertile


def stat(Y, strata, p):
    return float(np.mean([L.coh(Y, i[tertile(p[i]) == 2], L.WF) - L.coh(Y, i[tertile(p[i]) == 0], L.WF) for i in strata]))


def one(j, project):
    t, _ = S.load_tumours(project)
    Y = np.log((t[S.GENES].to_numpy(float) + 0.5) / t.depth.to_numpy(float)[:, None])
    pt = tertile(t.purity.to_numpy()); strata = [np.flatnonzero(pt == k) for k in range(3) if (pt == k).sum() >= 30]
    p = t.prolif.to_numpy(float); T = stat(Y, strata, p)
    rng = np.random.default_rng(np.random.SeedSequence([SEED, j])); nl = np.empty(a.nperm)
    for k in range(a.nperm):
        q = p.copy()
        for i in strata:
            q[i] = rng.permutation(q[i])
        nl[k] = stat(Y, strata, q)
    return dict(project=project, n=len(t), T=T, null_mean=nl.mean(), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1),
                p_one_lower=cs.stats.perm_p(nl, T, "less", tie=0.0))


res = pd.DataFrame(Parallel(n_jobs=a.jobs)(delayed(one)(j, p) for j, p in enumerate(S.PROJECTS)))
res.to_csv(HERE.parent / "output" / "tcga_all.csv", index=False); print(res.round(4).to_string(), flush=True)
