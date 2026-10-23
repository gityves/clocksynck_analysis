"""A104 (GTEx fast deaths): are weak clocks pulled more tightly? Own clock amplitude vs coupling across organs.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A104). Per tissue, on the host of the other tissues of its block
(every tissue is an anchor, leave-one-out): coupling c = z2 V / (z2 V + 2 psi) and own amplitude
A = sqrt(z2 V + 2 psi). Blocks: periphery (H104.1) and brain (S1). Statistic: Spearman rho(A, c) over tissues,
one-sided permutation p (rho < 0), donor-bootstrap interval (host refitted). Aggregates only (GTEx controlled).

    python compute_oscillator_gtex.py [--subsample N] [--jobs 4]
Outputs: paper/results/figure_2/output/backbone/a104_gtex.csv, a104_gtex_tissues.csv
Legacy: legacy/desync/validation/oscillator/scripts/a104_gtex.py (S2, timed amplitude, moved to organs).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "backbone"
N_PERM = 10_000
MIN_DONORS = 40


def per_tissue(ds: cs.Dataset) -> pd.DataFrame:
    t = ds.fit(n_boot=0, gate=False).table
    return t[["coupling", "own_amplitude", "z2", "V", "psi", "boundary", "n"]]


def rho(t: pd.DataFrame) -> float:
    t = t[["own_amplitude", "coupling"]].dropna()
    return spearmanr(t.own_amplitude, t.coupling).statistic if len(t) >= 5 else np.nan


def p_negative(t: pd.DataFrame, seed: int) -> float:
    t = t[["own_amplitude", "coupling"]].dropna()
    r0 = spearmanr(t.own_amplitude, t.coupling).statistic
    rng = np.random.default_rng(seed)
    x, y = t.own_amplitude.to_numpy(), t.coupling.to_numpy()
    null = [spearmanr(rng.permutation(x), y).statistic for _ in range(N_PERM)]
    return cs.stats.perm_p(null, r0, "less", tie=0.0)


def boot_one(ds: cs.Dataset, seed: int) -> float:
    limit_threads(1)
    try:
        return rho(per_tissue(ds.resample(np.random.default_rng(seed))))
    except (ValueError, np.linalg.LinAlgError):
        return np.nan


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    pb, _ = gtex.load(min_donors=MIN_DONORS)
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb = pb[pb.donor.isin(keep)]
    rows, tabs = [], []
    for i, block in enumerate(("periphery", "brain")):
        comps = sorted(c for c in pb.compartment.unique() if c.startswith("Brain") == (block == "brain"))
        ds = cs.Dataset(pb[pb.compartment.isin(comps)], anchors=comps, name=f"gtex_{block}")
        pt = per_tissue(ds).assign(block=block)
        tabs.append(pt)
        reps = np.array(Parallel(n_jobs=args.jobs)(delayed(boot_one)(ds, 1000 * (i + 1) + j)
                                                   for j in range(run["n_boot"])), float)
        reps = reps[np.isfinite(reps)]
        ok = pt[~pt.boundary.astype(bool)]
        rows.append(dict(block=block, analysis="H104.1" if block == "periphery" else "S1",
                         tissues=int(pt[["own_amplitude", "coupling"]].notna().all(axis=1).sum()), rho=rho(pt),
                         p_one_sided=p_negative(pt, seed=i), rho_lo=np.quantile(reps, 0.025),
                         rho_hi=np.quantile(reps, 0.975), n_boot=len(reps), rho_nonboundary=rho(ok),
                         tissues_nonboundary=len(ok)))
        print(pd.DataFrame(rows).round(3).to_string(index=False), flush=True)
    pd.concat(tabs).rename_axis("tissue").to_csv(out / "a104_gtex_tissues.csv")
    pd.DataFrame(rows).to_csv(out / "a104_gtex.csv", index=False)


if __name__ == "__main__":
    main()
