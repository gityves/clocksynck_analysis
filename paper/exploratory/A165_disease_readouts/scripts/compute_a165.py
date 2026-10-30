"""A165: phase, amplitude and noise readouts of the A153 disease contrasts (registered in
paper/docs/registered_experiments.md, A165).

Reruns the A153 stage 1 (`A153_corrected_contrasts/scripts/compute_a153.py`: same targets, tables, learned anchors,
seeds, blocks, per-arm host gate, anchor_readout="loo", 1,000 permutations, 200 bootstrap draws) with the current
`Dataset.contrast`, which adds the permutation p of the log z2 and log psi ratios (p_log_z2, p_log_psi). The
stromal-vascular companion (leave-self-out) is computed for every readout on the same permutation draws. No MDE stage:
A165 reads the MDE from the bootstrap SE (score_a165.py).

    python compute_a165.py [--subsample N] [--jobs 8] [--only NAME ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A153_corrected_contrasts" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.dataset import NULL_READOUTS  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import compute_a153 as C  # noqa: E402  (puts paper/ on sys.path)
from common import cli  # noqa: E402

OUT = HERE.parent / "output"


def companion_all(c, comp: list, n_boot_ok=None) -> pd.DataFrame:
    """Every readout of NULL_READOUTS minus the mean of the companion set (leave-self-out), on the permutation draws of
    the contrast. Offsets: the circular difference to the companion mean offset."""
    rows = []
    for k in c.table.index:
        if k not in c.null["d"].columns:
            continue
        others = [x for x in comp if x != k and x in c.null["d"].columns]
        r = dict(compartment=k, companion="+".join(others), companion_self_excluded=k in comp)
        for sc in NULL_READOUTS:
            if not others or not np.isfinite(c.observed[sc][k]):
                r.update({f"{sc}_minus_companion": np.nan, f"p_companion_{sc}": np.nan})
                continue
            if sc == "offset_diff_h":
                mean_h = lambda x: np.angle(np.nanmean(np.exp(1j * np.asarray(x, float) / HOURS), axis=-1)) * HOURS
                obs = cs.stats.wrap(c.observed[sc][k] - mean_h([c.observed[sc][x] for x in others]))
                null = cs.stats.wrap(c.null[sc][k].to_numpy(float) - mean_h(c.null[sc][others].to_numpy(float)))
            else:
                obs = c.observed[sc][k] - np.nanmean([c.observed[sc][x] for x in others])
                null = (c.null[sc][k] - c.null[sc][others].mean(axis=1, skipna=True)).to_numpy(float)
            r.update({f"{sc}_minus_companion": obs,
                      f"p_companion_{sc}": cs.stats.perm_p(null, obs, tie=0.0 if sc == "offset_diff_h" else cs.stats.TIE)})
        rows.append(r)
    return pd.DataFrame(rows)


def one(row, run, subsample, seed, out):
    """A153 stage 1 for one target, written under A165's output with the all-readout companion (set in the worker:
    joblib workers import compute_a153 afresh)."""
    C.OUT, C._own_companion = out, companion_all
    return C.one_target(row, run, subsample, seed)


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=None)
    args = ap.parse_args()
    run = cli.setup(args)
    out = OUT / "subsample" if run["quick"] else OUT
    (out / "targets").mkdir(parents=True, exist_ok=True)
    tg = C.sh.targets("contrast")
    tg = tg[tg.tissue != "brain"]
    if args.only:
        tg = tg[tg.name.isin(args.only)]
    C.say(f"{len(tg)} non-brain targets, jobs {args.jobs}, n_perm {run['n_perm']}, n_boot {run['n_boot']}")
    res = []
    for r in Parallel(n_jobs=args.jobs, return_as="generator_unordered")(
            delayed(one)(rw, run, args.subsample, args.seed, out) for rw in tg.itertuples()):
        res.append(r)
        C.say(f"  {r['target']}: {r['status']} {r.get('reason', '')} ({r.get('seconds', '?')} s)")
    pd.DataFrame(res).sort_values("target").to_csv(out / "target_status.csv", index=False)
    C.say("done")


if __name__ == "__main__":
    main()
