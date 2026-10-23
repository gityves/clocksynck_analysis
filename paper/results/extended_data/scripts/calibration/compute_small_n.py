"""A10 step 1: how many donors does a deep single-cell compartment need? Calibration by donor subsampling of healthy atlases.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A10, written after the A8 results, before this run; run note in
legacy/desync/common/small_n.md: no n meets the rule in both heart and kidney, the minimum stays at 15 donors).
Data: healthy heart (`heart_global`, 45 donors), healthy kidney (`kidney_sc` controls, classes with >= 20 of them, 26), healthy lung
(`lung_hlca` healthy donors). For n = 8, 10, 12, 15: 40 random donor subsets each (`--draws`); per subset the gate (200 permutations)
and the coupling of every compartment with >= 5 donors with a donor-bootstrap 95 % interval (200); every compartment is an anchor
(free host, each read on the host of the others). Reference value: the full-data point estimate.
Readouts per set and n: gate pass rate (p_acrophase < 0.05); median absolute error of coupling; coverage of the full-data coupling by
the interval; share of boundary estimates (coupling >= 0.999, or the package's boundary flag). Rule (A10): the smallest n with pass
rate >= 0.9, coverage >= 0.85 and median error <= 0.15 in both heart and kidney; otherwise 15.
The gate is run first; as in the legacy, estimates of a failing subset are kept (`force=True`) because the readouts compare them with the
full-data value.

    python compute_small_n.py [--subsample N] [--sets heart kidney lung] [--draws 40] [--sizes 8 10 12 15]
`--subsample N`: check mode (N donors per set, 2 draws, sizes 8 and 12, 5 bootstraps, 50 gate permutations).
Outputs: paper/results/extended_data/output/calibration/small_n_{draws,summary,full}.csv
Legacy: legacy/desync/common/small_n.py (common_output/small_n_summary.csv)
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "calibration"
MIN_DONORS_SET = 20


def load_set(name: str, subsample, seed):
    """(pseudobulk, donors) of the healthy part of one atlas; classes with >= 20 healthy donors (heart, kidney)."""
    ds = {"heart": "heart_global", "kidney": "kidney_sc", "lung": "lung_hlca"}[name]
    pb, d = data.load(ds, subsample=subsample, seed=seed)
    pb = pb[pb.donor.isin(d.index[d.healthy.astype(bool)])]
    if name != "lung":
        n = pb.groupby("compartment").donor.nunique()
        pb = pb[pb.compartment.isin(n.index[n >= (MIN_DONORS_SET if subsample is None else 5)])]
    return pb.reset_index(drop=True), d


def one(pb: pd.DataFrame, donors, n_boot: int, n_perm: int, seed: int):
    """Gate and bootstrap intervals of one donor subset (None when the subset cannot be fitted)."""
    cs.config.limit_threads(1)
    t = pb[pb.donor.isin(donors)]
    n = t.groupby("compartment").donor.nunique()
    t = t[t.compartment.isin(n[n >= 5].index)].reset_index(drop=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds = cs.Dataset(t, anchors=sorted(t.compartment.unique()))
            f = ds.fit(n_boot=n_boot, seed=seed, n_perm_gate=n_perm, force=True)
    except (ValueError, np.linalg.LinAlgError):
        return None
    return f.table[["coupling", "coupling_lo", "coupling_hi", "n", "boundary"]].assign(
        p_acrophase=f.gate.p_acrophase, clock_share=f.gate.clock_share).rename_axis("compartment").reset_index()


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--sets", nargs="*", default=["heart", "kidney", "lung"])
    ap.add_argument("--sizes", nargs="*", type=int, default=None)
    ap.add_argument("--draws", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args()
    opt = cli.setup(args)
    quick = opt["quick"]
    sizes = args.sizes or ([8, 12] if quick else [8, 10, 12, 15])
    draws = args.draws or (2 if quick else 40)
    n_boot, n_perm = (5, 50) if quick else (200, 200)
    n_boot, n_perm = args.n_boot or n_boot, args.n_perm or n_perm
    out = out_dir(FOLDER, quick)
    rows, summ, fulls = [], [], []
    for name in args.sets:
        pb, _ = load_set(name, args.subsample, args.seed)
        all_donors = np.array(sorted(pb.donor.unique()))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            full_fit = cs.Dataset(pb, anchors=sorted(pb.compartment.unique())).fit(n_boot=0, gate=False)
        full = full_fit.table.coupling
        fulls.append(full.rename("coupling").rename_axis("compartment").reset_index().assign(set=name, donors=len(all_donors)))
        rng = np.random.default_rng([0, ["heart", "kidney", "lung"].index(name)])
        tasks = [(n, i, rng.choice(all_donors, min(n, len(all_donors)), replace=False)) for n in sizes for i in range(draws)]
        res = Parallel(n_jobs=args.jobs)(delayed(one)(pb, d, n_boot, n_perm, 1000 * n + i) for n, i, d in tasks)
        mine = []
        for (n, i, _), r in zip(tasks, res):
            if r is not None:
                mine.append(r.assign(set=name, n=n, draw=i, full=r.compartment.map(full)))
        if not mine:
            continue
        D = pd.concat(mine)
        rows.append(D)
        for n, g in D.groupby("n"):
            gate = g.drop_duplicates("draw")
            summ.append(dict(set=name, n=n, draws=gate.draw.nunique(), gate_pass=(gate.p_acrophase < 0.05).mean(),
                             median_abs_error=(g.coupling - g.full).abs().median(),
                             coverage=((g.coupling_lo <= g.full) & (g.full <= g.coupling_hi)).mean(),
                             boundary=(g.coupling >= 0.999).mean(), boundary_flag=g.boundary.mean(), full_donors=len(all_donors)))
        print(pd.DataFrame(summ).round(3).to_string(), flush=True)
    s = pd.DataFrame(summ)
    s["meets_rule"] = (s.gate_pass >= 0.9) & (s.coverage >= 0.85) & (s.median_abs_error <= 0.15)
    pd.concat(rows).to_csv(out / "small_n_draws.csv", index=False)
    s.to_csv(out / "small_n_summary.csv", index=False)
    pd.concat(fulls).to_csv(out / "small_n_full.csv", index=False)
    print(s.round(3).to_string())
    piv = s[s.set.isin(["heart", "kidney"])].pivot(index="n", columns="set", values="meets_rule")
    ok = piv.index[piv.all(axis=1)] if piv.shape[1] == 2 else []
    print("rule (A10): smallest n meeting all criteria in heart and kidney:", int(ok.min()) if len(ok) else "none (the minimum stays at 15)")
