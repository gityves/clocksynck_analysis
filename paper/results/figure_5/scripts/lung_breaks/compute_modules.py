"""Module tests of the coupling breaks: is the loss of the named cell shared by its module, and specific to it?

Pre-registered as the "module test of the breaks" addenda of A18 (2026-10-13; exploratory, post hoc): addendum 1 (T2D
islets, ILD Natri, HLCA fibrotic), addendum 2 (atopic-dermatitis skin, lesional and non-lesional), addendum 3
(the same statistics at a common host spread, d*). Per contrast, on the anchored host of the A17 learned rule, with
d_k = c_case - c_control of every module compartment:
  1. module mean of d_k;  2. module mean minus d of the stroma;  3. d of the named cell minus the mean d of the rest
of its module.
One label permutation shared by all compartments, within the contrast's blocks; two-sided p; compartments with a
non-estimable d (both groups at the boundary, or < 5 donors) are left out of every mean; no multiplicity correction
(3 statistics x contrasts, reported in full). Estimand: the three statistics, on the d scale and on the d* scale
(common host spread of the pooled fit). Null: exchangeable case / control labels.

    python compute_modules.py [--subsample N] [--seed S] [--jobs 4] [--only NAME ...] [--n-perm N]
Outputs (paper/results/figure_5/output/lung_breaks/): modules.csv (one row per target x scale x statistic: value, p, n_perm_ok,
d[<compartment>])
Legacy: legacy/desync/common/module_breaks.py -> common_output/a18_module_breaks{,_common_v}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import pandas as pd  # noqa: E402

from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "lung_breaks"
MODULES = {
    "islets_t2d": dict(module=["islet:alpha", "islet:beta", "islet:delta", "islet:PP"], stroma=["islet:stellate"],
                       named="islet:beta"),
    "ild_natri": dict(module=["lung:T", "lung:B", "lung:NK"], stroma=["lung:fibroblast"], named="lung:T"),
    "hlca_fibrotic": dict(module=["lung:T", "lung:B", "lung:NK_ILC"], stroma=["lung:fibroblast"], named="lung:T"),
    "skin_ad": dict(module=["skin:T", "skin:NK"], stroma=["skin:fibroblast"], named="skin:T"),
    "skin_ad_nonlesional": dict(module=["skin:T", "skin:NK"], stroma=["skin:fibroblast"], named="skin:T"),
}


def one(name, n_perm, subsample, seed, jobs):
    spec = MODULES[name]
    pb, donors, lab, blocks, row = sh.load_target(name, subsample, seed)
    anc = [a for a in sh.anchor_set(name, "learned") if a in set(pb.compartment)]
    ds = sh.dataset(pb, donors, anc, lab, name)
    comps = [k for k in dict.fromkeys(spec["module"] + spec["stroma"]) if k in set(pb.compartment)]
    mod, stroma = [k for k in spec["module"] if k in comps], [k for k in spec["stroma"] if k in comps]
    named = spec["named"]
    rest = [k for k in mod if k != named]
    stats = {"module_mean": sh.mean_of(mod),
             "module_minus_stroma": sh.difference(sh.mean_of(mod), sh.mean_of(stroma)),
             "named_minus_rest": sh.difference(lambda dk: dk.get(named, float("nan")), sh.mean_of(rest))}
    blk = sh.block_series(pb, blocks, donors, lab)
    res = sh.module_test(ds, lab, blocks, blk, stats, comps, n_perm, seed=seed, jobs=jobs)
    return res.assign(target=name, donors=len(lab), anchors=";".join(anc))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=list(MODULES))
    ap.add_argument("--n-perm", type=int, default=None, help="permutations (default 2000 as pre-registered; 50 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = args.n_perm if args.n_perm is not None else (50 if run["quick"] else 2000)
    res = pd.concat([one(n, n_perm, args.subsample, args.seed, args.jobs) for n in args.only], ignore_index=True)
    first = ["target", "scale", "statistic", "value", "p", "n_perm_ok", "donors"]
    res = res[first + [c for c in res.columns if c not in first]]
    res.to_csv(out / "modules.csv", index=False)
    pd.set_option("display.width", 200)
    print(res[first].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
