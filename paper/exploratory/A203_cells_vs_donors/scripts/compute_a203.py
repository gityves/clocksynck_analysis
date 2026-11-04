"""A203: is the precision of internal-time readouts limited by cells or by donors? HLCA healthy lung, split halves.

Registered in paper/docs/registered_experiments.md (A203). Choices not in the registration: ../README.md.
Per-cell clock-gene counts are regrouped into the A17 `lung_hlca` healthy rows (acceptance: the full regrouping
reproduces the A17 counts), then
  cells arm   all donors, cells per pseudobulk capped at m (10, 25, 50, 100, 200, all), split halves;
  donors arm  all cells, donors subsampled to N (50, 100, 150, 198) within strata, split halves;
  thinning    full pseudobulks binomially thinned to the expected UMI of m cells (m = 25, 100), split halves.
Per half: `Dataset.fit(gate=False, n_boot=0)`; per readout compartment |offset_1 - offset_2| (h, circular) and |c_1 - c_2|.

    python compute_a203.py [--jobs 8] [--splits 50] [--subsample N]
Outputs: ../output/acceptance.csv, ../output/raw.csv (per arm x level x split x compartment), ../output/compartments.csv.
Public data (CELLxGENE Census); aggregates only.
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/atlas_fits"))
sys.path.insert(0, str(PAPER / "common/datasets/_atlas"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

from _shared import Target, frozen  # noqa: E402  (atlas_fits)
from lung import HLCA_TISSUE, L3  # noqa: E402  (the A17 class rule)

CELLS = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian")) / "cells/hlca_full/cells.parquet"
GENES = list(cs.gtex().genes)
M_LEVELS = (10, 25, 50, 100, 200, None)
N_LEVELS = (50, 100, 150, 198)
THIN_LEVELS = (25, 100)
MIN_DONORS_READOUT = 40
MIN_GROUP_CELLS = 10  # Census group filter of the A17 build
OUT = HERE.parent / "output"


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def regroup(pb: pd.DataFrame):
    """Cells of the A17 rows: (group index per cell, counts n_cells x 14, raw_sum per cell, rows table)."""
    cols = (["donor_id", "assay", "tissue", "cell_type", "study", "tissue_sampling_method", "ann_level_2", "ann_level_3",
             "raw_sum"] + GENES)
    c = pd.read_parquet(CELLS, columns=cols)
    # the A17 build summed the Census groups (donor x tissue x cell_type x assay x study x sampling method x ann levels
    # 2-3) with n_cells >= 10 (atlas/hlca_full/PROVENANCE.md): cells of smaller groups never entered the table
    grp = ["donor_id", "tissue", "cell_type", "assay", "study", "tissue_sampling_method", "ann_level_2", "ann_level_3"]
    c = c[c.groupby(grp, observed=True, dropna=False).raw_sum.transform("size").to_numpy() >= MIN_GROUP_CELLS]
    c["donor"] = "hlca_full:" + c.donor_id.astype(str)
    c = c[c.donor.isin(set(pb.donor))]
    cls = c.ann_level_3.astype(str).map(L3)
    cls[c.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
    c["compartment"] = "lung:" + cls.astype(str)
    assay = pb.drop_duplicates("donor").set_index("donor").stratum.str.rsplit("|", n=1).str[1]
    keep = (c.tissue_sampling_method.astype(str).isin(HLCA_TISSUE) & (c.tissue.astype(str) != "nose")
            & cls.notna() & (c.assay.astype(str).to_numpy() == assay.reindex(c.donor).to_numpy()))
    c = c[keep.to_numpy()]
    rows = pb[["donor", "compartment", "stratum"]].reset_index(drop=True)
    key = pd.Series(np.arange(len(rows)), index=pd.MultiIndex.from_frame(rows[["donor", "compartment"]]))
    g = key.reindex(pd.MultiIndex.from_frame(c[["donor", "compartment"]])).to_numpy()
    ok = ~np.isnan(g)
    return g[ok].astype(np.int64), c.loc[ok, GENES].to_numpy(np.int64), c.loc[ok, "raw_sum"].to_numpy(float), rows


def table(rows, g, X, depth, sel, n_rows):
    """Pseudobulk table of the cells `sel` (bool or index) summed per row."""
    gi = g[sel]
    cnt = np.zeros((n_rows, X.shape[1]))
    np.add.at(cnt, gi, X[sel])
    dep = np.bincount(gi, weights=depth[sel], minlength=n_rows)
    nc = np.bincount(gi, minlength=n_rows)
    t = rows.copy()
    t["n_cells"], t["depth"] = nc, dep
    t[GENES] = cnt
    return t[t.n_cells > 0]


def cap(g, m, rng):
    """Cells kept when at most m random cells per group are used (all for m None)."""
    if m is None:
        return np.ones(len(g), bool)
    o = np.lexsort((rng.random(len(g)), g))
    gs = g[o]
    start = np.r_[0, np.flatnonzero(np.diff(gs)) + 1]
    rank = np.arange(len(gs)) - np.repeat(start, np.diff(np.r_[start, len(gs)]))
    keep = np.zeros(len(g), bool)
    keep[o[rank < m]] = True
    return keep


def halves(donors: pd.Series, rng):
    """Donor -> 0 / 1, alternating within stratum after a random permutation."""
    lab = {}
    for _, idx in donors.groupby(donors).groups.items():
        p = rng.permutation(np.asarray(idx))
        for i, d in enumerate(p):
            lab[d] = i % 2
    return pd.Series(lab)


def fit(t, anchors):
    ds = cs.Dataset(t.reset_index(drop=True), anchors=[a for a in anchors if a in set(t.compartment)])
    return ds.fit(gate=False, n_boot=0).table[["offset_h", "coupling"]]


def one(arm, level, split, rows, g, X, depth, strata, anchors, readout):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(np.random.SeedSequence([203, ["cells", "donors", "thin"].index(arm),
                                                        0 if level is None else int(level), split]))  # 0 = all cells
    donors = strata
    if arm == "donors" and level < len(strata):
        frac = level / len(strata)  # within-stratum subsample, rounded per stratum, at least one donor per stratum
        pick = []
        for _, idx in strata.groupby(strata).groups.items():
            k = max(1, int(round(frac * len(idx))))
            pick += list(rng.choice(np.asarray(idx), k, replace=False))
        donors = strata.loc[pick]
    h = halves(donors, rng)
    out = []
    fits = []
    for half in (0, 1):
        dset = set(h.index[h == half])
        on = np.isin(rows.donor.to_numpy()[g], list(dset))
        if arm == "cells":
            sel = on & cap(g, level, rng)
            t = table(rows, g, X, depth, sel, len(rows))
        else:
            t = table(rows, g, X, depth, on, len(rows))
            if arm == "thin":
                f = np.minimum(1.0, level / t.n_cells.to_numpy(float))
                t[GENES] = rng.binomial(t[GENES].to_numpy(np.int64), f[:, None])
                t["depth"] = t.depth * f
        fits.append(fit(t, anchors))
    a, b = fits
    for k in readout:
        if k in a.index and k in b.index:
            do = (a.offset_h[k] - b.offset_h[k] + 12) % 24 - 12
            out.append(dict(arm=arm, level=-1 if level is None else level, split=split, compartment=k,
                            off_1=a.offset_h[k], off_2=b.offset_h[k], d_off=abs(do), c_1=a.coupling[k],
                            c_2=b.coupling[k], d_c=abs(a.coupling[k] - b.coupling[k])))
    return out


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--splits", type=int, default=50)
    args = ap.parse_args(argv)
    cli.setup(args)
    warnings.filterwarnings("ignore")
    OUT.mkdir(parents=True, exist_ok=True)

    T = Target("lung_hlca")
    pb = T.pb
    a = frozen()
    anchors = a.compartment[(a.target == "lung_hlca") & a["anchor"].astype(bool)].tolist()
    g, X, depth, rows = regroup(pb)
    full = table(rows, g, X, depth, np.ones(len(g), bool), len(rows))
    ref = pb.set_index(["donor", "compartment"])
    mine = full.set_index(["donor", "compartment"]).reindex(ref.index)
    diff = (mine[GENES] - ref[GENES]).abs()
    acc = pd.DataFrame([dict(rows_a17=len(ref), rows_regrouped=int(mine.n_cells.notna().sum()),
                             max_abs_count_diff=float(diff.max().max()),
                             rows_exact=int((diff.max(axis=1) == 0).sum()),
                             max_rel_depth_diff=float(((mine.depth - ref.depth).abs() / ref.depth).max()),
                             n_cells_equal=int((mine.n_cells == ref.n_cells).sum()), cells=len(g))])
    acc.to_csv(OUT / "acceptance.csv", index=False)
    say("acceptance\n", acc.T.to_string())
    if acc.rows_exact.iloc[0] != len(ref):
        say("ACCEPTANCE FAILED: the regrouping does not reproduce the A17 counts; no fit is run")
        sys.exit(2)

    nd = pb.groupby("compartment").donor.nunique()
    readout = sorted(nd.index[nd >= MIN_DONORS_READOUT])
    pd.DataFrame(dict(compartment=nd.index, donors=nd.values, readout=nd.index.isin(readout),
                      anchor=nd.index.isin(anchors))).to_csv(OUT / "compartments.csv", index=False)
    strata = pb.drop_duplicates("donor").set_index("donor").stratum
    if args.subsample:
        strata = strata.sample(min(args.subsample, len(strata)), random_state=args.seed)
    splits = 3 if args.subsample else args.splits
    jobs = ([("cells", m, s) for m in M_LEVELS for s in range(splits)]
            + [("donors", n, s) for n in N_LEVELS for s in range(splits)]
            + [("thin", m, s) for m in THIN_LEVELS for s in range(splits)])
    say(f"{len(jobs)} split fits (pairs of halves), {len(readout)} readout compartments, {len(strata)} donors")
    t0, res = time.time(), []
    gen = Parallel(n_jobs=args.jobs, return_as="generator")(
        delayed(one)(arm, lv, s, rows, g, X, depth, strata, anchors, readout) for arm, lv, s in jobs)
    for i, r in enumerate(gen, 1):
        res += r
        if i % 50 == 0 or i == len(jobs):
            say(f"{i}/{len(jobs)} done, {time.time() - t0:.0f} s")
    pd.DataFrame(res).to_csv(OUT / ("raw_subsample.csv" if args.subsample else "raw.csv"), index=False)


if __name__ == "__main__":
    main()
