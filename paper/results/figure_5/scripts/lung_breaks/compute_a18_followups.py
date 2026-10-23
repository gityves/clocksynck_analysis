"""Follow-up checks of the A17 re-analysis, non-cancer parts of A18 part 3 (exploratory): depth thinning, PD dorsal-motor-vagal neurons, HBCC.

Pre-registered as A18 part 3 (2026-10-12; docs/PLAN_internal_desync.md (archive); all exploratory, written after the A17 couplings had been seen).
Parts (`--parts`):
  thinning  (3.1) fibrotic increases of stromal / vascular coupling (DCM, ILD, fibrotic HLCA, HLCA IPF): case pseudobulks of the tested
            compartment thinned binomially to the control median depth, anchored contrast (learned rule) rerun, 200 permutations;
            also median depth and cells per pseudobulk per group. Estimand: d = c(case) - c(control), before and after thinning.
  pd        (3.3) Parkinson DMV compartments: donors, cells and depth per group; the contrast restricted to pseudobulks >= 50 cells;
            coupling by Braak Lewy-body stage (0, 1-4, 5-6) on the pooled host with the contrasts 1-4 vs 0 and 5-6 vs 0 (permutation p
            within stratum, donor bootstrap).
  hbcc      (3.4) BD and SCZ (HBCC) against each of two random halves of the shared controls (10 splits), anchored contrast, 200
            permutations. The comparison with the MSSM controls and the MSSM SCZ MDE are not run (as in the legacy run note).
Null: label permutation within the contrast's blocks (`Dataset.contrast`); interval: donor bootstrap within group (host fixed).
Boundary couplings (psi = 0) are upper bounds.

    python compute_a18_followups.py [--subsample N] [--seed S] [--jobs 4] [--parts thinning pd hbcc]
Outputs (paper/results/figure_5/output/lung_breaks/): a18_thinning.csv, a18_pd_dmv.csv, a18_pd_braak.csv, a18_hbcc.csv
Legacy: legacy/desync/common/followups.py -> common_output/a18_{thinning,pd_dmv,pd_braak,hbcc}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "lung_breaks"
THIN = {"dcm": ["heart:mural", "heart:endothelial", "heart:fibroblast"],
        "dcm_healthy": ["heart:mural", "heart:endothelial", "heart:fibroblast"],
        "ild_natri": ["lung:EC_capillary", "lung:ciliated", "lung:fibroblast"],
        "hlca_fibrotic": ["lung:fibroblast", "lung:T"], "lung_ipf_hlca": ["lung:fibroblast", "lung:multiciliated"]}


def _dataset(name, subsample, seed, labels=None):
    pb, donors, lab, blocks, row = sh.load_target(name, subsample, seed, labels=labels)
    anc = [a for a in sh.anchor_set(name, "learned") if a in set(pb.compartment)]
    return pb, donors, lab, blocks, anc


def thin_one(name, subsample, seed, n_perm, n_boot):
    limit_threads(1)
    pb, donors, lab, blocks, anc = _dataset(name, subsample, seed)
    genes = list(cs.gtex().genes)
    g = pb.donor.map(lab)
    rng = np.random.default_rng([seed, 1])
    rows = []
    for k in THIN[name]:
        sel = pb.compartment == k
        if not sel.any():
            continue
        desc = pb[sel].assign(g=g[sel]).groupby("g").agg(depth=("depth", "median"), cells=("n_cells", "median"))
        target = desc.loc["control", "depth"]
        out = pb.copy()
        case = (sel & (g == "case") & (pb.depth > target)).to_numpy()
        p = (target / pb.loc[case, "depth"]).to_numpy()
        out.loc[case, genes] = rng.binomial(np.round(pb.loc[case, genes].to_numpy(float)).astype(np.int64), p[:, None]).astype(float)
        out.loc[case, "depth"] = target
        res = []
        for t in (pb, out):
            ds = sh.dataset(t, donors, anc, lab, name)
            res.append(ds.contrast("lab", levels=sh.GROUPS, targets=[k], blocks=blocks, n_perm=n_perm, n_boot=n_boot, seed=seed).table.iloc[0])
        rows.append(dict(target=name, compartment=k, depth_control=target, depth_case=desc.loc["case", "depth"],
                         cells_control=desc.loc["control", "cells"], cells_case=desc.loc["case", "cells"], diff=res[0]["d"], p=res[0]["p"],
                         diff_thinned=res[1]["d"], lo_thinned=res[1]["lo"], hi_thinned=res[1]["hi"], p_thinned=res[1]["p"]))
    return rows


def pd_parts(subsample, seed, n_perm, n_boot):
    limit_threads(1)
    pb, donors, lab, blocks, anc = _dataset("pd", subsample, seed)
    dmv = sorted(k for k in set(pb.compartment) if "DMV_" in k)
    desc = pb[pb.compartment.isin(dmv)].assign(g=pb.donor.map(lab)).groupby(["compartment", "g"]).agg(
        donors=("donor", "nunique"), cells=("n_cells", "median"), depth=("depth", "median")).reset_index()
    rows = [dict(part="describe", **r) for r in desc.to_dict("records")]
    big = pb[~(pb.compartment.isin(dmv) & (pb.n_cells < 50))].reset_index(drop=True)
    ds = sh.dataset(big, donors, anc, lab, "pd")
    t = ds.contrast("lab", levels=sh.GROUPS, targets=dmv, blocks=blocks, n_perm=n_perm, n_boot=n_boot, seed=seed).table.reset_index()
    rows += [dict(part="contrast_cells>=50", **r) for r in t.to_dict("records")]
    stage = sh.braak_stage(donors.index)
    ds = sh.dataset(pb, donors, anc, lab, "pd")
    br = []
    for s in ("1-4", "5-6"):
        sub = stage[stage.isin(["0", s])]
        if sub.nunique() < 2:
            continue
        keep = sub.index[sub.index.isin(ds.donor_ids)]
        d2 = ds.subset(donors=keep)
        try:
            c = d2.contrast(sub, levels=("0", s), targets=dmv, blocks=blocks if isinstance(blocks, str) else None, n_perm=n_perm,
                            n_boot=n_boot, seed=seed).table.reset_index()
        except ValueError as e:
            br.append(dict(contrast=f"{s} minus 0", status=f"infeasible: {e}"))
            continue
        c = c.rename(columns=lambda x: x[:-len(f"_{s}")] + "_stage" if x.endswith(f"_{s}") else x[:-2] + "_ref" if x.endswith("_0") else x)
        c.insert(0, "contrast", f"{s} minus 0")
        br += c.to_dict("records")
    return pd.DataFrame(rows), pd.DataFrame(br)


def hbcc_one(name, split, subsample, seed, n_perm):
    limit_threads(1)
    pb, donors, lab, blocks, anc = _dataset(name, subsample, seed)
    rng = np.random.default_rng([seed, split])
    ctrl = lab.index[lab == "control"].to_numpy()
    half = rng.permutation(ctrl)[: len(ctrl) // 2]
    out = []
    for h, keep in (("A", half), ("B", np.setdiff1d(ctrl, half))):
        l2 = pd.concat([lab[lab == "case"], lab[keep]])
        ds = sh.dataset(pb[pb.donor.isin(l2.index)], donors, anc, l2, name)
        t = ds.contrast(l2, levels=sh.GROUPS, blocks=blocks, n_perm=n_perm, n_boot=0, seed=seed + split).table.reset_index()
        out.append(t.assign(target=name, split=split, half=h))
    return pd.concat(out)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--parts", nargs="+", default=["thinning", "pd", "hbcc"], choices=["thinning", "pd", "hbcc"])
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_perm = args.n_perm or (50 if run["quick"] else 200)
    n_boot = args.n_boot or (20 if run["quick"] else 200)
    out = out_dir(FOLDER, run["quick"])
    pd.set_option("display.width", 220)
    if "thinning" in args.parts:
        res = Parallel(n_jobs=args.jobs)(delayed(thin_one)(n, args.subsample, args.seed, n_perm, n_boot) for n in THIN)
        df = pd.DataFrame([r for x in res for r in x])
        df.to_csv(out / "a18_thinning.csv", index=False)
        print(df.round(3).to_string(index=False))
    if "pd" in args.parts:
        a, b = pd_parts(args.subsample, args.seed, 1000 if not run["quick"] else 50, 500 if not run["quick"] else 20)
        a.to_csv(out / "a18_pd_dmv.csv", index=False)
        b.to_csv(out / "a18_pd_braak.csv", index=False)
        print(a[a.part == "contrast_cells>=50"].round(3).head(12).to_string(index=False))
        print(b.round(3).head(12).to_string(index=False))
    if "hbcc" in args.parts:
        res = Parallel(n_jobs=args.jobs)(delayed(hbcc_one)(n, s, args.subsample, args.seed, n_perm)
                                          for n in ("BD_hbcc", "SCZ_hbcc") for s in range(10 if not run["quick"] else 3))
        df = pd.concat(res, ignore_index=True)
        df.to_csv(out / "a18_hbcc.csv", index=False)
        print(df.groupby(["target", "compartment"]).agg(diff=("d", "mean"), p_lt_05=("p", lambda x: (x < 0.05).mean())).round(3))


if __name__ == "__main__":
    main()
