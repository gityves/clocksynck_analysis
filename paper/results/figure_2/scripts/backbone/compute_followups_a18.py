"""A18 part 3 (exploratory follow-ups of the A17 re-analysis): items 1 (depth thinning), 3 (PD DMV neurons) and 4 (HBCC).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A18 part 3; run note 2026-10-12/13 for the deviations). Anchored contrasts
(`Dataset.contrast`, learned anchors of the frozen A17 table, permutation p within the target's blocks).
1. `thinning`: for the fibrotic increases of stromal/vascular coupling (dcm, dcm_healthy, ild_natri, hlca_fibrotic,
   lung_ipf_hlca) the case pseudobulks of compartment k are thinned binomially to the control median depth of k (only those
   above it) and the contrast is rerun (200 permutations, 200 bootstraps); also median depth and cells per pseudobulk per group.
3. `pd`: PD DMV compartments: donors, cells and depth per group; the contrast with the pseudobulks of >= 50 cells; the
   coupling by Braak Lewy-body stage group (0, 1-4, 5-6) and 5-6 minus 0 (`path_braak_lb`).
4. `hbcc`: BD and SCZ (HBCC) against each random half of the shared controls (10 splits, 200 permutations).
Null: label permutations (A17 contrast rules). Interval: donor bootstrap (host fixed). Exploratory, as labelled in A18.

    python compute_followups_a18.py [--subsample N] [--jobs 2] [--only thinning pd hbcc]
Outputs: paper/results/figure_2/output/backbone/a18_{thinning,pd_dmv,pd_braak,hbcc}.csv
Legacy: legacy/desync/common/followups.py (thin_one, pd_dmv, pd_braak, hbcc_one).
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

FOLDER = "backbone"
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
_spec = importlib.util.spec_from_file_location("f5_shared", group_dir("resilience") / "_shared.py")
F5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(F5)
GROUPS = F5.GROUPS
THIN = {"dcm": ["heart:mural", "heart:endothelial", "heart:fibroblast"],
        "dcm_healthy": ["heart:mural", "heart:endothelial", "heart:fibroblast"],
        "ild_natri": ["lung:EC_capillary", "lung:ciliated", "lung:fibroblast"],
        "hlca_fibrotic": ["lung:fibroblast", "lung:T"], "lung_ipf_hlca": ["lung:fibroblast", "lung:multiciliated"]}
BRAAK = ("0", "1-4", "5-6")


def target(name: str, subsample, seed):
    row = F5.targets().set_index("name").loc[name]
    row.name = name
    pb, donors, lab, blocks = F5.load_contrast(row, subsample, seed)
    anc = [a for a in F5.anchor_set(name, "learned") or [] if a in set(pb.compartment)]
    return pb, donors, lab, blocks, anc


def contrast(pb, donors, lab, blocks, anc, targets, n_perm, n_boot, seed, groups="lab", levels=GROUPS):
    ds = cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab) if groups == "lab" else donors)
    return ds.contrast(groups, levels=levels, targets=targets, blocks=blocks, n_perm=n_perm, n_boot=n_boot, seed=seed).table


def thin_one(name, seed, subsample, n_perm, n_boot):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    pb, donors, lab, blocks, anc = target(name, subsample, seed)
    if len(anc) < 3:
        return []
    g = pb.donor.map(lab)
    rng = np.random.default_rng(seed)
    genes = list(cs.gtex().genes)
    rows = []
    for k in THIN[name]:
        if k not in set(pb.compartment):
            continue
        sel = pb.compartment == k
        desc = pb[sel].assign(g=g[sel]).groupby("g").agg(depth=("depth", "median"), cells=("n_cells", "median"))
        depth = desc.loc["control", "depth"]
        out = pb.copy()
        case = sel & (g == "case") & (pb.depth > depth)
        p = (depth / pb.loc[case, "depth"]).to_numpy()
        cnt = np.round(pb.loc[case, genes].to_numpy(float)).astype(np.int64)
        out.loc[case, genes] = rng.binomial(cnt, p[:, None])
        out.loc[case, "depth"] = depth
        r0 = contrast(pb, donors, lab, blocks, anc, [k], n_perm, n_boot, seed).iloc[0]
        r1 = contrast(out, donors, lab, blocks, anc, [k], n_perm, n_boot, seed).iloc[0]
        rows.append(dict(target=name, compartment=k, depth_control=depth, depth_case=desc.loc["case", "depth"],
                         cells_control=desc.loc["control", "cells"], cells_case=desc.loc["case", "cells"],
                         diff=r0["d"], p=r0.get("p", np.nan), diff_thinned=r1["d"], lo_thinned=r1.get("lo", np.nan),
                         hi_thinned=r1.get("hi", np.nan), p_thinned=r1.get("p", np.nan)))
    return rows


def braak_stage(index) -> pd.Series:
    don = pd.read_csv(DATA / "atlas" / "parkinson_regions" / "donors.csv").set_index("donor_id")
    b = pd.to_numeric(don["path_braak_lb"], errors="coerce")
    b.index = "parkinson_regions:" + b.index.astype(str)
    return pd.cut(b.reindex(index), [-0.5, 0.5, 4.5, 6.5], labels=list(BRAAK)).astype(str)


def pd_parts(seed, subsample, n_perm, n_boot):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    pb, donors, lab, blocks, anc = target("pd", subsample, seed)
    dmv = sorted(k for k in set(pb.compartment) if "DMV_" in k)
    desc = pb[pb.compartment.isin(dmv)].assign(g=pb.donor.map(lab)).groupby(["compartment", "g"]).agg(
        donors=("donor", "nunique"), cells=("n_cells", "median"), depth=("depth", "median")).reset_index()
    rows = [dict(part="describe", **r) for r in desc.to_dict("records")]
    big = pb[~(pb.compartment.isin(dmv) & (pb.n_cells < 50))].reset_index(drop=True)
    t = contrast(big, donors, lab, blocks, anc, dmv, n_perm, n_boot, seed).reset_index()
    rows += [dict(part="contrast_cells>=50", **x) for x in t.to_dict("records")]
    # Braak groups: control = stage 0 against each other stage group, same compartments
    stage = braak_stage(donors.index)
    brk = []
    for other in BRAAK[1:]:
        s = stage[stage.isin(["0", other])]
        keep = s.index
        tb = contrast(pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[keep], s.replace({"0": "control", other: "case"}),
                      "none", anc, dmv, 0, n_boot, seed).reset_index()
        for r in tb.to_dict("records"):
            brk.append(dict(compartment=r["compartment"], g=other, donors=r.get("donors_case", np.nan), coupling=r["c_case"],
                            se=r.get("se_d", np.nan), c_stage0=r["c_control"], d_vs_0=r["d"], d_lo=r.get("lo", np.nan),
                            d_hi=r.get("hi", np.nan)))
    return pd.DataFrame(rows), pd.DataFrame(brk)


def hbcc_one(name, split, seed, subsample, n_perm):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    pb, donors, lab, blocks, anc = target(name, subsample, seed)
    if len(anc) < 3:
        return pd.DataFrame()
    rng = np.random.default_rng([seed, split])
    ctrl = lab.index[lab == "control"].to_numpy()
    half = rng.permutation(ctrl)[: len(ctrl) // 2]
    out = []
    for h, keep in (("A", half), ("B", np.setdiff1d(ctrl, half))):
        l2 = pd.concat([lab[lab == "case"], lab[keep]])
        t = contrast(pb[pb.donor.isin(l2.index)].reset_index(drop=True), donors.loc[l2.index], l2, blocks, anc,
                     sorted(set(pb.compartment)), n_perm, 0, seed).reset_index()
        out.append(t.assign(target=name, split=split, half=h))
    return pd.concat(out)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=["thinning", "pd", "hbcc"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm, n_boot = (30, 10) if run["quick"] else (200, 200)
    n_splits = 2 if run["quick"] else 10
    if "thinning" in args.only:
        res = Parallel(n_jobs=args.jobs)(delayed(thin_one)(n, i, args.subsample, n_perm, n_boot) for i, n in enumerate(THIN))
        df = pd.DataFrame([r for x in res for r in x])
        df.to_csv(out / "a18_thinning.csv", index=False)
        print(df.round(3).to_string(index=False))
    if "pd" in args.only:
        d, b = pd_parts(0, args.subsample, 30 if run["quick"] else 1000, 10 if run["quick"] else 500)
        d.to_csv(out / "a18_pd_dmv.csv", index=False)
        b.to_csv(out / "a18_pd_braak.csv", index=False)
        print(b.round(3).to_string(index=False))
    if "hbcc" in args.only:
        res = Parallel(n_jobs=args.jobs)(delayed(hbcc_one)(n, s, 0, args.subsample, n_perm)
                                         for n in ("BD_hbcc", "SCZ_hbcc") for s in range(n_splits))
        df = pd.concat(res)
        df.to_csv(out / "a18_hbcc.csv", index=False)
        print(df.groupby(["target", "compartment"]).agg(d=("d", "mean"), p=("p", "mean"), p_lt05=("p", lambda x: (x < .05).mean()))
              .round(3).to_string())


if __name__ == "__main__":
    main()
