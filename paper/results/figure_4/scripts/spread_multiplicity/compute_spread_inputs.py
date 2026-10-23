"""A20 / A22 / A23 inputs: per contrast and compartment the coupling difference at a common host spread, the log ratios
of slope and own noise, and the relative phase, with donor-bootstrap SEs.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A20 inputs; A22 note: own-d* detector; A23 estimand and inference).
The 22 contrasts of the A20 table (one per set of shared donors; fibrosis, inflammation, neurodegeneration,
psychiatric, metabolic) on the contrast's pooled anchored host (learned anchors; fixed rule as sensitivity):
  y1 = d* (A18 addendum 3), y2 = log(z2_case / z2_control), y3 = log(psi_case / psi_control), y4 = relative phase
  arg(z_case / z_control) in h, wrapped to (-12, 12] (y4 > 0: later relative to the anchors in cases).
Family "other" is excluded; a compartment needs >= 5 donors per group. A22 detector (note of 2026-10-13): an anchor
whose own d* has p* < 0.05 and |d*| >= 0.2 (on its leave-one-out host) is dropped from that contrast's host (kept
when fewer than 3 anchors would remain); the A23 inputs (y4) are computed on the confirmed anchors. SEs: donor
bootstrap within groups (n_boot); y4 draws are wrapped around the point estimate; MDE of y4 = 2.8 SE; rows with
SE > 4 h do not enter the pool.

    python compute_spread_inputs.py [--subsample N] [--seed S] [--jobs 4] [--only NAME ...]
Outputs (paper/results/figure_4/output/spread_multiplicity/): a20_inputs.csv, a20_inputs_fixed.csv, a22_flags.csv,
a23_inputs.csv, a23_ibd.csv (R23.3: IBDverse Crohn's ileum on its two stromal-vascular compartments)
Legacy: legacy/desync/common/{a20_inputs,common_v,a22_group_anchors,a23_phase}.py.
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
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "spread_multiplicity"
F5 = sh.F5
IBD_ANCHORS = ["ileum:stromal", "ileum:endothelial"]


def rows_of(row, anchors, pb, donors, lab, blocks, run, seed, fam, rule):
    """Contrast table (with p, p*, offset p) and bootstrap SEs of one contrast on given anchors."""
    ds = cs.Dataset(pb, anchors=anchors, donors=donors.assign(lab=lab), name=row.name)
    t = ds.contrast("lab", levels=sh.GROUPS, blocks=blocks, n_perm=run["n_perm"], n_boot=0, seed=seed).table
    keep = [k for k in t.index if fam.get(k, "other") != "other" and t.donors_control[k] >= 5 and t.donors_case[k] >= 5]
    t = t.loc[keep]
    se = sh.boot_se(pb, anchors, lab, keep, run["n_boot"], seed + 1, t)
    out = pd.DataFrame(dict(contrast=row.name, rule=rule, disease_class=sh.CLASSES[row.name],
                            cohort=sh.COHORT.get(row.name, row.name), compartment=keep,
                            family=[fam.get(k) for k in keep], anchor=[k in anchors for k in keep],
                            n_control=t.donors_control.to_numpy(), n_case=t.donors_case.to_numpy(),
                            c_control=t.c_control.to_numpy(), c_case=t.c_case.to_numpy(),
                            y1=t.d_star.to_numpy(), y2=t.log_z2_ratio.to_numpy(), y3=t.log_psi_ratio.to_numpy(),
                            y4=t.offset_diff_h.to_numpy(), p=t.p.to_numpy(), p_star=t.p_star.to_numpy(),
                            p_offset=t.p_offset.to_numpy()))
    se = se.reset_index(drop=True)
    for c in ("se1", "se2", "se3", "se4", "n_boot"):
        out[c] = se[c].to_numpy()
    out["mde4_h"] = 2.8 * out.se4
    return out


def one(row, run, subsample, seed):
    limit_threads(1)
    warnings.simplefilter("ignore", RuntimeWarning)
    pb, donors, lab, blocks = F5.load_contrast(row, subsample, seed)
    fam = F5.family_map()
    res = {}
    for rule in F5.RULES:
        anc = F5.anchor_set(row.name, rule)
        if anc is None:
            res[rule] = None
            continue
        anc = [a for a in anc if a in set(pb.compartment)]
        res[rule] = (anc, rows_of(row, anc, pb, donors, lab, blocks, run, seed, fam, rule))
    flags, a23 = [], None
    if res["learned"] is not None:
        anc, t = res["learned"]
        own = t[t.anchor]
        drop = own[(own.p_star < 0.05) & (own.y1.abs() >= 0.2)].compartment.tolist()
        for a in anc:
            r = own[own.compartment == a]
            flags.append(dict(contrast=row.name, anchor=a, d_star=r.y1.iloc[0] if len(r) else np.nan,
                              p_star=r.p_star.iloc[0] if len(r) else np.nan,
                              drop=a in drop and len(anc) - len(drop) >= 3))
        if drop and len(anc) - len(drop) >= 3:
            conf = [a for a in anc if a not in drop]
            a23 = rows_of(row, conf, pb, donors, lab, blocks, run, seed + 2, fam, "learned_confirmed")
        elif not drop:
            a23 = t.assign(rule="learned_confirmed")
    return res, flags, a23


def ibd(run, seed, subsample):
    """R23.3: IBDverse Crohn's ileum (all cases, inflamed only) vs healthy on a host of the two stromal-vascular
    compartments; y4 for every non-anchor compartment (secondary, 2-compartment host)."""
    limit_threads(1)
    fam = F5.family_map()
    out = []
    for name in ("ibd_crohn", "ibd_crohn_inflamed"):
        row = [r for r in F5.targets("contrast").itertuples() if r.name == name][0]
        pb, donors, lab, blocks = F5.load_contrast(row, subsample, seed)
        ds = cs.Dataset(pb, anchors=IBD_ANCHORS, donors=donors.assign(lab=lab), name=name)
        comps = [k for k in ds.compartments if k not in IBD_ANCHORS and fam.get(k, "other") != "other"]
        t = ds.contrast("lab", levels=sh.GROUPS, targets=comps, blocks=blocks, n_perm=run["n_perm"], n_boot=0).table
        se = sh.boot_se(pb, IBD_ANCHORS, lab, list(t.index), run["n_boot"], seed + 3, t)
        out.append(t.assign(contrast=name, se4=se.se4.reindex(t.index)).reset_index())
    return pd.concat(out, ignore_index=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=list(sh.CLASSES))
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    rows = [r for r in F5.targets("contrast").itertuples() if r.name in args.only and r.name in sh.CLASSES]
    res = Parallel(n_jobs=args.jobs)(delayed(one)(r, run, args.subsample, args.seed + i) for i, r in enumerate(rows))
    for rule, f in (("learned", "a20_inputs.csv"), ("fixed", "a20_inputs_fixed.csv")):
        parts = [r[0][rule][1] for r in res if r[0][rule] is not None]
        pd.concat(parts, ignore_index=True).to_csv(out / f, index=False)
    pd.DataFrame([x for r in res for x in r[1]]).to_csv(out / "a22_flags.csv", index=False)
    a23 = pd.concat([r[2] for r in res if r[2] is not None], ignore_index=True)
    a23.to_csv(out / "a23_inputs.csv", index=False)
    pd.concat([r[0]["learned"][1] for r in res if r[0]["learned"] is not None]).pipe(
        lambda t: print(t.groupby("contrast").size().to_string()))
    if not args.only or set(args.only) == set(sh.CLASSES):
        ibd(run, args.seed, args.subsample).to_csv(out / "a23_ibd.csv", index=False)
    ok = a23[a23.se4 <= sh.SE_MAX_H]
    print(f"A23: rows {len(a23)}, y4 estimable {a23.y4.notna().sum()}, SE <= {sh.SE_MAX_H} h: {len(ok)}; median SE "
          f"{a23.se4.median():.2f} h, median MDE {a23.mde4_h.median():.2f} h")


if __name__ == "__main__":
    main()
