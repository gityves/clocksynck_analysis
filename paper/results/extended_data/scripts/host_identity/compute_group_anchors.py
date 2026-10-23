"""A22: anchors confirmed within each group (disease contrasts), the revised detector and its validation.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A22 and its note of 2026-10-13). Steps (`--step`, default all):
- `null` (V17.1): for each of the 22 A20 contrasts the control donors are split at random into two pseudo-groups (pseudo-
  case size = n_control x the case share, capped at half; 20 splits); the group confirmation rule is applied with the
  learned anchors. Pass: anchors dropped in <= 5 % of (split x anchor).
- `diagnostic` (Part 1): per contrast, rule (learned, fixed) and group: pair coherence among the anchors only (200
  shuffled pairings, 200 donor bootstraps within the group), e_kg = mean pair excess of anchor k with the other anchors,
  delta_k = e_case - e_control, SE = hypot; drop if delta < -1.96 SE and e_case <= 0 (or delta > 1.96 SE and
  e_control <= 0); fewer than 3 anchors left: "no anchor (A22)". Reported whether or not V17.1 passes.
- `flags` (revised detector, replaces the coherence rule for Part 2): an anchor is dropped from a contrast's host if its
  own d* (common host spread, on its leave-one-out pooled host) has p* < 0.05 and |d*| >= 0.2 (1,000 permutations within
  the contrast's blocks; `Dataset.contrast(targets=anchors)` with the learned anchors as host).
- `sim` (V17.3, positive control of the revised detector): two anchors per contrast (seeded draw), host-locked clock part
  damped to s = 0 and 0.5 in the case donors by binomial thinning (`clocksync.power.damped`), 3 replicates; the anchor's
  d* and p* (200 permutations); detected = p* < 0.05 and d* <= -0.2; rate by s and by contrast size.
Not here: V17.2 (A21 mouse colon, species) and Part 2 (A20 inputs recomputed on the group-confirmed host and
the A20 model refit, resilience / spread_multiplicity).

    python compute_group_anchors.py [--subsample N] [--jobs 4] [--step null diagnostic flags sim] [--only NAME ...]
Outputs: paper/results/extended_data/output/host_identity/a22_{null,diagnostic,flags,sim}.csv (columns of the legacy files).
Legacy: legacy/desync/common/a22_group_anchors.py.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import power  # noqa: E402
from _shared import FOLDER_S04, Target, frozen  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

CONTRASTS = ["ild_natri", "hlca_fibrotic", "dcm", "ckd_sc", "skin_ssc", "skin_ad", "myocarditis", "ia_crohn",
             "endometriosis", "AD_mssm", "AD_radc", "AD_seaad_dlpfc", "DLBD_mssm", "PD_mssm", "VaD_mssm", "pd",
             "SCZ_hbcc", "BD_hbcc", "SCZ_mssm", "islets_t2d", "DM_radc", "ASCVD_hbcc"]  # A20 table (legacy a20_inputs.CLASSES)
N_PERM, N_BOOT, Z, D_MIN = 200, 200, 1.96, 0.2


def anchor_sets(name: str) -> dict:
    a = frozen()
    a = a[a.target == name]
    out = {}
    for rule, col in (("learned", "anchor"), ("fixed", "anchor_fixed")):
        comps = a.compartment[a[col].astype(bool)].tolist()
        ok = len(comps) >= 3 and not (rule == "learned" and str(a.status.iloc[0]).startswith("no anchor"))
        out[rule] = comps if ok else None
    return out


def group_excess(T: Target, anchors, donors, seed: int, n_perm: int, n_boot: int) -> pd.DataFrame:
    """Mean pair excess (and bootstrap SE) of every anchor with the other anchors, over `donors`."""
    t = T.pb[T.pb.compartment.isin(anchors) & T.pb.donor.isin(donors)]
    empty = pd.DataFrame(dict(compartment=list(anchors), excess=np.nan, se=np.nan))
    if t.compartment.nunique() < 2:
        return empty
    try:
        ds = cs.Dataset(t, anchors=None, min_donors=3)
        _, cls = ds.coherence(n_perm=n_perm, n_boot=n_boot, seed=seed)
    except ValueError:  # no anchor pair shares >= 5 donors in this group
        return empty
    return cls.set_index("compartment").reindex(list(anchors))[["excess", "se", "n_donors"]].reset_index()


def confirm(T: Target, anchors, lab: pd.Series, seed: int, n_perm: int, n_boot: int) -> pd.DataFrame:
    """Per anchor: e in each group, delta (case - control), its SE and the drop flag (A22 rule)."""
    e = {g: group_excess(T, anchors, lab.index[lab == g], seed + j, n_perm, n_boot).set_index("compartment")
         for j, g in enumerate(("control", "case"))}
    c, k = e["control"], e["case"]
    out = pd.DataFrame(dict(e_control=c.excess, se_control=c.se, e_case=k.excess, se_case=k.se))
    out["delta"] = out.e_case - out.e_control
    out["se_delta"] = np.sqrt(out.se_case ** 2 + out.se_control ** 2)
    fall = (out.delta < -Z * out.se_delta) & (out.e_case <= 0)
    rise = (out.delta > Z * out.se_delta) & (out.e_control <= 0)
    out["drop"] = (fall | rise).fillna(False)
    return out.rename_axis("anchor").reset_index()


def _target(name, args):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    return Target(name, args.subsample, args.seed)


def _null_one(name, split, seed, args, n_perm, n_boot):
    T = _target(name, args)
    anc = anchor_sets(name)["learned"]
    if anc is None:
        return None
    lab = T.labels
    ctrl = lab.index[lab == "control"].to_numpy()
    n1 = int(min(round(len(ctrl) * (lab == "case").mean()), len(ctrl) // 2))
    rng = np.random.default_rng([seed, split])
    pick = set(rng.choice(ctrl, n1, replace=False))
    pseudo = pd.Series(np.where(np.isin(ctrl, list(pick)), "case", "control"), index=ctrl)
    return confirm(T, anc, pseudo, int(rng.integers(1 << 30)), n_perm, n_boot).assign(contrast=name, split=split)


def _diag_one(name, rule, seed, args, n_perm, n_boot):
    T = _target(name, args)
    anc = anchor_sets(name)[rule]
    if anc is None:
        return pd.DataFrame([dict(contrast=name, rule=rule, status="no anchor")])
    r = confirm(T, anc, T.labels, seed, n_perm, n_boot)
    kept = [a for a, d in zip(r.anchor, r["drop"]) if not d]
    return r.assign(contrast=name, rule=rule, n_anchors=len(anc), n_kept=len(kept),
                    status="ok" if len(kept) >= 3 else "no anchor (A22)")


def _flags_one(name, args, n_perm):
    T = _target(name, args)
    anc = anchor_sets(name)["learned"] or []
    if not anc:
        return []
    c = T.dataset(anchors=anc).contrast(T.labels, levels=("control", "case"), targets=anc, n_perm=n_perm, n_boot=0,
                                        blocks=T.blocks, seed=0).table
    f = ((c.p_star < 0.05) & (c.d_star.abs() >= D_MIN)).reindex(anc).fillna(False)
    kept = int((~f).sum())
    return [dict(contrast=name, anchor=a, c_control=c.get("c_control", pd.Series()).get(a),
                 c_case=c.get("c_case", pd.Series()).get(a), d_star=c.d_star.get(a), p_star=c.p_star.get(a),
                 drop=bool(f[a]), n_anchors=len(anc), n_kept=kept,
                 status="ok" if kept >= 3 else "no anchor (A22)") for a in anc]


def _sim_one(name, k, s, rep, i, args, n_perm):
    T = _target(name, args)
    anc = anchor_sets(name)["learned"]
    ds = T.dataset(anchors=anc)
    rng = np.random.default_rng([i, int(100 * s), rep])
    try:
        tab = power.damped(ds.study, T.labels, anc, k, s, rng, case="case")
        if tab is None:
            return dict(contrast=name, anchor=k, s=s, rep=rep, d_star=np.nan, p_star=np.nan)
        ds2 = cs.Dataset(tab.rename(columns={"resection": "donor"}), anchors=anc, donors=T.donors)
        c = ds2.contrast(T.labels, levels=("control", "case"), targets=[k], n_perm=n_perm, n_boot=0, blocks=T.blocks,
                         seed=int(rng.integers(1 << 30))).table.loc[k]
    except (ValueError, np.linalg.LinAlgError) as e:
        return dict(contrast=name, anchor=k, s=s, rep=rep, d_star=np.nan, p_star=np.nan, error=str(e)[:80])
    return dict(contrast=name, anchor=k, s=s, rep=rep, d_star=c.d_star, p_star=c.p_star,
                n_control=c.get("donors_control"), n_case=c.get("donors_case"))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--step", nargs="*", default=["null", "diagnostic", "flags", "sim"])
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER_S04, run["quick"])
    names = [n for n in CONTRASTS if args.only is None or n in args.only]
    n_perm, n_boot = (30, 20) if run["quick"] else (N_PERM, N_BOOT)
    n_split = 2 if run["quick"] else 20
    n_perm_flag = 50 if run["quick"] else 1000
    n_perm_sim = 30 if run["quick"] else 200
    if "null" in args.step:
        res = Parallel(n_jobs=args.jobs)(delayed(_null_one)(n, s, 1000 + CONTRASTS.index(n), args, n_perm, n_boot)
                                         for n in names for s in range(n_split))
        t = pd.concat([r for r in res if r is not None], ignore_index=True)
        t.to_csv(out / "a22_null.csv", index=False)
        print(f"V17.1 dropped share (split x anchor): {t['drop'].mean():.3f} over {len(t)}; pass <= 0.05", flush=True)
    if "diagnostic" in args.step:
        jobs = [(n, r) for n in names for r in ("learned", "fixed")]
        res = Parallel(n_jobs=args.jobs)(delayed(_diag_one)(n, r, 2000 + CONTRASTS.index(n) * 2 + j % 2, args, n_perm,
                                                           n_boot) for j, (n, r) in enumerate(jobs))
        t = pd.concat(res, ignore_index=True)
        t.to_csv(out / "a22_diagnostic.csv", index=False)
        print(t.groupby(["rule", "status"]).contrast.nunique().to_string(), flush=True)
    if "flags" in args.step:
        res = Parallel(n_jobs=args.jobs)(delayed(_flags_one)(n, args, n_perm_flag) for n in names)
        t = pd.DataFrame([r for x in res for r in x])
        t.to_csv(out / "a22_flags.csv", index=False)
        print(t[t["drop"]].round(3).to_string(index=False), flush=True)
    if "sim" in args.step:
        jobs = []
        for i, n in enumerate(CONTRASTS):
            anc = anchor_sets(n)["learned"]
            if anc is None or n not in names:
                continue
            pick = np.random.default_rng(4000 + i).choice(anc, min(2, len(anc)), replace=False)
            jobs += [(n, k, s, r, i) for k in pick for s in (0.0, 0.5) for r in range(3)]
        t = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(_sim_one)(*j, args, n_perm_sim) for j in jobs))
        t["detected"] = (t.p_star < 0.05) & (t.d_star <= -D_MIN)
        t.to_csv(out / "a22_sim.csv", index=False)
        print(t.groupby("s").detected.mean().round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
