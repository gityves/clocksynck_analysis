"""A167: advanced liver fibrosis (F3-4 vs F0-2) and the stellate / endothelial clocks, living-donor controls excluded
(registered in paper/docs/registered_experiments.md, A167).

Stage `power` (run and recorded first): MDE on d* of the two primaries by `power.mde`, with the A153 conversion.
Stage `contrast`: per-arm host gate, contrast with every readout (d*, phase offset, log z2, log psi), and the liver-wide
companion (each readout minus the mean of all other estimable compartments, leave-self-out) on the same permutations,
for the primary design D2 and the secondaries (D1, leave-one-study-out). Labels are permuted within study.

    python compute_a167.py --stage power|contrast [--subsample N] [--jobs 4]
Writes output/: power.csv (power stage); contrast_<analysis>.csv, gates.csv, analyses.csv (contrast stage).
"""
from __future__ import annotations

import sys
import time
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A165_disease_readouts" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import compute_a165 as A165  # noqa: E402  (puts paper/ on sys.path through compute_a153)
from clocksync import power  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402

import clocksync as cs  # noqa: E402

OUT = HERE.parent / "output"
SEED = zlib.crc32(b"A167") % (2 ** 31)
ANCHORS = ["liver:endothelial", "liver:stellate", "liver:macrophage"]
PRIMARIES = ["liver:stellate", "liver:endothelial"]
STUDIES = ("GSE202379", "GSE289173", "GSE244832", "GSE256398")
MDE_REPS, MDE_PERM, POWER = 100, 200, 0.8


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def designs(subsample=None, seed=0) -> dict:
    """Analysis name -> case / control labels (living-donor controls excluded); D2 is the primary."""
    pb, donors = data.load("masld_pooled")
    keep = ~donors.living_donor.astype(bool)
    d1, d2 = donors.arm_d1[keep].dropna(), donors.arm_d2[keep].dropna()
    if subsample:  # per study, so that every study keeps both arms
        rng = np.random.default_rng(seed)
        ids = [i for s in STUDIES for i in rng.permutation(d2.index[donors.study.reindex(d2.index) == s])[:subsample // 4]]
        d1, d2 = d1.reindex(ids).dropna(), d2.reindex(ids).dropna()
    out = {"d2": d2, "d1": d1}
    for s in STUDIES:
        out[f"d2_loso_{s}"] = d2[donors.study.reindex(d2.index) != s]
    return pb, donors, out


def dataset(pb, donors, lab):
    lab = lab[lab.isin(("control", "case"))]
    pb = pb[pb.donor.isin(lab.index)].reset_index(drop=True)
    return cs.Dataset(pb, anchors=[a for a in ANCHORS if a in set(pb.compartment)], donors=donors.assign(lab=lab)), lab


def mde_one(pb, donors, lab, k, n_rep, n_perm):
    """A153 conversion: s* interpolated linearly on the damping grid; MDE = c0 - s*^2 c0 / (s*^2 c0 + 1 - c0)."""
    limit_threads(1)
    t0 = time.time()
    ds, lab = dataset(pb, donors, lab)
    m = power.mde(ds.study, lab, ds.anchors, k, n_rep=n_rep, n_perm=n_perm, alpha=0.05, power=POWER, seed=SEED,
                  within_stratum=True)
    cur = m["curve"].sort_values("s", ascending=False).reset_index(drop=True)
    c0, reach = float(m["c_pooled"]), cur[cur.power >= POWER]
    if reach.empty:
        s_star = np.nan
    elif reach.index[0] == 0:
        s_star = float(cur.s.iloc[0])
    else:
        a, b = cur.iloc[reach.index[0] - 1], cur.iloc[reach.index[0]]
        s_star = float(a.s + (POWER - a.power) * (b.s - a.s) / (b.power - a.power))
    mde = c0 - s_star ** 2 * c0 / (s_star ** 2 * c0 + 1 - c0) if np.isfinite(s_star) else np.inf
    return dict(compartment=k, c0=c0, s_star=s_star, mde_dstar=mde, n_rep=n_rep, n_perm=n_perm,
                donors_control=int((lab == "control").sum()), donors_case=int((lab == "case").sum()),
                curve=m["curve"].to_json(orient="records"), seconds=round(time.time() - t0))


def contrast_one(name, pb, donors, lab, run):
    limit_threads(1)
    t0 = time.time()
    ds, lab = dataset(pb, donors, lab)
    head = dict(analysis=name, donors_control=int((lab == "control").sum()), donors_case=int((lab == "case").sum()),
                anchors=";".join(ds.anchors), status="ok", reason="")
    gates = []
    g = ds.gate(n_perm=run["n_perm"], seed=SEED)
    gates.append(dict(analysis=name, arm="pooled", p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                      clock_share=g.clock_share, passes=bool(g.passes)))
    if not g.passes:
        head.update(status="not estimable", reason="pooled gate failed")
        return head, gates, None
    c = ds.contrast("lab", levels=("control", "case"), n_perm=run["n_perm"], n_boot=run["n_boot"], seed=SEED,
                    require_host_gate_each_arm=True, n_gate=run["n_perm"], anchor_readout="loo")
    for arm, gt in c.gate.items():
        gates.append(dict(analysis=name, arm=arm, **{k: gt.get(k, np.nan) for k in
                                                     ("p_pairing", "p_acrophase", "clock_share", "passes")}))
    t = c.table.reset_index()
    if t.status.astype(str).str.startswith("not estimable").all():
        head.update(status="not estimable", reason=str(t.status.iloc[0]))
        return head, gates, None
    allk = [k for k in t.compartment if k in c.null["d"].columns]
    t = t.merge(A165.companion_all(c, allk), on="compartment", how="left")
    t.insert(0, "analysis", name)
    t["is_anchor"], t["is_primary"] = t.compartment.isin(ds.anchors), t.compartment.isin(PRIMARIES)
    head["seconds"] = round(time.time() - t0)
    return head, gates, t


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "contrast"), required=True)
    args = ap.parse_args()
    run = cli.setup(args)
    out = OUT / "subsample" if run["quick"] else OUT
    out.mkdir(parents=True, exist_ok=True)
    pb, donors, des = designs(args.subsample, args.seed)
    if args.stage == "power":
        n_rep, n_perm = (10, 50) if run["quick"] else (MDE_REPS, MDE_PERM)
        say(f"power: {PRIMARIES}, D2 {des['d2'].value_counts().to_dict()}, reps {n_rep}, perms {n_perm}")
        res = Parallel(n_jobs=min(args.jobs, len(PRIMARIES)))(
            delayed(mde_one)(pb, donors, des["d2"], k, n_rep, n_perm) for k in PRIMARIES)
        pd.DataFrame(res).to_csv(out / "power.csv", index=False)
        say(pd.DataFrame(res).drop(columns="curve").to_string(index=False))
        return
    say(f"contrast: {list(des)}, n_perm {run['n_perm']}, n_boot {run['n_boot']}, jobs {args.jobs}")
    res = Parallel(n_jobs=args.jobs, return_as="generator_unordered")(
        delayed(contrast_one)(n, pb, donors, lab, run) for n, lab in des.items())
    heads, gates = [], []
    for head, g, t in res:
        heads.append(head)
        gates += g
        if t is not None:
            t.to_csv(out / f"contrast_{head['analysis']}.csv", index=False)
        say(f"  {head['analysis']}: {head['status']} {head['reason']} ({head.get('seconds', '?')} s)")
    pd.DataFrame(heads).to_csv(out / "analyses.csv", index=False)
    pd.DataFrame(gates).to_csv(out / "gates.csv", index=False)
    say("done")


if __name__ == "__main__":
    main()
