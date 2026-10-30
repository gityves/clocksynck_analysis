"""A153: corrected Figure 4 disease contrasts (per-arm host gate, untruncated c, d*, companion) and the MDE on d*.

Registered in paper/docs/registered_experiments.md (A153 and its note of 2026-10-30, PI Q15 a). Stage 1: per non-brain
target of the A17 registry, the pooled model-free gate, then `Dataset.contrast(require_host_gate_each_arm=True,
anchor_readout="loo")` on the learned (relearned A17) anchors. Stage 2: the MDE on d* (`clocksync.power.mde`, n_rep 100,
n_perm 200, default damping grid) per estimable compartment, in three tiers: (1) family fibroblast_stromal, (2)
endothelial, (3) the rest. Everything is resumable: one file per target (output/targets) and per target x compartment
(output/mde). Scoring: score_a153.py. Choices not in the registration: README.md.

    python compute_a153.py [--subsample N] [--seed S] [--jobs 8] [--tiers 1 2 3] [--only NAME ...] [--no-mde]
"""

from __future__ import annotations

import json
import sys
import threading
import time
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_4/scripts/resilience"))
sys.path.insert(0, str(PAPER / "results/extended_data/scripts/host_identity"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import power  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

import _shared as sh  # noqa: E402
from anchor_rule import family  # noqa: E402

OUT = HERE.parent / "output"
MDE_REPS, MDE_PERM, POWER = 100, 200, 0.8


def seed_of(name: str, seed: int = 0) -> int:
    return (zlib.crc32(name.encode()) + seed) % (2 ** 31)


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def _dataset(row, pb, donors, lab):
    anc = sh.anchor_set(row.name, "learned")
    if anc is None:
        return None, None
    anc = [a for a in anc if a in set(pb.compartment)]
    if len(anc) < 3:
        return None, None
    return cs.Dataset(pb, anchors=anc, donors=donors.assign(lab=lab), name=row.name), anc


def _own_companion(c, comp: list, n_boot_ok=None) -> pd.DataFrame:
    """Companion (target minus the mean of the stromal-vascular companion set) from the permutation draws of the
    contrast, one row per compartment. A compartment of the set is read against the other members (leave-self-out)."""
    rows = []
    for k in c.table.index:
        if k not in c.null["d"].columns:
            continue
        others = [x for x in comp if x != k and x in c.null["d"].columns]
        r = dict(compartment=k, companion="+".join(others), companion_self_excluded=k in comp)
        for sc, name in (("d", "d"), ("d_star", "d_star")):
            if not others or not np.isfinite(c.observed[sc][k]):
                r.update({f"{name}_minus_companion": np.nan, f"p_companion_{name}": np.nan})
                continue
            obs = c.observed[sc][k] - np.nanmean([c.observed[sc][x] for x in others])
            null = (c.null[sc][k] - c.null[sc][others].mean(axis=1, skipna=True)).to_numpy(float)
            r.update({f"{name}_minus_companion": obs, f"p_companion_{name}": cs.stats.perm_p(null, obs)})
        rows.append(r)
    return pd.DataFrame(rows)


def one_target(row, run, subsample, seed):
    """Gate and contrast of one target; written to output/targets/<name>.*; returns the status line."""
    limit_threads(1)
    f = OUT / "targets" / f"{row.name}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    sd = seed_of(row.name, seed)
    pb, donors, lab, blocks = sh.load_contrast(row, subsample, sd)
    info = dict(target=row.name, dataset=row.dataset, tissue=row.tissue, label=row.label, donors=int(lab.size),
                donors_control=int((lab == "control").sum()), donors_case=int((lab == "case").sum()), seed=sd,
                blocks=str(blocks), status="ok", reason="")
    ds, anc = _dataset(row, pb, donors, lab)
    if ds is None:
        info.update(status="not estimable", reason="no anchor (learned rule gives fewer than 3 anchors)")
    else:
        info["anchors"] = ";".join(anc)
        g = ds.gate(n_perm=run["n_perm"], seed=sd)
        info.update(p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share, gate_passed=bool(g.passes))
        if not g.passes:
            info.update(status="not estimable", reason=f"pooled gate failed (p_pairing {g.p_pairing:.3g}, "
                                                      f"p_acrophase {g.p_acrophase:.3g})")
    if info["status"] == "ok":
        fixed = sh.anchor_set(row.name, "fixed") or []
        comp = [a for a in fixed if a in set(pb.compartment)]
        c = ds.contrast("lab", levels=sh.GROUPS, companion=comp or None, blocks=blocks, n_perm=run["n_perm"],
                        n_boot=run["n_boot"], seed=sd, require_host_gate_each_arm=True, n_gate=run["n_perm"],
                        anchor_readout="loo")
        for arm, gt in c.gate.items():
            for key in ("p_pairing", "p_acrophase", "clock_share"):
                info[f"arm_{arm}_{key}"] = gt.get(key, np.nan)
            info[f"arm_{arm}_passes"] = bool(gt.get("passes", False))
        t = c.table.reset_index()
        if t.status.astype(str).str.startswith("not estimable").all():
            info.update(status="not estimable", reason=str(t.status.iloc[0]))
        else:
            t = t.merge(_own_companion(c, comp), on="compartment", how="left")
            t["companion_set"] = ";".join(comp)
            t["family"] = t.compartment.map(family)
            t.insert(0, "target", row.name)
            t["is_anchor"] = t.compartment.isin(anc)
            t.to_csv(OUT / "targets" / f"{row.name}_table.csv", index=False)
    info["seconds"] = round(time.time() - t0)
    f.write_text(json.dumps(info, default=str))
    return info


def mde_one(row, k, subsample, seed, n_rep, n_perm):
    """MDE on d* of compartment k: power curve (`power.mde`) and the registered conversion (s* interpolated linearly
    on the damping grid; MDE = c0 - s*^2 c0 / (s*^2 c0 + 1 - c0), c0 the pooled coupling at the common spread)."""
    limit_threads(1)
    f = OUT / "mde" / f"{row.name}__{k.replace(':', '_')}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    sd = seed_of(row.name, seed)
    pb, donors, lab, blocks = sh.load_contrast(row, subsample, sd)
    ds, _ = _dataset(row, pb, donors, lab)
    if blocks == "none":
        within, blk = False, None
    elif blocks is None:
        within, blk = pb.get("stratum", pd.Series(dtype=object)).nunique() >= 2, None
    else:
        within, blk = True, donors[blocks]
    out = dict(target=row.name, compartment=k)
    try:
        m = power.mde(ds.study, lab, ds.anchors, k, n_rep=n_rep, n_perm=n_perm, alpha=0.05, power=POWER, seed=sd,
                      case="case", control="control", within_stratum=bool(within), blocks=blk)
        cur = m["curve"].sort_values("s", ascending=False).reset_index(drop=True)  # s = 1 first (drop ascending)
        c0 = float(m["c_pooled"])
        reach = cur[cur.power >= POWER]
        if reach.empty:
            s_star = np.nan
        else:
            j = reach.index[0]
            if j == 0:
                s_star = float(cur.s.iloc[0])
            else:
                a, b = cur.iloc[j - 1], cur.iloc[j]
                s_star = float(a.s + (POWER - a.power) * (b.s - a.s) / (b.power - a.power))
        mde_d = (np.nan if not np.isfinite(c0) else "inf" if not np.isfinite(s_star) else
                 float(c0 - s_star ** 2 * c0 / (s_star ** 2 * c0 + 1 - c0)))
        out.update(c0=c0, s_star=None if not np.isfinite(s_star) else s_star, mde_dstar=mde_d,
                   mde_power_mde=None if not np.isfinite(m["mde"]) else float(m["mde"]),
                   curve=m["curve"].to_json(orient="records"), n_rep=n_rep, n_perm=n_perm, error="")
    except Exception as e:  # reported in the table; the run goes on
        out.update(mde_dstar=np.nan, error=f"{type(e).__name__}: {e}")
    out["seconds"] = round(time.time() - t0)
    f.write_text(json.dumps(out, default=str))
    return out


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=None, help="targets (default: all non-brain contrasts)")
    ap.add_argument("--tiers", nargs="+", type=int, default=[1, 2, 3], help="MDE tiers: 1 fibroblast_stromal, 2 endothelial, 3 rest")
    ap.add_argument("--no-mde", action="store_true")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    global OUT
    if run["quick"]:
        OUT = OUT / "subsample"
    for d in ("targets", "mde"):
        (OUT / d).mkdir(parents=True, exist_ok=True)
    n_rep, n_perm = (6, 50) if run["quick"] else (MDE_REPS, MDE_PERM)
    tg = sh.targets("contrast")
    brain = tg[tg.tissue == "brain"]
    tg = tg[tg.tissue != "brain"]
    if args.only:
        tg = tg[tg.name.isin(args.only)]
    pd.DataFrame(dict(target=brain.name, status="not eligible", reason="brain target: Extended Data (A141), not in the "
                      "A153 family; not run")).to_csv(OUT / "brain_targets.csv", index=False)
    rows = {r.name: r for r in tg.itertuples()}
    progress = dict(stage="start", done=0, total=0)

    def beat():
        while True:
            time.sleep(600)
            say(f"alive: {progress['stage']} {progress['done']}/{progress['total']}")

    threading.Thread(target=beat, daemon=True).start()

    progress.update(stage="contrast", done=0, total=len(rows))
    say(f"stage 1: {len(rows)} non-brain targets ({len(brain)} brain targets skipped); jobs {args.jobs}")
    res = []
    for r in Parallel(n_jobs=args.jobs, return_as="generator_unordered")(
            delayed(one_target)(rw, run, args.subsample, args.seed) for rw in rows.values()):
        res.append(r)
        progress["done"] += 1
        say(f"  {r['target']}: {r['status']} {r.get('reason', '')} ({r.get('seconds', '?')} s)")
    pd.DataFrame(res).sort_values("target").to_csv(OUT / "target_status.csv", index=False)
    if args.no_mde:
        return
    tabs = [pd.read_csv(OUT / "targets" / f"{n}_table.csv") for n in rows if (OUT / "targets" / f"{n}_table.csv").exists()]
    tab = pd.concat(tabs, ignore_index=True)
    est = tab[(tab.status == "ok") & tab.p_star.notna() & tab.se_d_star.notna()]
    tier = lambda fam: 1 if fam == "fibroblast_stromal" else 2 if fam == "endothelial" else 3
    for tr in args.tiers:
        todo = [(t, k) for t, k, fam in est[["target", "compartment", "family"]].itertuples(index=False) if tier(fam) == tr]
        progress.update(stage=f"mde tier {tr}", done=0, total=len(todo))
        say(f"stage 2, tier {tr}: {len(todo)} compartments (n_rep {n_rep}, n_perm {n_perm})")
        for r in Parallel(n_jobs=args.jobs, return_as="generator_unordered")(
                delayed(mde_one)(rows[t], k, args.subsample, args.seed, n_rep, n_perm) for t, k in todo):
            progress["done"] += 1
            say(f"  {r['target']} {r['compartment']}: mde_dstar {r.get('mde_dstar')} ({r.get('seconds', '?')} s) {r.get('error', '')}")
    say("done")


if __name__ == "__main__":
    main()
