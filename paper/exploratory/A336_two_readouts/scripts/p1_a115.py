"""A336 part 1, Figure 6: A115 D115.1 (uninvolved-minus-healthy offset of the clock states against the cue states).

NAMING: in A115 'core' means the central clock genes (the GTEx waveform) against the 'cue' gene set; here the waveform of the
clock side is swapped for the PER-free core loop, the cue set (A81b) is unchanged. Original: `A115_A116_psoriasis/scripts/score_a115.py`
(imported, not edited). Same data, seed 115, 2,000 bootstraps, 10,000 permutations; the pooled row is reproduced first. There are
no companions in A115.

    python p1_a115.py [--check]     # --check: tiny run, assertions only, nothing printed about values nor written
Output: ../output/part1/a115.csv.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A115_A116_psoriasis/scripts"))
sys.path.insert(0, str(HERE.parent))
import score_a115 as R  # noqa: E402
import _lib as L  # noqa: E402
import a336_common as H  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a115"
SCRIPT = "paper/exploratory/A336_two_readouts/scripts/p1_a115.py"
ST = pd.read_csv(PAPER / "exploratory/A115_A116_psoriasis/output/a115_d115_1.csv").query("series == 'pooled'").iloc[0]
TRAIT_H = 0.7


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def load_series():
    acc = pd.read_csv(R.OUT / "a115_acceptance.csv")
    ok = acc.groupby(["series", "disease"]).accepted.all()
    dropped = [f"{s}|{d}" for (s, d), v in ok.items() if not v]
    pso = [s for s in ("GSE13355", "GSE121212", "GSE201827", "GSE268958") if f"{s}|psoriasis" not in dropped]
    raw = {}
    for s in pso:
        d = L.read_series(R.DIR / f"{s}_psoriasis.parquet", R.genes)
        counts = bool((d.value_type == "counts").all())
        x = pd.concat([d[d.status == "nonlesional"].drop_duplicates("patient"), d[d.status == "healthy"].drop_duplicates("patient")], ignore_index=True)
        raw[s] = dict(x=x, counts=counts, grp=(x.status == "healthy").to_numpy(), cue=L.states_for(x, R.sets["cue"], counts)[0])
    return pso, raw


def clock_states(r, wf):
    x = r["x"]
    g = [y for y in wf.genes if y in x.columns]
    w = wf.subset(g)
    Lp = P.projector(x[g].astype(float), w, counts=r["counts"])
    return P.states(x[g].astype(float), Lp, w, depth=x.depth.to_numpy(float) if r["counts"] else None)


def items_for(raw, wf):
    out = {}
    for s, r in raw.items():
        bc, vc = clock_states(r, wf)
        bk, vk = r["cue"]
        out[s] = dict(grp=r["grp"], st=(bc, vc, bk, vk), n_u=int((~r["grp"]).sum()), n_h=int(r["grp"].sum()))
    return out


def pooled(keys, items, n_boot, n_perm):
    """the pooled row of the original `analyse` (stream default_rng(115): bootstraps, then permutations)."""
    rng = np.random.default_rng(R.SEED)
    base = {k: (np.where(~items[k]["grp"])[0], np.where(items[k]["grp"])[0]) for k in keys}
    obs = R.stat(items, base, keys)
    boot = np.array([R.stat(items, {k: (rng.choice(base[k][0], len(base[k][0])), rng.choice(base[k][1], len(base[k][1]))) for k in keys}, keys)[0]
                     for _ in range(n_boot)])
    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = {}
        for k in keys:
            p = rng.permutation(len(items[k]["grp"]))
            nu = len(base[k][0])
            idx[k] = (p[:nu], p[nu:])
        null[i] = R.stat(items, idx, keys)[0]
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return dict(D=obs[0], lo=lo, hi=hi, p=(1 + np.sum(null <= obs[0])) / (1 + n_perm), null_sd=null.std(ddof=1), u=obs[1], h=obs[2]), boot


def d_only(keys, items):
    base = {k: (np.where(~items[k]["grp"])[0], np.where(items[k]["grp"])[0]) for k in keys}
    return R.stat(items, base, keys)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t0 = time.time()
    n_boot, n_perm = (5, 20) if chk else (R.NB, R.NP)
    pso, raw = load_series()
    WF = H.WF
    items14 = items_for(raw, WF)
    r14, boot14 = pooled(pso, items14, n_boot, n_perm)
    log("stage A done")
    exp = dict(D=ST.D115_h, lo=ST.ci_lo, hi=ST.ci_hi, p=ST.p_one_sided_le, null_sd=ST.perm_null_sd, u=ST.uninvolved_offset_h, h=ST.healthy_offset_h)
    bad = {k: (r14[k], v) for k, v in exp.items() if not np.isclose(r14[k], v, rtol=0, atol=1e-6)}
    if not chk:
        log("reproduced" if not bad else f"NOT reproduced {bad}")
    if not chk and bad:
        H.write_rows(ENTRY, [dict(id="A115_D115.1", figure_panel="Fig 6 (psoriasis trait delay in uninvolved skin)", statistic="D115.1 uninvolved minus healthy core-minus-cue offset (h)",
                                  direction="-1", value_14_stored=ST.D115_h, value_14=r14["D"], reproduced=False,
                                  notes=f"14 genes did not reproduce a115_d115_1.csv (pooled row); no core value computed. {bad}", script=SCRIPT)])
        return
    core, trio = H.core(WF), H.trio(WF)
    itemsc = items_for(raw, core)
    rc, bootc = pooled(pso, itemsc, n_boot, n_perm)
    dd = (boot14 - bootc + 12) % 24 - 12
    dl, tl = H.drops(WF)
    wfs = dict(dl)
    vals = {}
    for lab, w in (dl[:3] if chk else dl):
        try:
            vals[lab] = d_only(pso, items_for(raw, w))
        except Exception:
            vals[lab] = np.nan
    tg = {}
    for s, r in raw.items():
        x = r["x"]
        g = [y for y in trio.genes if y in x.columns]
        try:
            yr = P.log_expr(x[g].astype(float), x.depth.to_numpy(float) if r["counts"] else None).to_numpy()
            tg[s] = P.reference_gate(yr, trio.subset(g), 200, np.random.default_rng([R.SEED, 1]))["p"]
        except Exception:
            tg[s] = np.nan
    if chk:
        assert tl in wfs and H.drops(WF, n=3)[1] == tl
        print("check ok")
        return
    assert np.isclose(vals[tl], rc["D"]), (vals[tl], rc["D"])
    pc = H.drop_percentiles(vals, tl, wfs, -1)
    conds = [rc["D"] <= -TRAIT_H, rc["p"] < 0.05]
    st = "holds" if all(conds) else "falls" if not any(conds) else "at the line"
    gate_ok = bool(all(np.isfinite(list(tg.values()))) and all(v < 0.05 for v in tg.values())) if np.isfinite(list(tg.values())).any() else "not computable"
    note = (f"core: D115.1 {rc['D']:.4f} h, 95% bootstrap [{rc['lo']:.3f}, {rc['hi']:.3f}], one-sided permutation p {rc['p']:.4g}, perm null SD {rc['null_sd']:.3f} h, uninvolved {rc['u']:.3f} h, healthy {rc['h']:.3f} h; "
            f"conditions [D <= -{TRAIT_H} h, one-sided p < 0.05] = {conds}; 'core' = the clock-side waveform (PER-free core loop), cue set unchanged, no companions in A115; "
            f"trio gate: reference gate p per series {({k: round(v, 3) for k, v in tg.items()})} (3 genes allow 6 acrophase arrangements: p cannot fall below about 1/6)")
    H.write_rows(ENTRY, [dict(id="A115_D115.1", figure_panel="Fig 6 (psoriasis trait delay in uninvolved skin)",
                              statistic="D115.1 uninvolved-minus-healthy clock-minus-cue offset (h; clock side = waveform read; 4 psoriasis series, 185 vs 138)", direction="-1",
                              value_14_stored=ST.D115_h, value_14=r14["D"], reproduced=True, value_core=rc["D"],
                              core_criterion="descriptive entry; rule used: trait component of at least about 0.7 h (D <= -0.7 h) with one-sided permutation p < 0.05",
                              core_status=st, diff_14_minus_core=r14["D"] - rc["D"], diff_lo90=float(np.quantile(dd, .05)), diff_hi90=float(np.quantile(dd, .95)),
                              diff_draws="own: the 2000 bootstrap draws of A115 (samples within group and series, seed 115), shared by both waveforms", **pc,
                              trio_gate_passes=gate_ok, trio_alone_value=np.nan, reading=H.reading(st, pc["trio_drop_below_q10"]),
                              blind="core blind", notes=note, script=SCRIPT)])
    log("done", f"{(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
