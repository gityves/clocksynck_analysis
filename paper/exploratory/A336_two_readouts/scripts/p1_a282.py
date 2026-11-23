"""A336 part 1 wrapper: A282 pooled living D, living T_SV and post-mortem D (arm donors only, A283 form) on the 14 genes,
the core loop, the trio contribution, the gene-drop null and the trio alone.

Original: A282_pooled_living/scripts/score_a282.py (imported, not edited; its `arm_data` and `draw` are reused, seeds
[282, arm, d]). Same stacked-state design as p1_a283.py (ct_lib.py; the atlas loader of p1_a283 is reused). A282 saved no
null, so the placebo null (5000 draws per arm) is redrawn for BOTH the 14-gene and the core projector in one pass (the
draw only permutes donors). Core p_cal are informative only if V283 passed on the core projector (output/part1/a283.csv).
Gene drops: 100 random + trio (clock only). Trio contribution: 200 pair-level donor bootstraps within stratum, rng [336, 1].
Gate: trio gate (1000 permutations) on the arm's donors in every atlas entering the statistic.

    python p1_a282.py [--check] [--draws 5000] [--jobs 4]
Output: ../output/part1/a282.csv
"""
import argparse
import importlib.util
import sys

import ct_lib as L
from ct_lib import *  # noqa: F401,F403
import a336_common as C
import p1_a283 as P3

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a282.py"
spec = importlib.util.spec_from_file_location("p1_a282_orig", L.PAPER / "exploratory/A282_pooled_living/scripts/score_a282.py")
M = importlib.util.module_from_spec(spec)
sys.modules["p1_a282_orig"] = M
spec.loader.exec_module(M)
CAL, A272, A150 = P3.CAL, P3.A272, P3.A150
V283_ROW = L.EXP / "A336_two_readouts/output/part1/a283.csv"
STORED = L.PAPER / "exploratory/A282_pooled_living/output/tests.csv"
ROWS = {"A282_living_D": ("living", "D", "pooled D = T_SV - T_LY, living atlases"),
        "A282_living_T_SV": ("living", "T_SV", "pooled T_SV, living atlases"),
        "A282_postmortem_D": ("post_mortem", "D", "pooled D = T_SV - T_LY, post-mortem atlases")}
PI = P3.pi


def arm_groups(a):
    body = a[~a.atlas.isin(list(CAL.A274.CORTEX))]
    return (list(body.atlas[body.SV_X >= M.MIN_PAIRS]), list(body.atlas[(body.SV_X >= M.MIN_PAIRS) & (body.LY_X >= M.MIN_PAIRS)]),
            list(a.atlas[a.atlas.isin(list(CAL.A274.CORTEX)) & (a.CX_diff >= M.MIN_PAIRS)]))


def null_task(task):
    arm, d = task
    cm = M.draw([M.SEED, M.ARMS[arm], d], G["d"][arm])
    st = M.statistics(cm, *G["grp"][arm])
    return {k: {f: float(PI(v, sl)[0]) for f, sl in FAMILIES.items()} for k, v in st.items()}


def boot_arm(d, grp, rng):
    cm = {}
    for atl in d.values():
        for n, (key, S, meta) in atl.items():
            pairs = A272.pair_index(key)
            cm.update({(n, k): v for (_, k), v in CAL.class_medians(L.boot_abs(S, pairs, rng), meta).items()})
    return M.statistics(cm, *grp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    data, aux, groups = P3.load14(2 if chk else None)
    cls = M.source_classes()
    ref = pd.read_csv(STORED)

    def arm_inputs(dat):
        d, grp = {}, {}
        for arm in M.ARMS:
            d[arm], tab = M.arm_data(dat, cls, arm)
            grp[arm] = arm_groups(tab)
        return d, grp

    # ---- 14 genes: observed values, null redrawn, reproduction (before any core state exists) -----------------------
    d14, grp = arm_inputs(data)
    G["grp"], G["d"] = grp, d14
    obs14 = {arm: M.statistics(M.draw(None, d14[arm]), *grp[arm]) for arm in M.ARMS}
    not chk and log("arm groups", grp)
    ndraw = 10 if chk else a.draws
    tasks = [(arm, k) for arm in M.ARMS for k in range(ndraw)]

    def split(res):
        return {arm: res[i * ndraw:(i + 1) * ndraw] for i, arm in enumerate(M.ARMS)}

    t0 = time.time()
    null14 = split(L.run_pool(null_task, tasks, a.jobs))
    log(f"14-gene null: {len(tasks)} draws in {(time.time() - t0) / 60:.1f} min")

    def pcal(null, arm, stat, fam, o):
        nl = np.array([r[stat][fam] for r in null[arm]])
        return float((1 + np.sum(nl >= o)) / (1 + len(nl)))

    bad = []
    if not chk:
        for rid, (arm, stat, _) in ROWS.items():
            v = obs14[arm][stat]
            row = ref[(ref.arm == arm) & (ref.statistic == stat)].iloc[0]
            chk_ = [("clock", v[0], row["clock"])]
            for f, sl in FAMILIES.items():
                o = float(PI(v, sl)[0])
                chk_ += [(f"{f}_median", float(np.nanmedian(v[sl])), row[f"{f}_median"]), (f"{f}_pi", o, row[f"{f}_pi"]),
                         (f"{f}_p_cal", pcal(null14, arm, stat, f, o), row[f"{f}_p_cal"])]
            bad += [(rid, nm, g, w) for nm, g, w in chk_ if not np.isclose(g, w, atol=1e-6, rtol=0)]
        log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        C.write_rows("a282", [dict(id=r, figure_panel="Fig 2/S", statistic=ROWS[r][2], reproduced=False,
                                   notes=f"14-gene values not reproduced (tolerance 1e-6): {bad}", script=SCRIPT) for r in ROWS])
        return

    # ---- core, drops, trio states, core null ---------------------------------------------------------------------
    dl, tl = C.drops(n=100)
    if chk:
        dl = dl[:3] + [x for x in dl if x[0] == tl]
    lay = L.layout(len(dl))
    for study, atl in data.items():
        for n, (key, S, meta) in list(atl.items()):
            atl[n] = (key, np.vstack([S, L.stack_states(aux[(study, n)]["pb"], aux[(study, n)]["wfs"], dl, P3.STATEFN)]), meta)
    d, grp = arm_inputs(data)
    G["grp"] = grp
    G["d"] = {arm: {s: {n: (k, S[IC], m) for n, (k, S, m) in atl.items()} for s, atl in d[arm].items()} for arm in M.ARMS}
    t0 = time.time()
    nullc = split(L.run_pool(null_task, tasks, a.jobs))
    log(f"core null: {len(tasks)} draws in {(time.time() - t0) / 60:.1f} min")
    obs = {arm: M.statistics(M.draw(None, d[arm]), *grp[arm]) for arm in M.ARMS}

    # ---- trio contribution and gates ---------------------------------------------------------------------------------
    brng = np.random.default_rng([336, 1])
    d2 = {arm: {s: {n: (k, S[[0, IC.start]], m) for n, (k, S, m) in atl.items()} for s, atl in d[arm].items()} for arm in M.ARMS}
    diffs = {arm: {} for arm in M.ARMS}
    for b in range(3 if chk else 200):
        for arm in M.ARMS:
            v = boot_arm(d2[arm], grp[arm], brng)
            for k in v:
                diffs[arm].setdefault(k, []).append(v[k][0] - v[k][1])
    gates = {}
    for rid, (arm, stat, _) in ROWS.items():
        sv, both, cx = grp[arm]
        names = both if stat in ("D", "T_LY") else sv if stat == "T_SV" else cx
        res_g = []
        for n in names:
            study = "A272" if n in P3.CAL.A272_ATLASES else "A274"
            pb = aux[(study, n)]["pb"]
            lab = np.array([cls.get((M.ALIAS.get(n, n), str(x)), "missing") for x in pb.donor])
            try:
                res_g.append(L.gate_atlas(pb[lab == arm].reset_index(drop=True), C.trio(), 50 if chk else 1000, n))
            except Exception as e:  # noqa: BLE001
                res_g.append((False, np.nan, np.nan, type(e).__name__))
        gates[rid] = res_g
        not chk and log("gates", rid, [(n, g[0]) for n, g in zip(names, res_g)])
    dwf = dict(dl)
    v283 = pd.read_csv(V283_ROW).set_index("id").loc["A283_V283"] if V283_ROW.exists() else None
    v283_ok = (v283 is not None and v283.core_status == "holds") or chk

    rows = []
    for rid, (arm, stat, nm) in ROWS.items():
        v = obs[arm][stat]
        v14, vc = v[I14], v[IC]
        dv = {lab: float(v[lay["drops"]][j]) for j, (lab, _) in enumerate(dl)}
        assert np.isclose(dv[tl], vc[0], atol=1e-9)
        pr = C.drop_percentiles(dv, tl, dwf, +1)
        gp = gates[rid]
        gate_all = bool(len(gp) and all(g[0] for g in gp))
        lo, hi = L.ci90(diffs[arm][stat])
        out = {}
        for f, sl in FAMILIES.items():
            o14, oc = float(PI(v14, sl)[0]), float(PI(vc, sl)[0])
            out[f] = (o14, pcal(null14, arm, stat, f, o14), oc, pcal(nullc, arm, stat, f, oc))
        met = [out[f][3] < 0.05 for f in FAMILIES]
        st = L.status(met) if v283_ok else "uninformative (V283 fails on the core projector or is not available)"
        rows.append(dict(
            id=rid, figure_panel="Fig 2 / S (living vs post-mortem)", statistic=nm + ", calibrated p", direction=1,
            value_14_stored=float(ref[(ref.arm == arm) & (ref.statistic == stat)].iloc[0]["clock"]), value_14=float(v14[0]),
            reproduced=True, value_core=float(vc[0]),
            core_criterion="descriptive; reference reading: p_cal < 0.05 against both families (core, 5000 redrawn placebo draws), given V283 passes",
            core_status=st, diff_14_minus_core=float(v14[0] - vc[0]), diff_lo90=lo, diff_hi90=hi,
            diff_draws="200 pair-level donor bootstraps within stratum, rng [336,1], shared (clock statistic)",
            **{k: pr[k] for k in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
            trio_gate_passes=gate_all, trio_alone_value=float(v[lay["trio"]]) if gate_all else np.nan,
            reading=C.reading(st, pr["trio_drop_below_q10"]) if v283_ok else st, blind="core blind",
            notes=(f"value = clock statistic (the p_cal are here): matched p_cal 14 {out['matched'][1]:.4f} / core {out['matched'][3]:.4f} "
                   f"(pi {out['matched'][2]:.3f}); responsive p_cal 14 {out['responsive'][1]:.4f} / core {out['responsive'][3]:.4f} "
                   f"(pi {out['responsive'][2]:.3f}); V283 core {'passes' if v283_ok else 'fails or missing'}; "
                   f"trio gate passes in {sum(g[0] for g in gp)} of {len(gp)} atlases"), script=SCRIPT))
    if chk:
        assert len(rows) == 3
        print("check ok")
        return
    df = C.write_rows("a282", rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes"]).round(4).T.to_string(), flush=True)
    log("done")


if __name__ == "__main__":
    main()
