"""A336 part 1 wrapper: A280 H280.1 (D = T_SV - T_LY, IBDverse rectum + ileum) and H280.2 (T_CX, PsychAD MSSM) in the A283
form on the 14 genes, the core loop, the trio contribution, the gene-drop null and the trio alone.

Original: A280_celltype_replication2/scripts/calibrated_a280.py (imported, not edited; its `draw` is reused, seeds
[280, 9, d]). Same stacked-state design as p1_a283.py (ct_lib.py). The placebo null (5000 draws) is redrawn for the core
projector; the 14-gene p_cal comes from the stored null.npz after 10 recomputed draws are shown identical. Core p_cal are
read only as informative if V283 passed on the core projector (output/part1/a283.csv, row A283_V283). Gene drops: 100
random + trio (clock only). Trio contribution: 200 pair-level donor bootstraps within stratum, rng [336, 1].

    python p1_a280.py [--check] [--draws 5000] [--jobs 4]
Output: ../output/part1/a280.csv, a280_core_null_pi.npz
"""
import argparse
import importlib.util
import sys

import ct_lib as L
from ct_lib import *  # noqa: F401,F403
import a336_common as C

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a280.py"
spec = importlib.util.spec_from_file_location("p1_a280_orig", L.PAPER / "exploratory/A280_celltype_replication2/scripts/calibrated_a280.py")
M = importlib.util.module_from_spec(spec)
sys.modules["p1_a280_orig"] = M
spec.loader.exec_module(M)
A150, A272, A274, C276 = M.A150, M.A272, M.A274, M.C276
META = A272.META
STORED = L.PAPER / "exploratory/A280_celltype_replication2/output/calibrated"
V283_ROW = L.EXP / "A336_two_readouts/output/part1/a283.csv"
TESTS = {"H280.1 D": ("A280_H280.1", "D = T_SV - T_LY (IBDverse rectum + ileum)"), "H280.2 T_CX": ("A280_H280.2", "T_CX (PsychAD MSSM)")}
STATEFN = lambda pb, w: A272.atlas_states(pb[META + list(w.genes)], w)[0]  # noqa: E731


def load14(subset):
    A150.ATLASES.update(A274.ALL)
    data, aux = {}, {}
    for n in list(A274.ALL)[:subset]:
        try:
            A = A150.load_atlas(n, None, 0)
        except Exception as e:  # noqa: BLE001
            log(n, "left out:", type(e).__name__)
            continue
        if A["rep"]["rows_with_count_diff"] != 0:
            log(n, "left out: twin rule")
            continue
        pb = A["pb"].reset_index(drop=True)
        A["pb"] = pb
        wfs = [cs.gtex()]
        for kind in A272.KINDS:
            wfs += A274.pseudo_sets(A, list(A274.ALL).index(n), N_COMP, kind)
        S = np.array([STATEFN(pb, w) for w in wfs])
        key = pb[["donor", "compartment", "stratum"]].copy()
        pairs = A272.pair_index(key)
        meta = pd.DataFrame([dict(atlas=n, a=p["a"], b=p["b"], family_a=A150.FAM.family(p["a"]),
                                  family_b=A150.FAM.family(p["b"])) for p in pairs])
        cortex = n in A274.CORTEX
        meta["klass"] = [("CX_diff" if fa != fb else "CX_same") if cortex else A272.klass(A272.group(fa), A272.group(fb))
                         for fa, fb in zip(meta.family_a, meta.family_b)]
        data[n] = (key, S, meta)
        cols = META + sorted({g for w in wfs for g in w.genes} | set(cs.gtex().genes))
        aux[n] = dict(pb=pb[cols].copy(), wfs=wfs)
        log(n, "rows", len(key), "pairs", len(pairs))
    cnt = pd.concat([m for _, _, m in data.values()]).groupby(["atlas", "klass"]).size().unstack(fill_value=0)
    col = lambda c: cnt[c] if c in cnt else pd.Series(0, index=cnt.index)  # noqa: E731
    body = [n for n in data if n in A274.BODY and col("SV_X")[n] >= A274.MIN_PAIRS and col("LY_X")[n] >= A274.MIN_PAIRS]
    cx = [n for n in data if n in A274.CORTEX and col("CX_diff")[n] >= A274.MIN_PAIRS]
    return data, aux, dict(body=body, cx=cx)


def boot_draw(data, groups, rng):
    rows, metas = [], []
    for n, (key, S, meta) in data.items():
        rows.append(L.boot_abs(S, A272.pair_index(key), rng))
        metas.append(meta)
    return M.statistics(C276.class_medians(np.concatenate(rows, 1), pd.concat(metas, ignore_index=True)), groups)


def null_task(d):
    return M.draw([A274.SEED, 9, d], G["data"], G["groups"])[0]


def pi_(v, sl):
    return float(M.pi(v, sl)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    data, aux, groups = load14(2 if chk else None)
    not chk and log("groups", groups)
    ref = pd.read_csv(STORED / "tests.csv").set_index("test")

    # ---- 14 genes, reproduced first -------------------------------------------------------------------------------
    obs14 = M.draw(None, data, groups)[0]
    z = np.load(STORED / "null.npz")
    for d in range(10):
        got = M.draw([A274.SEED, 9, d], data, groups)[0]
        for t in z.files:
            assert chk or np.allclose(got[t], z[t][d], atol=1e-9, equal_nan=True), (t, d, "stored null differs")
    pis14 = {t: {k: M.pi(z[t], sl) for k, sl in FAMILIES.items()} for t in z.files}
    bad = []
    if not chk:
        for t in TESTS:
            if t not in obs14:
                bad.append((t, "absent"))
                continue
            row = ref.loc[t]
            checks = [("clock", obs14[t][0], row["clock"])]
            for k, sl in FAMILIES.items():
                o = pi_(obs14[t], sl)
                checks += [(f"{k}_pi", o, row[f"{k}_pi"]), (f"{k}_median", float(np.nanmedian(obs14[t][sl])), row[f"{k}_median"]),
                           (f"{k}_p_cal", float((1 + np.sum(pis14[t][k] >= o)) / (1 + len(pis14[t][k]))), row[f"{k}_p_cal"])]
            bad += [(t, nm, g, w) for nm, g, w in checks if not np.isclose(g, w, atol=1e-6, rtol=0)]
        log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        C.write_rows("a280", [dict(id=TESTS[t][0], figure_panel="Fig 2", statistic=TESTS[t][1], value_14=obs14.get(t, [np.nan])[0],
                                   reproduced=False, notes=f"14-gene values not reproduced (tolerance 1e-6): {bad}", script=SCRIPT)
                              for t in TESTS])
        return

    # ---- core, drops, trio states -----------------------------------------------------------------------------------
    dl, tl = C.drops(n=100)
    if chk:
        dl = dl[:3] + [x for x in dl if x[0] == tl]
    lay = L.layout(len(dl))
    for n, (key, S, meta) in list(data.items()):
        data[n] = (key, np.vstack([S, L.stack_states(aux[n]["pb"], aux[n]["wfs"], dl, STATEFN)]), meta)
    allv = M.draw(None, data, groups)[0]
    G["groups"] = groups
    G["data"] = {n: (k, S[IC], m) for n, (k, S, m) in data.items()}
    ndraw = 20 if chk else a.draws
    t0 = time.time()
    res = L.run_pool(null_task, range(ndraw), a.jobs)
    log(f"core null: {ndraw} draws in {(time.time() - t0) / 60:.1f} min")
    pisc = {t: {k: M.pi(np.array([r[t] for r in res]), sl) for k, sl in FAMILIES.items()} for t in res[0] if t in TESTS}
    if not chk:
        np.savez_compressed(L.EXP / "A336_two_readouts/output/part1/a280_core_null_pi.npz",
                            **{f"{t}|{k}": v for t in pisc for k, v in pisc[t].items()})
    v283 = pd.read_csv(V283_ROW).set_index("id").loc["A283_V283"] if V283_ROW.exists() else None
    v283_ok = (v283 is not None and v283.core_status == "holds") or chk

    # ---- trio contribution and gates ---------------------------------------------------------------------------------
    brng = np.random.default_rng([336, 1])
    D2 = {n: (k, S[[0, IC.start]], m) for n, (k, S, m) in data.items()}
    diffs = {t: [] for t in TESTS}
    for b in range(3 if chk else 200):
        v = boot_draw(D2, groups, brng)
        for t in diffs:
            if t in v:
                diffs[t].append(v[t][0] - v[t][1])
    used = {"H280.1 D": groups["body"], "H280.2 T_CX": groups["cx"]}
    gates = {n: L.gate_atlas(aux[n]["pb"], C.trio(), 50 if chk else 1000, n) for n in sorted({x for v in used.values() for x in v})}
    not chk and log("gates", {n: (v[0], round(v[1], 3), round(v[2], 3)) for n, v in gates.items()})
    dwf = dict(dl)

    rows = []
    for t, (rid, nm) in TESTS.items():
        if t not in allv:
            continue
        v = allv[t]
        v14, vc = v[I14], v[IC]
        dv = {lab: float(v[lay["drops"]][j]) for j, (lab, _) in enumerate(dl)}
        assert np.isclose(dv[tl], vc[0], atol=1e-9)
        pr = C.drop_percentiles(dv, tl, dwf, +1)
        gp = [gates[n] for n in used[t]]
        gate_all = bool(all(g[0] for g in gp))
        lo, hi = L.ci90(diffs[t])
        out = {}
        for k, sl in FAMILIES.items():
            o14, oc = pi_(v14, sl), pi_(vc, sl)
            out[k] = (o14, float((1 + np.sum(pis14[t][k] >= o14)) / (1 + len(pis14[t][k]))), oc,
                      float((1 + np.sum(pisc[t][k] >= oc)) / (1 + len(pisc[t][k]))))
        met = [out[k][3] < 0.05 for k in FAMILIES]
        st = L.status(met) if v283_ok else "uninformative (V283 fails on the core projector or is not available)"
        rows.append(dict(
            id=rid, figure_panel="Fig 2 (replication)", statistic=nm + ", calibrated p", direction=1,
            value_14_stored=float(ref.loc[t, "clock"]), value_14=float(v14[0]), reproduced=True, value_core=float(vc[0]),
            core_criterion="p_cal < 0.05 against both families (core, 5000 redrawn placebo draws), given V283 passes on the core projector",
            core_status=st, diff_14_minus_core=float(v14[0] - vc[0]), diff_lo90=lo, diff_hi90=hi,
            diff_draws="200 pair-level donor bootstraps within stratum, rng [336,1], shared (clock statistic)",
            **{k: pr[k] for k in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
            trio_gate_passes=gate_all, trio_alone_value=float(v[lay["trio"]]) if gate_all else np.nan,
            reading=C.reading(st, pr["trio_drop_below_q10"]) if v283_ok else st, blind="core blind",
            notes=(f"value = clock statistic (the p_cal are here): matched p_cal 14 {out['matched'][1]:.4f} (stored {ref.loc[t, 'matched_p_cal']:.4f}) / core "
                   f"{out['matched'][3]:.4f} (pi {out['matched'][2]:.3f}); responsive p_cal 14 {out['responsive'][1]:.4f} (stored "
                   f"{ref.loc[t, 'responsive_p_cal']:.4f}) / core {out['responsive'][3]:.4f} (pi {out['responsive'][2]:.3f}); V283 core read from a283.csv: "
                   f"{'passes' if v283_ok else 'fails or missing'}; 14-gene p_cal from the stored null (10 draws recomputed identical); "
                   f"trio gate passes in {sum(g[0] for g in gp)} of {len(gp)} atlases; atlases {used[t]}"), script=SCRIPT))
    if chk:
        assert len(rows) == len(TESTS) or len(rows) >= 1
        print("check ok")
        return
    df = C.write_rows("a280", rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes"]).round(4).T.to_string(), flush=True)
    log("done")


if __name__ == "__main__":
    main()
