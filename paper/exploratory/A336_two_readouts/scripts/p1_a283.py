"""A336 part 1 wrapper: A272 H272.1-2, A274 H274.1-3, A283 H283.1-5 and V283 (Figure 2 cell-type ranking, anchor-free
pairwise readout) on the 14 genes, the core loop, the trio contribution, the gene-drop null and the trio alone.

Originals (imported, not edited): A276_pairwise_calibration/scripts/{calibrate_a276,corrected_a277,calibrated_a278,
calibrated_a283}.py with A272 / A274 (`score_a272`, `score_a274`). One stacked state matrix per atlas (ct_lib.py): 14-gene
block, core block (companions drawn as originally with the 14-gene waveform, restricted to the core positions), the 101
clock-only drop rows (drops(n=100)) and the trio. The A283 placebo null (5000 draws, seeds [283, d]) is REDRAWN for the
core projector (the draw only permutes donors; the same permutation serves every row). The 14-gene p_cal is read from the
stored null after checking that 10 recomputed draws equal the stored ones to 1e-9 (the null of the 14-gene projector is
deterministic, so a redraw would be identical; deviation from "redraw per projector" for the 14-gene block only).
Trio contribution: 200 pair-level donor bootstraps within stratum (A272's offset bootstrap), rng [336, 1], the same
resamples for both waveforms. Gate: model-free gate (1000 permutations) of the trio on every atlas entering the statistic;
trio_gate_passes is True only if all pass.

    python p1_a283.py [--check] [--draws 5000] [--jobs 4]
Outputs: ../output/part1/{a272,a274,a283}.csv, a283_core_null_pi.npz
"""
import argparse
import sys
import time

import ct_lib as L
from ct_lib import *  # noqa: F401,F403
import a336_common as C

sys.path.insert(0, str(L.PAPER / "exploratory/A276_pairwise_calibration/scripts"))
import calibrate_a276 as CAL  # noqa: E402
from calibrated_a278 import p_cal, pi  # noqa: E402

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a283.py"
A272, A274, A150, META = CAL.A272, CAL.A274, CAL.A150, CAL.A272.META
STORED = L.PAPER / "exploratory/A276_pairwise_calibration/output/a283"
TESTS = ["H272.1 T_SV", "H272.2 D", "H274.1 D", "H274.2 T_SV", "H274.3 T_CX"]
H283 = {"H283.1": ("H272.1 T_SV", "both"), "H283.2": ("H272.2 D", "both"), "H283.3": ("H274.1 D", "both"),
        "H283.4": ("H274.2 T_SV", "matched"), "H283.5": ("H274.3 T_CX", "both")}
ROW = {"H272.1 T_SV": ("A272_H272.1", "a272", "A272"), "H272.2 D": ("A272_H272.2", "a272", "A272"),
       "H274.1 D": ("A274_H274.1", "a274", "A274"), "H274.2 T_SV": ("A274_H274.2", "a274", "A274"),
       "H274.3 T_CX": ("A274_H274.3", "a274", "A274")}
NAMES = {"H272.1 T_SV": "T_SV (5 atlases)", "H272.2 D": "D = T_SV - T_LY (A272 atlases with lymphocytes)",
         "H274.1 D": "D = T_SV - T_LY (A274 body atlases)", "H274.2 T_SV": "T_SV (A274 body atlases)",
         "H274.3 T_CX": "T_CX (cortex cell classes)"}
STATEFN = lambda pb, w: A272.atlas_states(pb[META + list(w.genes)], w)[0]  # noqa: E731


def load14(subset):
    """Atlases, 14-gene waveforms and states (as calibrate_a276.load), pairs, classes and groups (as corrected_a277)."""
    A150.ATLASES.update(A274.ALL)
    data, aux, groups = {}, {}, {}
    for study, names in (("A272", CAL.A272_ATLASES), ("A274", CAL.A274_ATLASES)):
        data[study] = {}
        for n in names[:subset]:
            A = A150.load_atlas(n, None, 0)
            assert study == "A272" or A["rep"]["rows_with_count_diff"] == 0, n
            A["pb"] = A["pb"].reset_index(drop=True)
            wfs = [cs.gtex()]
            for kind in CAL.KINDS:
                wfs += (A272.pseudo_sets(A, CAL.A272_ATLASES.index(n), N_COMP, kind) if study == "A272"
                        else A274.pseudo_sets(A, list(A274.ALL).index(n), N_COMP, kind))
            S = np.array([STATEFN(A["pb"], w) for w in wfs])
            key = A["pb"][["donor", "compartment", "stratum"]].copy()
            pairs = A272.pair_index(key)
            meta = pd.DataFrame([dict(atlas=n, a=p["a"], b=p["b"], family_a=A150.FAM.family(p["a"]),
                                      family_b=A150.FAM.family(p["b"])) for p in pairs])
            meta["klass"] = [CAL.klass(study, n, fa, fb) for fa, fb in zip(meta.family_a, meta.family_b)]
            data[study][n] = (key, S, meta)
            cols = META + sorted({g for w in wfs for g in w.genes} | set(cs.gtex().genes))
            aux[(study, n)] = dict(pb=A["pb"][cols].copy(), wfs=wfs)
            log(study, n, "rows", len(key), "pairs", len(pairs))
            del A
        cnt = pd.concat([m for _, _, m in data[study].values()]).groupby(["atlas", "klass"]).size().unstack(fill_value=0)
        col = lambda c: cnt[c] if c in cnt else pd.Series(0, index=cnt.index)  # noqa: E731
        nm = list(data[study])
        if study == "A272":
            groups[study] = dict(sv=nm, ly=[n for n in nm if col("LY_X")[n] > 0])
        else:
            body = [n for n in nm if n in A274.BODY]
            sv = [n for n in body if col("SV_X")[n] >= A274.MIN_PAIRS]
            groups[study] = dict(sv=sv, both=[n for n in sv if col("LY_X")[n] >= A274.MIN_PAIRS],
                                 cx=[n for n in nm if n in A274.CORTEX and col("CX_diff")[n] >= A274.MIN_PAIRS])
    return data, aux, groups


def stats_of(cm, groups):
    full = {}
    for t in TESTS:
        study, g = ("A272", "sv" if t == "H272.1 T_SV" else "ly") if t.startswith("H272") else (
            "A274", {"H274.1 D": "both", "H274.2 T_SV": "sv", "H274.3 T_CX": "cx"}[t])
        if groups[study][g]:
            full[t] = CAL.statistics(cm, groups)[t]
    return full


def draw(seed, data, groups, boot=None):
    """{test: vector over the sets of the stacked matrix} on the real pairing (seed None) or a placebo draw (as
    corrected_a277.draw; same rng order)."""
    warnings.filterwarnings("ignore")
    rng = None if seed is None else np.random.default_rng(seed)
    cm = {}
    for study, atl in data.items():
        rows, metas = [], []
        for n, (key, S, meta) in atl.items():
            k = key if rng is None else A272.placebo(key, rng)
            pairs = A272.pair_index(k)
            if boot is None:
                rows.append(CAL.pair_abs(S, pairs))
            else:
                rows.append(L.boot_abs(S, pairs, boot))
            metas.append(meta)
        cm[study] = CAL.class_medians(np.concatenate(rows, 1), pd.concat(metas, ignore_index=True))
    return stats_of(cm, groups)


def null_task(d):
    return draw([283, d], G["data"], G["groups"])


def close(a, b):
    return bool(np.isclose(a, b, atol=1e-6, rtol=0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    subset = 2 if chk else None
    data, aux, groups = load14(subset)
    not chk and log("groups", groups)
    ent = {"a272": [], "a274": [], "a283": []}

    # ---- 14 genes, reproduced first ---------------------------------------------------------------------------------
    obs14 = draw(None, data, groups)
    z = np.load(STORED / "null_draws.npz")
    for d in range(10):
        got = draw([283, d], data, groups)
        for t in got:
            assert np.allclose(got[t], z[t][d], atol=1e-9, equal_nan=True) or chk, (t, d, "stored null differs")
    ref = {"A272": pd.read_csv(L.PAPER / "exploratory/A272_backbone_pairwise/output/tests.csv"),
           "A274": pd.read_csv(L.PAPER / "exploratory/A274_celltype_map/output/tests.csv")}
    ref3 = pd.read_csv(STORED / "tests.csv").set_index("test")
    cal_ref = pd.read_csv(STORED / "calibration.csv").set_index(["test", "family"])
    stored, bad = {}, []
    pis14 = {t: {k: pi(z[t], sl) for k, sl in FAMILIES.items()} for t in z.files}
    cal14 = {}
    for t in z.files:
        for k in FAMILIES:
            p = pis14[t][k]
            h = len(p) // 2
            cal14[(t, k)] = (float(np.mean(p_cal(p[h:], p[:h]) < 0.05)), float(np.mean(p_cal(p[:h], p[h:]) < 0.05)))
    if not chk:
        for t in TESTS:
            study = ROW[t][2]
            r = ref[study].set_index(["test", "family"])
            for k, sl in FAMILIES.items():
                q95 = float(np.nanquantile(obs14[t][sl], 0.95))
                for nm_, got, want in (("clock", obs14[t][0], r.loc[(t, k), "clock"]),
                                       ("pseudo_q95", q95, r.loc[(t, k), "pseudo_q95"])):
                    if not close(got, want):
                        bad.append((t, k, nm_, got, want))
            h3 = [h for h, (tt, _) in H283.items() if tt == t][0]
            for k, sl in FAMILIES.items():
                o = pi(obs14[t], sl)[0]
                pc = float(p_cal(o, pis14[t][k])[0])
                for nm_, got, want in ((f"{k}_pi", o, ref3.loc[h3, f"{k}_pi"]), (f"{k}_p_cal", pc, ref3.loc[h3, f"{k}_p_cal"]),
                                       (f"{k}_median", float(np.nanmedian(obs14[t][sl])), ref3.loc[h3, f"{k}_median"])):
                    if not close(got, want):
                        bad.append((h3, nm_, got, want))
        for (t, k), (r1, r2) in cal14.items():
            if not (close(r1, cal_ref.loc[(t, k), "rate_null_first"]) and close(r2, cal_ref.loc[(t, k), "rate_null_second"])):
                bad.append((t, k, "V283 rates", r1, r2))
        log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        for e in ent:
            C.write_rows(e, [dict(id=ROW[t][0], figure_panel="Fig 2a", statistic=NAMES[t], value_14=obs14[t][0], reproduced=False,
                                  notes=f"14-gene values not reproduced (tolerance 1e-6): {bad}", script=SCRIPT)
                             for t in TESTS if ROW[t][1] == e] or [dict(id="A283_V283", reproduced=False, notes=f"not reproduced: {bad}", script=SCRIPT)])
        return

    # ---- core, drops, trio states (after the reproduction check) ----------------------------------------------------
    dl, tl = C.drops(n=100)
    if chk:
        dl = dl[:3] + [x for x in dl if x[0] == tl]
    lay = L.layout(len(dl))
    for study, atl in data.items():
        for n, (key, S, meta) in list(atl.items()):
            ex = L.stack_states(aux[(study, n)]["pb"], aux[(study, n)]["wfs"], dl, STATEFN)
            atl[n] = (key, np.vstack([S, ex]), meta)
    not chk and log("stacked states", {s: {n: v[1].shape for n, v in atl.items()} for s, atl in data.items()})
    allv = draw(None, data, groups)

    # ---- core placebo null (redrawn) ---------------------------------------------------------------------------------
    G["groups"] = groups
    G["data"] = {s: {n: (k, S[IC], m) for n, (k, S, m) in atl.items()} for s, atl in data.items()}
    ndraw = 20 if chk else a.draws
    t0 = time.time()
    res = L.run_pool(null_task, range(ndraw), a.jobs)
    log(f"core null: {ndraw} draws in {(time.time() - t0) / 60:.1f} min")
    nullc = {t: np.array([r[t] for r in res]) for t in res[0]}
    pisc = {t: {k: pi(nullc[t], sl) for k, sl in FAMILIES.items()} for t in nullc}
    if not chk:
        np.savez_compressed(L.EXP / "A336_two_readouts/output/part1/a283_core_null_pi.npz",
                            **{f"{t}|{k}": v for t in pisc for k, v in pisc[t].items()})
    half = ndraw // 2
    calc = {}
    for t in pisc:
        for k in FAMILIES:
            p = pisc[t][k]
            calc[(t, k)] = (float(np.mean(p_cal(p[half:], p[:half]) < 0.05)), float(np.mean(p_cal(p[:half], p[half:]) < 0.05)))
    v283_ok = all(max(r) <= 0.0685 for r in calc.values())
    not chk and log("V283 core:", "PASS" if v283_ok else "FAIL", {f"{t}|{k}": [round(x, 4) for x in r] for (t, k), r in calc.items()})

    # ---- trio contribution, gates ---------------------------------------------------------------------------------
    brng = np.random.default_rng([336, 1])
    D2 = {s: {n: (k, S[[0, IC.start]], m) for n, (k, S, m) in atl.items()} for s, atl in data.items()}
    nb = 3 if chk else 200
    diffs = {t: [] for t in allv}
    for b in range(nb):
        v = draw(None, D2, groups, boot=brng)
        for t in v:
            diffs[t].append(v[t][0] - v[t][1])
    log("bootstraps done")
    used = {"H272.1 T_SV": [("A272", n) for n in groups["A272"]["sv"]], "H272.2 D": [("A272", n) for n in groups["A272"]["ly"]],
            "H274.1 D": [("A274", n) for n in groups["A274"]["both"]], "H274.2 T_SV": [("A274", n) for n in groups["A274"]["sv"]],
            "H274.3 T_CX": [("A274", n) for n in groups["A274"]["cx"]]}
    gates = {}
    for sn in sorted({x for v in used.values() for x in v}):
        gates[sn] = L.gate_atlas(aux[sn]["pb"], C.trio(), 50 if chk else 1000, sn[1])
    not chk and log("gates", {f"{k[1]}": (v[0], round(v[1], 3), round(v[2], 3)) for k, v in gates.items()})
    dwf = dict(dl)

    # ---- rows -------------------------------------------------------------------------------------------------------
    rows = {"a272": [], "a274": [], "a283": []}
    info = {}
    for t in allv:
        v = allv[t]
        vc, v14 = v[IC], v[I14]
        q = {k: (float(np.nanquantile(vc[sl], 0.95)), float(np.nanquantile(v14[sl], 0.95))) for k, sl in FAMILIES.items()}
        dv = {lab: float(v[lay["drops"]][j]) for j, (lab, _) in enumerate(dl)}
        assert np.isclose(dv[tl], vc[0], atol=1e-9)
        pr = C.drop_percentiles(dv, tl, dwf, +1)
        gp = [gates[s] for s in used[t]]
        gate_all = bool(all(g[0] for g in gp))
        lo, hi = L.ci90(diffs[t])
        info[t] = dict(vc=vc, v14=v14, q=q, pr=pr, gate_all=gate_all, ngate=f"{sum(g[0] for g in gp)} of {len(gp)}",
                       trio=float(v[lay["trio"]]), lo=lo, hi=hi, diff=float(v14[0] - vc[0]))
    for t, (rid, entry, study) in ROW.items():
        if t not in info:
            continue
        i = info[t]
        r = ref[study].set_index(["test", "family"])
        met = [i["vc"][0] > i["q"][k][0] for k in FAMILIES]
        st = L.status(met)
        rows[entry].append(dict(
            id=rid, figure_panel="Fig 2a", statistic=NAMES[t], direction=1, value_14_stored=float(r.loc[(t, "matched"), "clock"]),
            value_14=float(i["v14"][0]), reproduced=True, value_core=float(i["vc"][0]),
            core_criterion="clock above the core q95 of both companion families (matched and responsive, core-restricted draws)",
            core_status=st, diff_14_minus_core=i["diff"], diff_lo90=i["lo"], diff_hi90=i["hi"],
            diff_draws="200 pair-level donor bootstraps within stratum, rng [336,1], shared", **{k: i["pr"][k] for k in (
                "n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
            trio_gate_passes=i["gate_all"], trio_alone_value=i["trio"] if i["gate_all"] else np.nan,
            reading=C.reading(st, i["pr"]["trio_drop_below_q10"]), blind="core blind",
            notes=(f"core q95 matched {i['q']['matched'][0]:.4f} responsive {i['q']['responsive'][0]:.4f} "
                   f"(14-gene {i['q']['matched'][1]:.4f} / {i['q']['responsive'][1]:.4f}); core clock above matched {met[0]}, responsive {met[1]}; "
                   f"trio gate passes in {i['ngate']} atlases; drops: {len(dl) - 1} random + trio"), script=SCRIPT))
    for h, (t, rule) in H283.items():
        if t not in info:
            continue
        i = info[t]
        out = {}
        for k, sl in FAMILIES.items():
            o14, oc = float(pi(i["v14"], sl)[0]), float(pi(i["vc"], sl)[0])
            out[k] = (o14, float(p_cal(o14, pis14[t][k])[0]), oc, float(p_cal(oc, pisc[t][k])[0]))
        met = [out["matched"][3] < 0.05] + ([out["responsive"][3] < 0.05] if rule == "both" else [])
        st = L.status(met) if v283_ok else "uninformative (V283 fails on the core projector)"
        rd = C.reading(st, i["pr"]["trio_drop_below_q10"]) if v283_ok else st
        rows["a283"].append(dict(
            id=f"A283_{h}", figure_panel="Fig 2a (calibrated)", statistic=f"{NAMES[t]}, calibrated p (matched decides: {rule == 'matched'})",
            direction=1, value_14_stored=float(ref3.loc[h, "clock"]), value_14=float(i["v14"][0]), reproduced=True,
            value_core=float(i["vc"][0]),
            core_criterion=("p_cal < 0.05 against " + ("the matched family (core, 5000 redrawn placebo draws; responsive reported)"
                                                       if rule == "matched" else "both families (core, 5000 redrawn placebo draws)")) + ", given V283 passes on the core projector",
            core_status=st, diff_14_minus_core=i["diff"], diff_lo90=i["lo"], diff_hi90=i["hi"],
            diff_draws="200 pair-level donor bootstraps within stratum, rng [336,1], shared (clock statistic)", **{k: i["pr"][k] for k in (
                "n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
            trio_gate_passes=i["gate_all"], trio_alone_value=i["trio"] if i["gate_all"] else np.nan, reading=rd, blind="core blind",
            notes=(f"value = clock statistic (the p_cal are here): matched p_cal 14 {out['matched'][1]:.4f} (stored {ref3.loc[h, 'matched_p_cal']:.4f}) / core "
                   f"{out['matched'][3]:.4f} (pi {out['matched'][2]:.3f}); responsive p_cal 14 {out['responsive'][1]:.4f} (stored "
                   f"{ref3.loc[h, 'responsive_p_cal']:.4f}) / core {out['responsive'][3]:.4f} (pi {out['responsive'][2]:.3f}); V283 core "
                   f"{'passes' if v283_ok else 'FAILS'}; 14-gene p_cal from the stored null (10 draws recomputed identical); trio gate passes in {i['ngate']} atlases"),
            script=SCRIPT))
    worst = max(max(r) for r in calc.values())
    rows["a283"].append(dict(
        id="A283_V283", figure_panel="Fig 2a (calibration)", statistic="max rejection rate of p_cal < 0.05 over null-only halves, ten test x family combinations",
        direction=-1, value_14_stored=float(cal_ref[["rate_null_first", "rate_null_second"]].to_numpy().max()),
        value_14=float(max(max(r) for r in cal14.values())), reproduced=True, value_core=worst,
        core_criterion="each of the ten rejection rates <= 0.0685 in both directions (core projector, 5000 redrawn placebo draws)",
        core_status="holds" if v283_ok else "falls", reading="calibration holds" if v283_ok else
        "calibration FAILS: the core cell-type p_cal (A283, A280, A282, A288) are uninformative",
        blind="core blind", trio_gate_passes="not computable",
        notes="core rates (first, second half) per test|family: " + "; ".join(f"{t}|{k} {r[0]:.4f}/{r[1]:.4f}" for (t, k), r in calc.items()),
        script=SCRIPT))
    if chk:
        assert all(len(v) for v in rows.values())
        print("check ok")
        return
    for e, rr in rows.items():
        C.write_rows(e, rr)
    pd.set_option("display.width", 250)
    for e in rows:
        print(pd.DataFrame(rows[e]).drop(columns=["script", "notes"], errors="ignore").round(4).T.to_string(), flush=True)
    log("done")


if __name__ == "__main__":
    main()
