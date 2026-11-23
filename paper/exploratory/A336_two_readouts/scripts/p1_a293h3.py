"""A336 part 1 wrapper: A293 H293.3 (Part 2: endothelial minus fibroblast untimed offset across human atlases; circular
mean over atlases) on the 14 genes, the core loop, the trio contribution, the gene-drop null and the trio alone.

Original: A293_organisation_levels/scripts/part2_phase_order.py (imported, not edited: `circ_mean`, `wrap`, `mouse_offset`).
The A272 / A274 pair tables are recomputed here from the stacked states of p1_a283.py (14 genes: compared with the stored
`pairs.csv` of A272 and A274 to 1e-6; core: clock core and the matched companions restricted to the core positions give
|r|, the offset and the matched q95; readable = |r| > matched q95). The MOUSE offsets (A273 colon, A279 aorta) are NOT
rerun on the core: the stored 14-gene untimed mouse offsets give the reference sign (deviation, mouse core reruns skipped).
Registered reading used here (user's brief): equivalence, the 90 % atlas-bootstrap interval of the human circular mean
inside +-1 h (A335 part 3); the original sign test is reported too. Drops: the offset of the drop clock on the pairs
readable on the core (no companions are drawn for drops; the trio drop equals the core). Trio contribution: 200 atlas
bootstraps (rng [336, 1]) of the circular mean, the same resamples for both waveforms. Gate: trio gate (1000
permutations) on every atlas readable on the core. NOT BLIND (the 14-gene values).

    python p1_a293h3.py [--check]
Output: ../output/part1/a293h3.csv, a293h3_core_pairs.csv (class pair rows only)
"""
import argparse
import importlib.util
import sys

import ct_lib as L
from ct_lib import *  # noqa: F401,F403
import a336_common as C
import p1_a283 as P3
from clocksync.waveform import HOURS
from scipy.stats import binomtest

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a293h3.py"
spec = importlib.util.spec_from_file_location("p1_a293_p2", L.EXP / "A293_organisation_levels/scripts/part2_phase_order.py")
P2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P2)
A272, CAL = P3.A272, P3.CAL
C1, C2 = "endothelial", "fibroblast_stromal"
W = P2.W


def pair_table(data, rows):
    """Pairs of every atlas with the complex clock r of state row `rows[0]` and the matched q95 over `rows[1]`."""
    out = []
    for study, atl in data.items():
        for n, (key, S, meta) in atl.items():
            pairs = A272.pair_index(key)
            for k, pr in enumerate(pairs):
                r = A272.pair_r(S[rows[0]], None, pr)
                q = float(np.nanquantile(CAL.pair_abs(S[rows[1]], [pr])[:, 0], 0.95))
                out.append(dict(study=study, atlas=n, a=pr["a"], b=pr["b"], family_a=meta.family_a[k], family_b=meta.family_b[k],
                                r_abs=float(abs(r)), offset_h=float(np.angle(r) * HOURS), matched_q95=q, r=r))
    return pd.DataFrame(out)


def oriented(pt):
    sel = pt[((pt.family_a == C1) & (pt.family_b == C2)) | ((pt.family_a == C2) & (pt.family_b == C1))].copy()
    sel["oriented_h"] = np.where(sel.family_a == C1, sel.offset_h, -sel.offset_h)
    return sel


def per_atlas(sel, readable):
    """Weighted circular mean per atlas over its readable pairs (weights |r|^2)."""
    s = sel[readable(sel)]
    return {n: P2.circ_mean(g.oriented_h, g.r_abs.to_numpy() ** 2) for n, g in s.groupby("atlas")}


def summary(per, sign, rng):
    v = np.array(list(per.values()))
    k = int((np.sign(v) == sign).sum())
    boot = [P2.circ_mean(v[i], np.ones(len(v))) for i in (rng.integers(0, len(v), len(v)) for _ in range(2000))]
    return dict(n=len(v), k=k, p_sign=float(binomtest(k, len(v), 0.5, alternative="greater").pvalue),
                mean=P2.circ_mean(v, np.ones(len(v))), lo=float(np.quantile(P2.wrap(boot), 0.05)), hi=float(np.quantile(P2.wrap(boot), 0.95)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(1)
    data, aux, groups = P3.load14(2 if chk else None)
    # ---- 14 genes, reproduced first ------------------------------------------------------------------------------------
    pt14 = pair_table(data, (0, MATCHED))
    stored = pd.concat([pd.read_csv(L.PAPER / "exploratory" / f / "output/pairs.csv").assign(study=s) for f, s in
                        (("A272_backbone_pairwise", "A272"), ("A274_celltype_map", "A274"))])
    mg = pt14.merge(stored, on=["atlas", "a", "b"], suffixes=("", "_s"))
    mc, ms = P2.mouse_offset("colon", C1, C2), P2.mouse_offset("aorta", C1, C2)
    sign = np.sign(mc[0])
    assert sign == np.sign(ms[0]), "mouse signs disagree"
    sel14 = oriented(pt14)
    per14 = per_atlas(sel14, lambda s: s.r_abs > s.matched_q95)
    rng = np.random.default_rng(293)
    s14 = summary(per14, sign, rng)
    ref = pd.read_csv(L.EXP / "A293_organisation_levels/output/part2/tests.csv")
    ref = ref[ref.test == "H293.3"].iloc[0]
    bad = []
    if not chk:
        d_off = np.abs(P2.wrap(mg.offset_h - mg.offset_h_s))
        for nm, g, w in (("pairs matched", len(mg), len(stored[stored.atlas.isin(pt14.atlas.unique())])),
                         ("max |d r_abs|", float(np.abs(mg.r_abs - mg.r_abs_s).max()), 0.0),
                         ("max |d offset|", float(d_off.max()), 0.0), ("max |d matched_q95|", float(np.abs(mg.matched_q95 - mg.matched_q95_s).max()), 0.0),
                         ("atlases", s14["n"], ref.atlases), ("atlases_with_mouse_sign", s14["k"], ref.atlases_with_mouse_sign),
                         ("p_sign", s14["p_sign"], ref.p_sign_one_sided), ("circular_mean", s14["mean"], ref.human_circular_mean_h),
                         ("boot_lo90", s14["lo"], ref.boot_lo90), ("boot_hi90", s14["hi"], ref.boot_hi90)):
            if not np.isclose(g, w, atol=1e-6, rtol=0):
                bad.append((nm, g, w))
        log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        C.write_rows("a293h3", [dict(id="A293_H293.3", figure_panel="Fig 4/5", statistic="endothelial - fibroblast untimed offset, circular mean over human atlases",
                                     value_14_stored=float(ref.human_circular_mean_h), value_14=s14["mean"], reproduced=False,
                                     notes=f"14-gene values not reproduced (tolerance 1e-6): {bad}", script=SCRIPT)])
        return

    # ---- core, drops, trio -------------------------------------------------------------------------------------------
    dl, tl = C.drops(n=100)
    if chk:
        dl = dl[:3] + [x for x in dl if x[0] == tl]
    lay = L.layout(len(dl))
    for study, atl in data.items():
        for n, (key, S, meta) in list(atl.items()):
            atl[n] = (key, np.vstack([S, L.stack_states(aux[(study, n)]["pb"], aux[(study, n)]["wfs"], dl, P3.STATEFN)]), meta)
    ptc = pair_table(data, (IC.start, slice(IC.start + 1, IC.start + 1 + N_COMP)))
    selc = oriented(ptc)
    readable_c = lambda s: s.r_abs > s.matched_q95  # noqa: E731
    perc = per_atlas(selc, readable_c)
    sc = summary(perc, sign, np.random.default_rng(293))
    if not chk:
        selc.drop(columns=["r"]).to_csv(L.EXP / "A336_two_readouts/output/part1/a293h3_core_pairs.csv", index=False)

    # trio contribution: 200 atlas bootstraps, shared resamples
    atl_all = sorted(set(sel14.atlas) | set(selc.atlas))
    brng = np.random.default_rng([336, 1])
    diffs = []
    for _ in range(3 if chk else 200):
        pick = list(brng.choice(atl_all, len(atl_all)))
        m = []
        for per in (per14, perc):
            v = np.array([per[x] for x in pick if x in per])
            m.append(P2.circ_mean(v, np.ones(len(v))) if len(v) else np.nan)
        diffs.append(float(P2.wrap(m[0] - m[1])))

    # drops: offset of the drop clock on the pairs readable on the core, per atlas, circular mean over atlases
    ready = selc[readable_c(selc)][["atlas", "a", "b"]]
    drop_val, trio_val = {}, None

    def stat_rows(row):
        out = []
        for study, atl in data.items():
            for n, (key, S, meta) in atl.items():
                g = ready[ready.atlas == n]
                if not len(g):
                    continue
                pairs = A272.pair_index(key)
                offs, ws = [], []
                for pr in pairs:
                    if ((g.a == pr["a"]) & (g.b == pr["b"])).any():
                        r = A272.pair_r(S[row], None, pr)
                        fa = meta[(meta.a == pr["a"]) & (meta.b == pr["b"])].family_a.iloc[0]
                        h = np.angle(r) * HOURS * (1 if fa == C1 else -1)
                        offs.append(h)
                        ws.append(abs(r) ** 2)
                out.append(P2.circ_mean(offs, np.array(ws)))
        return P2.circ_mean(out, np.ones(len(out))) if out else np.nan

    for j, (lab, _) in enumerate(dl):
        drop_val[lab] = stat_rows(lay["drops"].start + j)
    trio_val = stat_rows(lay["trio"])
    assert np.isclose(P2.wrap(drop_val[tl] - sc["mean"]), 0, atol=1e-6), (drop_val[tl], sc["mean"])
    pr_ = C.drop_percentiles(drop_val, tl, dict(dl), +1)
    gate_atl = sorted(perc)
    gates = {}
    for n in gate_atl:
        study = "A272" if n in CAL.A272_ATLASES else "A274"
        gates[n] = L.gate_atlas(aux[(study, n)]["pb"], C.trio(), 50 if chk else 1000, n)
    gate_all = bool(len(gates) and all(g[0] for g in gates.values()))
    inside = (sc["lo"] >= -1) and (sc["hi"] <= 1)
    inside14 = (s14["lo"] >= -1) and (s14["hi"] <= 1)
    st = "holds" if inside else "falls"
    row = dict(
        id="A293_H293.3", figure_panel="Fig 4/5", statistic="endothelial minus fibroblast untimed offset (h), circular mean over human atlases",
        direction=1, value_14_stored=float(ref.human_circular_mean_h), value_14=s14["mean"], reproduced=True, value_core=sc["mean"],
        core_criterion="equivalence: 90 % atlas-bootstrap interval (2000, rng 293) of the core circular mean inside +-1 h (A335 part 3); sign test reported",
        core_status=st, diff_14_minus_core=float(P2.wrap(s14["mean"] - sc["mean"])), diff_lo90=L.ci90(diffs)[0], diff_hi90=L.ci90(diffs)[1],
        diff_draws="200 atlas bootstraps (rng [336,1]), shared", **{k: pr_[k] for k in (
            "n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
        trio_gate_passes=gate_all, trio_alone_value=trio_val if gate_all else np.nan,
        reading=C.reading(st, pr_["trio_drop_below_q10"]), blind="core not blind (A335 / A329: 14-gene values)",
        notes=(f"14 genes: {s14['n']} atlases, mean {s14['mean']:.4f} h, 90% [{s14['lo']:.4f}, {s14['hi']:.4f}] (inside +-1 h: {inside14}), mouse-sign {s14['k']} of {s14['n']} "
               f"p {s14['p_sign']:.4f}; core: {sc['n']} readable atlases, mean {sc['mean']:.4f} h, 90% [{sc['lo']:.4f}, {sc['hi']:.4f}], mouse-sign {sc['k']} of {sc['n']} "
               f"(one-sided sign p {sc['p_sign']:.4f}); mouse reference sign from the stored 14-gene mouse offsets (colon {mc[0]:.3f} h, aorta {ms[0]:.3f} h; A273/A279 not rerun on the core); "
               f"gate per atlas {{{', '.join(f'{n}: {g[0]}' for n, g in gates.items())}}}; drops use the core-readable pairs (no companions for drops)"),
        script=SCRIPT)
    if chk:
        print("check ok")
        return
    df = C.write_rows("a293h3", [row])
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script"]).T.to_string(), flush=True)
    log("done")


if __name__ == "__main__":
    main()
