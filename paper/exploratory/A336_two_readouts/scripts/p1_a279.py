"""A336 part 1 wrapper: A279 V279.1, H279.1 and the Bmal1-KO descriptive row (mouse aorta, Auerbach 2026) on the 14 genes, the
core loop (Per1, Per2, Nr1d1 removed; stored under human symbols), the trio contribution, the gene-drop null and the trio alone.

Original: A279_aorta_pairwise_truth/scripts/score_a279.py (imported, not edited; `tables`, `wrap`, constants). The two
functions that read the waveform (`untimed_pairs`, `timed`) are copied with a `wf` argument; the aligned primary set (both
timed couplings >= 0.5) is rebuilt for each waveform with the original rule, pair counts reported. No companions exist.
14 genes: tests.csv rows and pairs.csv reproduced to 1e-6 first. Criteria: V279.1 median |error| <= 1 h; H279.1 as the
original (diff < 0 and jackknife 90 % upper bound < 0); Bmal1 KO is descriptive: operationalised here as "holds" if every
KO timed coupling <= 0.20 and the three stromal-vascular pairs keep |r| >= 0.5, "falls" otherwise. Trio contribution: V279.1
200 pair bootstraps over the aligned all-pairs table (6 pairs), H279.1 and the KO row 200 library bootstraps (libraries
relabelled; rng [336, 1]), same resamples for both waveforms. Gene drops: 100 random + trio, point statistics (H279.1: the
coupling difference; KO: median |r| of all pairs). Trio alone: Dataset gate on the arm, statistic only if it passes.

    python p1_a279.py [--check]
Output: ../output/part1/a279.csv; ../output/part1_tables/a279_core_pairs.csv, a279_core_compartments.csv
"""
import argparse
import importlib.util
import itertools
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402

ENTRY = "a279"
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a279.py"
OUTA = PAPER / "exploratory/A279_aorta_pairwise_truth/output"
spec = importlib.util.spec_from_file_location("p1_a279_orig", PAPER / "exploratory/A279_aorta_pairwise_truth/scripts/score_a279.py")
M = importlib.util.module_from_spec(spec)
sys.modules["p1_a279_orig"] = M
spec.loader.exec_module(M)
A272, wrap, SV, ARMS = M.A272, M.wrap, M.SV, M.ARMS
WF = cs.gtex()


def untimed_pairs(t, need, wf, leave=None):
    if leave is not None:
        t = t[t.donor != leave].reset_index(drop=True)
    z, v = A272.atlas_states(t, wf)
    rows = {c: g.set_index("donor")["_i"] for c, g in t.assign(_i=np.arange(len(t))).groupby("compartment")}
    out = {}
    for a, b in itertools.combinations(sorted(rows), 2):
        d = rows[a].index.intersection(rows[b].index)
        if len(d) >= need - (leave is not None):
            out[(a, b)] = A272.pair_r(z, v, dict(ia=rows[a].loc[d].to_numpy(), ib=rows[b].loc[d].to_numpy(), s=np.zeros(len(d), int)))
    return out


def timed(t, don, need, wf, leave=None):
    if leave is not None:
        t = t[t.donor != leave]
    return cs.Dataset(t, anchors=None, donors=don, waveform=wf, min_donors=need - (leave is not None)).timed_loadings("zt")


def arm_tables(arm, t, don, need, wf, jackknife=True):
    """(compartments, pairs, H279.1 (diff, jackknife sd) or None) of one arm for waveform wf."""
    T = timed(t, don, need, wf)
    comp = pd.DataFrame([dict(arm=arm, compartment=c, n=r["n"], timed_phase_h=r["clock_offset_h"], timed_coupling=r["coupling"])
                         for c, r in T.iterrows()])
    U = untimed_pairs(t, need, wf)
    jk, jh = {k: [] for k in U}, []
    if jackknife:
        for m in sorted(t.donor.unique()):
            Tm, Um = timed(t, don, need, wf, leave=m), untimed_pairs(t, need, wf, leave=m)
            for (a, b) in U:
                if (a, b) in Um and a in Tm.index and b in Tm.index:
                    jk[(a, b)].append(wrap(np.angle(Um[(a, b)]) * HOURS - wrap(Tm.clock_offset_h[a] - Tm.clock_offset_h[b])))
            if "Macrophage" in Tm.index and all(c in Tm.index for c in SV):
                jh.append(Tm.coupling["Macrophage"] - min(Tm.coupling[c] for c in SV))
    rows = []
    for (a, b), r in U.items():
        td, ud = float(wrap(T.clock_offset_h[a] - T.clock_offset_h[b])), float(np.angle(r) * HOURS)
        e = np.asarray(jk[(a, b)], float)
        n = len(e)
        rows.append(dict(arm=arm, a=a, b=b, r_abs=abs(r), untimed_h=ud, timed_h=td, error_h=float(wrap(ud - td)),
                         jackknife_sd_h=float(np.sqrt((n - 1) / n * np.sum((e - e.mean()) ** 2))) if n > 1 else np.nan,
                         min_timed_coupling=float(min(T.coupling[a], T.coupling[b]))))
    pairs = pd.DataFrame(rows)
    pairs["rhythmic"] = pairs.min_timed_coupling >= M.MIN_COUPLING
    h1 = None
    if "Macrophage" in T.index and all(c in T.index for c in SV):
        d = float(T.coupling["Macrophage"] - min(T.coupling[c] for c in SV))
        if jackknife and jh:
            e = np.asarray(jh)
            n = len(e)
            h1 = (d, float(np.sqrt((n - 1) / n * np.sum((e - e.mean()) ** 2))))
        else:
            h1 = (d, np.nan)
    return comp, pairs, h1


def prim(pairs):
    g = pairs[(pairs.arm == "aligned") & pairs.rhythmic]
    return g


def err(g):
    return float(g.error_h.abs().median()) if len(g) else np.nan


def ko_stat(pairs):
    return float(pairs[pairs.arm == "bmal1_ko"].r_abs.median())


def run_all(tabs, wf, jackknife=True):
    out = {a: arm_tables(a, *tabs[a], wf, jackknife) for a in ARMS}
    comp = pd.concat([o[0] for o in out.values()], ignore_index=True)
    pairs = pd.concat([o[1] for o in out.values()], ignore_index=True)
    return comp, pairs, out["aligned"][2]


def boot_tabs(tabs, arm, rng):
    """Library bootstrap of one arm (libraries relabelled so duplicates are distinct donors)."""
    t, don, need = tabs[arm]
    libs = sorted(t.donor.unique())
    pick = rng.choice(libs, len(libs))
    parts, rows = [], []
    for i, l in enumerate(pick):
        x = t[t.donor == l].copy()
        x["donor"] = f"b{i}"
        parts.append(x)
        d = don.loc[[l]].copy()
        d.index = [f"b{i}"]
        rows.append(d)
    return pd.concat(parts, ignore_index=True), pd.concat(rows), need


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(1)
    warnings.filterwarnings("ignore")
    tabs = M.tables()
    st = pd.read_csv(OUTA / "tests.csv").set_index("scope")
    sp = pd.read_csv(OUTA / "pairs.csv")
    prim_row = st.loc["primary: aligned, both timed-rhythmic"]
    ko_row = st.loc["bmal1_ko, all pairs"]
    # ---- 14 genes, reproduced first ---------------------------------------------------------------------------------
    comp14, pairs14, h14 = run_all(tabs, WF)
    mm = pairs14.merge(sp, on=["arm", "a", "b"], suffixes=("", "_s"))
    e14, k14 = err(prim(pairs14)), ko_stat(pairs14)
    d14, sd14 = h14
    ok = (len(mm) == len(sp) and all(np.allclose(mm[c], mm[c + "_s"], atol=1e-6) for c in ("r_abs", "untimed_h", "timed_h"))
          and np.isclose(e14, prim_row.median_abs_error_h, atol=1e-6, rtol=0) and len(prim(pairs14)) == int(prim_row.pairs)
          and np.isclose(k14, ko_row.median_r_abs, atol=1e-6, rtol=0)
          and np.isclose(d14, prim_row.H279_1_diff, atol=1e-6, rtol=0) and np.isclose(sd14, prim_row.H279_1_jackknife_sd, atol=1e-6, rtol=0))
    if not ok and not chk:
        rows = [dict(id=i, figure_panel="Fig 2c", statistic=s, value_14_stored=v, value_14=x, reproduced=False,
                     notes="14-gene values not reproduced to 1e-6", script=SCRIPT)
                for i, s, v, x in (("A279_V279.1", "median |error| aligned", prim_row.median_abs_error_h, e14),
                                   ("A279_H279.1", "macrophage minus smallest SV coupling", prim_row.H279_1_diff, d14),
                                   ("A279_Bmal1KO", "KO median |r|", ko_row.median_r_abs, k14))]
        C.write_rows(ENTRY, rows)
        raise SystemExit("A279 not reproduced")
    # ---- core ---------------------------------------------------------------------------------------------------------
    core, trio = C.core(WF), C.trio(WF)
    compc, pairsc, hc = run_all(tabs, core)
    ec, kc = err(prim(pairsc)), ko_stat(pairsc)
    dc, sdc = hc
    # ---- trio contribution -------------------------------------------------------------------------------------
    key = ["arm", "a", "b"]
    m = pairs14[key + ["error_h", "rhythmic"]].merge(pairsc[key + ["error_h", "rhythmic"]], on=key, suffixes=("_14", "_c"))
    ma = m[m.arm == "aligned"].reset_index(drop=True)
    rng = np.random.default_rng([336, 1])
    nb = 3 if chk else 200
    dv1 = []
    for _ in range(nb):
        b = ma.iloc[rng.integers(0, len(ma), len(ma))]
        dv1.append(err(b[b.rhythmic_14].rename(columns={"error_h_14": "error_h"})) - err(b[b.rhythmic_c].rename(columns={"error_h_c": "error_h"})))
    dh, dk = [], []
    for _ in range(nb):
        ba = boot_tabs(tabs, "aligned", rng)
        bk = boot_tabs(tabs, "bmal1_ko", rng)
        h = [arm_tables("aligned", *ba, w, jackknife=False)[2] for w in (WF, core)]
        dh.append(h[0][0] - h[1][0])
        k = [ko_stat(arm_tables("bmal1_ko", *bk, w, jackknife=False)[1]) for w in (WF, core)]
        dk.append(k[0] - k[1])
    # ---- gene drops --------------------------------------------------------------------------------------------
    dl, tl = C.drops(WF, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    ve, vh, vk = {}, {}, {}
    for lab, w in dl:
        _, pr, h = run_all(tabs, w, jackknife=False)
        ve[lab], vh[lab], vk[lab] = err(prim(pr)), (h[0] if h else np.nan), ko_stat(pr)
    assert np.isclose(ve[tl], ec, equal_nan=True) and np.isclose(vh[tl], dc) and np.isclose(vk[tl], kc)
    dwf = dict(dl)
    # ---- trio alone --------------------------------------------------------------------------------------------
    gate, gp = "not computable", ""
    tv = dict(e=np.nan, h=np.nan, k=np.nan)
    try:
        g = cs.Dataset(tabs["aligned"][0], anchors=None, donors=tabs["aligned"][1], waveform=trio, name="trio",
                       min_donors=tabs["aligned"][2]).gate(n_perm=50 if chk else 1000, seed=336)
        gate = bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05)
        gp = f"trio gate (aligned arm) p_pairing {g.p_pairing:.4f}, p_acrophase {g.p_acrophase:.4f}"
        if gate:
            _, pt, ht = run_all(tabs, trio, jackknife=False)
            tv = dict(e=err(prim(pt)), h=ht[0] if ht else np.nan, k=ko_stat(pt))
    except Exception as e:  # noqa: BLE001
        gp = f"trio gate not computable: {type(e).__name__}: {str(e)[:100]}"
    if chk:
        print("check ok")
        return
    T = HERE.parents[1] / "output" / "part1_tables"
    T.mkdir(parents=True, exist_ok=True)
    pairsc.to_csv(T / "a279_core_pairs.csv", index=False)
    compc.to_csv(T / "a279_core_compartments.csv", index=False)
    n14, nc = len(prim(pairs14)), len(prim(pairsc))
    # criteria on the core
    h_pass = bool(dc < 0 and dc + 1.645 * sdc < 0)
    h14_pass = bool(d14 < 0 and d14 + 1.645 * sd14 < 0)
    ko_c = compc[compc.arm == "bmal1_ko"].timed_coupling
    sv_ko = pairsc[(pairsc.arm == "bmal1_ko") & pairsc.a.isin(SV) & pairsc.b.isin(SV)].r_abs
    ko_ok = bool(ko_c.max() <= 0.20 and (sv_ko >= 0.5).all() and len(sv_ko) == 3)
    ko_line = bool(round(ko_c.max(), 2) <= 0.20 and (sv_ko >= 0.5).all() and len(sv_ko) == 3)  # added after the run (rounding rule); see notes of the row
    ko_ok14 = bool(comp14[comp14.arm == "bmal1_ko"].timed_coupling.max() <= 0.20)
    common = dict(figure_panel="Fig 2c", trio_gate_passes=gate, blind="core blind", script=SCRIPT)
    defs = [
        ("A279_V279.1", "mouse aorta median |untimed - timed| offset (h), aligned timed-rhythmic pairs", -1, e14, prim_row.median_abs_error_h, ec,
         "median |error| <= 1.0 h over aligned pairs with timed coupling >= 0.5", "holds" if (np.isfinite(ec) and ec <= 1.0) else "falls", dv1, ve,
         "200 pair bootstraps over the 6 aligned pairs (rng [336, 1]), shared", tv["e"],
         f"primary pairs: 14 genes {n14}, core {nc}"),
        ("A279_H279.1", "aorta macrophage timed coupling minus smallest stromal-vascular coupling", -1, d14, prim_row.H279_1_diff, dc,
         "diff < 0 and jackknife 90 % upper bound < 0 (original: 14 genes FAIL)", "holds" if h_pass else "falls", dh, vh,
         "200 library bootstraps of the aligned arm (rng [336, 1]), shared",
         tv["h"], f"core diff {dc:.4f}, jackknife SD {sdc:.4f}, 90 % [{dc - 1.645 * sdc:.4f}, {dc + 1.645 * sdc:.4f}]; 14 genes pass={h14_pass} (as stored)"),
        ("A279_Bmal1KO", "Bmal1 KO (descriptive): median |r| of all pairs, no timed coupling", 1, k14, ko_row.median_r_abs, kc,
         "descriptive; operationalised: KO timed coupling <= 0.20 everywhere and the three stromal-vascular pairs |r| >= 0.5",
         "holds" if ko_ok else ("at the line" if ko_line else "falls"), dk, vk, "200 library bootstraps of the KO arm (rng [336, 1]), shared", tv["k"],
         f"KO max timed coupling core {ko_c.max():.3f} (14 genes ok {ko_ok14}); SV-pair |r| core {', '.join(f'{x:.3f}' for x in sorted(sv_ko))}")]
    rows = []
    for sid, name, direc, v14, vs, vc, crit, status, dd, dv, draws, tvv, note in defs:
        pr = C.drop_percentiles(dv, tl, dwf, direction=direc)
        rows.append(dict(id=sid, statistic=name, direction=direc, value_14_stored=vs, value_14=v14, reproduced=True, value_core=vc,
                         core_criterion=crit, core_status=status, diff_14_minus_core=v14 - vc,
                         diff_lo90=float(np.nanquantile(dd, .05)), diff_hi90=float(np.nanquantile(dd, .95)), diff_draws=draws, **pr,
                         trio_alone_value=tvv, reading=C.reading(status, pr["trio_drop_below_q10"]),
                         notes=f"{note}; {gp}; drops with an undefined statistic excluded", **common))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws", "figure_panel"]).round(4).T.to_string())
    print(df.notes.tolist())


if __name__ == "__main__":
    main()
