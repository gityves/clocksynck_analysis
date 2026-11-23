"""A336 part 1 wrapper: A273 V273.1 and V273.2 (mouse colon, untimed against timed pairwise offsets) on the 14 genes, the core
loop (Per1, Per2, Nr1d1 removed; stored under human symbols), the trio contribution, the gene-drop null and the trio alone.

Original: A273_pairwise_offset_truth/scripts/score_a273.py (imported, not edited; `arm_table`, `wrap`, constants). The two
functions that read the waveform (`untimed_pairs`, `timed`) are copied here with a `wf` argument. The primary set is
rebuilt for every waveform with the original rule (both compartments' timed coupling >= 0.5, arms C-DC + R-DC); the pair
counts are reported. No companions exist for these entries. 14 genes: tests.csv (primary row) and pairs.csv are
reproduced to 1e-6 before the core is computed. Trio contribution: 200 pair bootstraps over the all-pairs table of the two
arms (the statistic's unit; rng [336, 1]; the primary set re-filtered per waveform on each resample). Gene drops: 100
random + trio, point statistic only (no jackknife); V273.2 is ranked by |slope - 1| (direction -1). Trio alone: the
Dataset gate on the trio (donor x compartment, both arms stacked), the statistics only if it passes.

    python p1_a273.py [--check]
Output: ../output/part1/a273.csv; ../output/part1_tables/a273_core_pairs.csv, a273_core_compartments.csv
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
from common import data  # noqa: E402

ENTRY = "a273"
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a273.py"
OUTA = PAPER / "exploratory/A273_pairwise_offset_truth/output"
spec = importlib.util.spec_from_file_location("p1_a273_orig", PAPER / "exploratory/A273_pairwise_offset_truth/scripts/score_a273.py")
M = importlib.util.module_from_spec(spec)
sys.modules["p1_a273_orig"] = M
spec.loader.exec_module(M)
A272, wrap = M.A272, M.wrap
ARMS, SENS = M.ARMS, M.SENS
WF = cs.gtex()


def untimed_pairs(t, wf, leave=None):
    if leave is not None:
        t = t[t.donor != leave].reset_index(drop=True)
    z, v = A272.atlas_states(t, wf)
    rows = {c: g.set_index("donor")["_i"] for c, g in t.assign(_i=np.arange(len(t))).groupby("compartment")}
    out = {}
    for a, b in itertools.combinations(sorted(rows), 2):
        d = rows[a].index.intersection(rows[b].index)
        if len(d) >= M.MIN_MICE - (leave is not None):
            pr = dict(ia=rows[a].loc[d].to_numpy(), ib=rows[b].loc[d].to_numpy(), s=np.zeros(len(d), int))
            out[(a, b)] = A272.pair_r(z, v, pr)
    return out


def timed(t, donors, wf, leave=None):
    if leave is not None:
        t = t[t.donor != leave]
    return cs.Dataset(t, anchors=None, donors=donors, waveform=wf, min_donors=M.MIN_MICE - (leave is not None)).timed_loadings("zt")


def tables(wf, tabs, donors, jackknife=True):
    """(compartments, pairs) of every arm for waveform wf, as the original (jackknife SD optional)."""
    comp_rows, pair_rows = [], []
    for arm, t in tabs.items():
        T = timed(t, donors, wf)
        for c, r in T.iterrows():
            comp_rows.append(dict(arm=arm, compartment=c, n=r["n"], timed_phase_h=r["clock_offset_h"], timed_coupling=r["coupling"]))
        U = untimed_pairs(t, wf)
        jk = {k: [] for k in U}
        if jackknife:
            for m in sorted(t.donor.unique()):
                Tm, Um = timed(t, donors, wf, leave=m), untimed_pairs(t, wf, leave=m)
                for (a, b) in U:
                    if (a, b) in Um and a in Tm.index and b in Tm.index:
                        jk[(a, b)].append(wrap(np.angle(Um[(a, b)]) * HOURS - wrap(Tm.clock_offset_h[a] - Tm.clock_offset_h[b])))
        for (a, b), r in U.items():
            td, ud = float(wrap(T.clock_offset_h[a] - T.clock_offset_h[b])), float(np.angle(r) * HOURS)
            e = np.asarray(jk[(a, b)], float)
            n = len(e)
            pair_rows.append(dict(arm=arm, a=a, b=b, r_abs=abs(r), untimed_h=ud, timed_h=td, error_h=float(wrap(ud - td)),
                                  jackknife_sd_h=float(np.sqrt((n - 1) / n * np.sum((e - e.mean()) ** 2))) if n > 1 else np.nan,
                                  min_timed_coupling=float(min(T.coupling[a], T.coupling[b]))))
    comp, pairs = pd.DataFrame(comp_rows), pd.DataFrame(pair_rows)
    pairs["rhythmic"] = pairs.min_timed_coupling >= M.MIN_COUPLING
    return comp, pairs


def primary(pairs):
    return pairs[pairs.arm.isin(ARMS) & pairs.rhythmic]


def stats(g):
    """(median |error|, slope of untimed on timed) of a primary pair set; NaN when empty."""
    if not len(g):
        return np.nan, np.nan
    return float(g.error_h.abs().median()), float(np.sum(g.untimed_h * g.timed_h) / np.sum(g.timed_h ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(1)
    warnings.filterwarnings("ignore")
    pb, donors = data.load("colon_chrono")
    donors = donors.set_axis(donors.index.astype(str))
    tabs = {a: M.arm_table(pb, donors, a) for a in ARMS + SENS}
    stored_t = pd.read_csv(OUTA / "tests.csv").iloc[0]
    stored_p = pd.read_csv(OUTA / "pairs.csv")
    # ---- 14 genes, reproduced first ---------------------------------------------------------------------------------
    comp14, pairs14 = tables(WF, tabs, donors)
    e14, s14 = stats(primary(pairs14))
    mm = pairs14.merge(stored_p, on=["arm", "a", "b"], suffixes=("", "_s"))
    ok = (np.isclose(e14, stored_t.median_abs_error_h, atol=1e-6, rtol=0) and np.isclose(s14, stored_t.slope_untimed_on_timed, atol=1e-6, rtol=0)
          and len(primary(pairs14)) == int(stored_t.pairs) and len(mm) == len(stored_p)
          and all(np.allclose(mm[c], mm[c + "_s"], atol=1e-6) for c in ("r_abs", "untimed_h", "timed_h")))
    if not ok and not chk:
        rows = [dict(id=i, figure_panel="Fig 4a", statistic=s, direction=-1 if i.endswith("1") else 1, value_14_stored=v, value_14=x,
                     reproduced=False, notes="14-gene primary row not reproduced to 1e-6", script=SCRIPT)
                for i, s, v, x in (("A273_V273.1", "median |untimed - timed| offset", stored_t.median_abs_error_h, e14),
                                   ("A273_V273.2", "slope of untimed on timed", stored_t.slope_untimed_on_timed, s14))]
        C.write_rows(ENTRY, rows)
        raise SystemExit("A273 not reproduced")
    # ---- core -------------------------------------------------------------------------------------------------------
    core, trio = C.core(WF), C.trio(WF)
    compc, pairsc = tables(core, tabs, donors)
    ec, sc = stats(primary(pairsc))
    # ---- trio contribution: pair bootstrap, shared ---------------------------------------------------------------
    key = ["arm", "a", "b"]
    m = pairs14[key + ["error_h", "untimed_h", "timed_h", "rhythmic"]].merge(
        pairsc[key + ["error_h", "untimed_h", "timed_h", "rhythmic"]], on=key, suffixes=("_14", "_c"))
    m = m[m.arm.isin(ARMS)].reset_index(drop=True)
    rng = np.random.default_rng([336, 1])
    d_e, d_s = [], []
    for _ in range(5 if chk else 200):
        b = m.iloc[rng.integers(0, len(m), len(m))]
        x14 = b[b.rhythmic_14].rename(columns={"error_h_14": "error_h", "untimed_h_14": "untimed_h", "timed_h_14": "timed_h"})
        xc = b[b.rhythmic_c].rename(columns={"error_h_c": "error_h", "untimed_h_c": "untimed_h", "timed_h_c": "timed_h"})
        (a1, b1), (a2, b2) = stats(x14), stats(xc)
        d_e.append(a1 - a2), d_s.append(b1 - b2)
    # ---- gene drops (point statistics) ---------------------------------------------------------------------------
    dl, tl = C.drops(WF, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    dv_e, dv_s = {}, {}
    for lab, w in dl:
        _, pr = tables(w, tabs, donors, jackknife=False)
        a, b = stats(primary(pr))
        dv_e[lab], dv_s[lab] = a, (abs(b - 1) if np.isfinite(b) else np.nan)
    assert np.isclose(dv_e[tl], ec, equal_nan=True) and np.isclose(dv_s[tl], abs(sc - 1), equal_nan=True)
    dwf = dict(dl)
    # ---- trio alone ----------------------------------------------------------------------------------------------
    gate, gp, tv_e, tv_s = "not computable", "", np.nan, np.nan
    try:
        ds = cs.Dataset(pd.concat([tabs[a] for a in ARMS]), anchors=None, donors=donors, waveform=trio, name="trio", min_donors=M.MIN_MICE)
        g = ds.gate(n_perm=50 if chk else 1000, seed=336)
        gate = bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05)
        gp = f"trio gate p_pairing {g.p_pairing:.4f}, p_acrophase {g.p_acrophase:.4f}"
        if gate:
            _, pt = tables(trio, tabs, donors, jackknife=False)
            tv_e, tv_s = stats(primary(pt))
    except Exception as e:  # noqa: BLE001
        gp = f"trio gate not computable: {type(e).__name__}: {str(e)[:100]}"
    if chk:
        print("check ok")
        return
    # ---- tables for A293 H293.3 ----------------------------------------------------------------------------------
    T = HERE.parents[1] / "output" / "part1_tables"
    T.mkdir(parents=True, exist_ok=True)
    pairsc.to_csv(T / "a273_core_pairs.csv", index=False)
    compc.to_csv(T / "a273_core_compartments.csv", index=False)
    n14, nc = len(primary(pairs14)), len(primary(pairsc))
    common = dict(figure_panel="Fig 4a", species="mouse", n_drops=None, trio_gate_passes=gate, blind="core blind", script=SCRIPT)
    rows = []
    for sid, name, v14, vs, vc, dv, direc, crit, ok_core, dd in (
            ("A273_V273.1", "mouse colon median |untimed - timed| pairwise offset (h), primary timed-rhythmic pairs", e14, stored_t.median_abs_error_h, ec, dv_e, -1,
             "median |error| over the primary set <= 1.0 h", bool(ec <= 1.0) if np.isfinite(ec) else False, d_e),
            ("A273_V273.2", "slope of untimed on timed offsets (mouse colon)", s14, stored_t.slope_untimed_on_timed, sc, dv_s, -1,
             "slope in [0.7, 1.3]", bool(0.7 <= sc <= 1.3) if np.isfinite(sc) else False, d_s)):
        pr = C.drop_percentiles(dv, tl, dwf, direction=direc)
        status = "holds" if ok_core else "falls"
        rows.append(dict(id=sid, statistic=name, direction=direc, value_14_stored=vs, value_14=v14, reproduced=True, value_core=vc,
                         core_criterion=crit, core_status=status, diff_14_minus_core=v14 - vc,
                         diff_lo90=float(np.nanquantile(dd, .05)), diff_hi90=float(np.nanquantile(dd, .95)),
                         diff_draws="200 pair bootstraps over the C-DC + R-DC all-pairs table (rng [336, 1]), shared", **pr,
                         trio_alone_value=(tv_e if sid.endswith("1") else tv_s), reading=C.reading(status, pr["trio_drop_below_q10"]),
                         notes=(f"primary pairs: 14 genes {n14}, core {nc} (rule: both timed couplings >= 0.5 per waveform); "
                                f"core shares within 1 h: {float((primary(pairsc).error_h.abs() <= 1).mean()):.3f}; "
                                f"drops ranked by {'median |error|' if sid.endswith('1') else '|slope - 1|'} (direction -1: drop_q10/median shown on that scale, "
                                f"n_drops with an empty primary set excluded); {gp}"),
                         **{k: v for k, v in common.items() if k not in ("n_drops", "species")}))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws"]).round(4).T.to_string())
    print(df.notes.tolist())


if __name__ == "__main__":
    main()
