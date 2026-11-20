"""A336 part 1, mouse feeding panels (Figure 5c): A297 H297.1, A322 H322.1, A324 H324.1-2 on 14 genes, the core loop
(Per1, Per2, Nr1d1 removed after the A111 QC; stored under human symbols) and the three-gene drops.

The waveform enters through `A111_feeding_controls/scripts/_common.waveform_without`, patched here to remove the QC
genes plus an extra set. 14 genes and core: the original `main()`s, unchanged, written to a scratch directory and
compared with / read beside the stored tables. Gene-drop null and trio contribution: the point statistics (A297 removal
amplitude ratio; A322 / A324 chord Q with the covariance of 200 library bootstraps, rng [336, 1], shared by every
waveform); no parametric bootstrap per drop. No companions exist for these entries. Trio alone: no gate exists for
timed mouse designs ("not computable").

    python p1_mouse_feeding.py [--check]
Output: ../output/part1/mouse_feeding.csv
"""
import argparse
import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(EXP / "A111_feeding_controls" / "scripts"))
import _common as CM  # noqa: E402
import a336_common as C  # noqa: E402


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


A297 = load("p1m_a297", EXP / "A297_feeding_model/scripts/score_a297.py")
A322 = load("p1m_a322", EXP / "A322_conflicting_drives/scripts/score_a322.py")
A324 = load("p1m_a324", EXP / "A322_conflicting_drives/scripts/score_a324.py")
ORIG = CM.waveform_without
EXTRA = set()
KEEP = None  # if set: use only these genes (trio alone)


def patched(genes):
    g = [x for x in CM.GENES if x not in set(genes) and x not in EXTRA and (KEEP is None or x in KEEP)]
    return CM.WF.subset(g)


CM.waveform_without = patched


def set_wave(extra=(), keep=None):
    global EXTRA, KEEP
    EXTRA, KEEP = set(extra), (set(keep) if keep is not None else None)


def mains(tmp):
    """Run A297, A322, A324 mains into tmp subfolders; return their test tables."""
    out = {}
    for key, m, sub in (("a297", A297, "a297"), ("a322", A322, "a322"), ("a324", A324, "a324")):
        m.OUT = Path(tmp) / sub
        m.main()
        out[key] = pd.read_csv(Path(tmp) / sub / "tests.csv")
    return out


def point_stats(n_boot=200, seed=(336, 1), reps=True):
    """A297 removal ratio and A322 / A324 chord Q (covariance from n_boot shared library bootstraps); also the
    bootstrap replicates of each (for the trio contribution): dict name -> (value, replicates)."""
    rng = np.random.default_rng(list(seed))
    A297.N_BOOT = n_boot
    res = {}
    s, _ = A297.liver_states("greenwell2019", "liver", "arm")
    a0, a1 = s[s.arm == "ad_libitum"], s[s.arm == "arrhythmic"]
    b0, b1 = A297.boot_amp(a0, rng), A297.boot_amp(a1, rng)
    res["A297_H297.1"] = (abs(A297.amplitude(a1) / A297.amplitude(a0)), np.abs(b1 / b0))
    for key, data, tissues, arms in (("A322_H322.1", "manella2021", ("Liver", "WAT", "Lung"), ("AlbCre_AL", "AlbCre_DF")),
                                     ("A324_H324.1", "xin2021", ("liver", "heart", "kidney", "VAT"), ("NRF", "DRF"))):
        r, R = [], []
        for t in tissues:
            s, _ = A297.liver_states(data, t, "arm")
            x0, x1 = (s[s.arm == a] for a in arms)
            r.append(A297.amplitude(x1) / A297.amplitude(x0))
            R.append(A297.boot_amp(x1, rng) / A297.boot_amp(x0, rng))
        r, R = np.array(r), np.array(R)
        P = np.linalg.inv(np.array([np.cov(A322.vec(R[k]).T) for k in range(len(tissues))]))
        Q = A322.fit(r, P)[0]
        res[key] = (Q, np.array([A322.fit(R[:, j], P)[0] for j in range(R.shape[1])]) if reps else None)
    A297.N_BOOT = 2000
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    A322.N_PAR = 50 if chk else 2000
    tmp = tempfile.mkdtemp(prefix="p1m_", dir="/tmp/claude-1000/-home-maxine-Documents-paychere-scCircadianOverview/45a469bb-a9f6-4040-a46a-2d3b8d13ecea/scratchpad")
    # ---- 14 genes reproduced ---------------------------------------------------------------------------------------
    set_wave()
    t14 = mains(tmp + "/f")
    st = {"a297": pd.read_csv(EXP / "A297_feeding_model/output/tests.csv"),
          "a322": pd.read_csv(EXP / "A322_conflicting_drives/output/tests.csv"),
          "a324": pd.read_csv(EXP / "A322_conflicting_drives/output/xin/tests.csv")}
    obs = lambda t: float(t[t.contrast.str.startswith("observed")].amp_ratio.iloc[0])  # noqa: E731
    v14 = {"A297_H297.1": obs(t14["a297"]), "A322_H322.1": float(t14["a322"].Q.iloc[0]), "A324_H324.1": float(t14["a324"].Q.iloc[0])}
    vst = {"A297_H297.1": obs(st["a297"]), "A322_H322.1": float(st["a322"].Q.iloc[0]), "A324_H324.1": float(st["a324"].Q.iloc[0])}
    rep = {k: bool(np.isclose(v14[k], vst[k], atol=1e-6)) for k in v14}
    if not chk:
        rep["A322 p"] = bool(np.isclose(t14["a322"].p_parametric.iloc[0], st["a322"].p_parametric.iloc[0]))
        rep["A324 H324.2"] = str(t14["a324"].all_above.iloc[1]) == str(st["a324"].all_above.iloc[1])
    if not all(rep.values()):
        raise SystemExit(f"mouse feeding 14-gene values not reproduced: {rep}")
    # ---- core ------------------------------------------------------------------------------------------------------
    set_wave(C.TRIO)
    tc = mains(tmp + "/c")
    vc = {"A297_H297.1": obs(tc["a297"]), "A322_H322.1": float(tc["a322"].Q.iloc[0]), "A324_H324.1": float(tc["a324"].Q.iloc[0])}
    # ---- point statistics: 14, core, drops, trio -------------------------------------------------------------------
    set_wave()
    p14 = point_stats()
    set_wave(C.TRIO)
    pc = point_stats()
    dl, tl = C.drops()
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    dv = {k: {} for k in p14}
    for lab, w in dl:
        set_wave(lab.split("+"))
        ps = point_stats(reps=False)
        for k in dv:
            dv[k][lab] = ps[k][0]
    set_wave(keep=C.TRIO)
    try:
        ptrio = point_stats()
    except Exception as e:  # noqa: BLE001
        ptrio = {k: (np.nan, None) for k in p14}
        print("trio alone not estimable:", type(e).__name__)
    set_wave()
    if chk:
        print("check ok")
        return
    wfs = dict(dl)
    a297 = tc["a297"]
    pred = a297.iloc[2]
    o = a297[a297.contrast.str.startswith("observed")].iloc[0]
    side = "inside" if pred.amp_lo90 <= o.amp_ratio <= pred.amp_hi90 else ("above" if o.amp_ratio > pred.amp_hi90 else "below")
    v297_14 = str(st["a297"][st["a297"].contrast.str.startswith("observed")].H297_1.iloc[0])
    verdict = {"A297_H297.1": (f"{'PASS' if side == 'inside' else 'FAIL'} ({side})", v297_14),
               "A322_H322.1": ("rejected" if not bool(tc["a322"].retained.iloc[0]) else "retained",
                               "rejected" if not bool(st["a322"].retained.iloc[0]) else "retained"),
               "A324_H324.1": ("rejected" if bool(tc["a324"].rejected.iloc[0]) else "retained",
                               "rejected" if bool(st["a324"].rejected.iloc[0]) else "retained")}
    rows = []
    for k, stat in (("A297_H297.1", "Greenwell arrhythmic / AL liver amplitude ratio vs the additive prediction"),
                    ("A322_H322.1", "chord-model Q, Manella liver / WAT / lung"), ("A324_H324.1", "chord-model Q, Xin four tissues")):
        same = verdict[k][0].split(" ")[0] == verdict[k][1].split(" ")[0] and (k != "A297_H297.1" or side in verdict[k][1])
        status = "holds" if same else "falls"
        assert np.isclose(dv[k][tl], pc[k][0])
        pr = C.drop_percentiles(dv[k], tl, wfs, direction=1)
        d = p14[k][1] - pc[k][1]
        rows.append(dict(id=k, figure_panel="Fig 5c", statistic=stat, direction=1, value_14_stored=vst[k], value_14=v14[k],
                         reproduced=True, value_core=vc[k], core_criterion="registered verdict unchanged",
                         core_status=status, diff_14_minus_core=p14[k][0] - pc[k][0], diff_lo90=float(np.quantile(d, 0.05)),
                         diff_hi90=float(np.quantile(d, 0.95)), diff_draws="200 shared library bootstraps within ZT (rng [336, 1])",
                         **pr, trio_gate_passes="not computable", trio_alone_value=np.nan,
                         reading=C.reading(status, pr["trio_drop_below_q10"]), blind="core blind",
                         notes=f"verdict core: {verdict[k][0]}; 14 genes: {verdict[k][1]}; trio-alone point value {ptrio[k][0]:.3f} (no gate)"
                         + (f"; core p_parametric {tc[k[:4].lower()].p_parametric.iloc[0]:.4f}" if k != "A297_H297.1" else
                            f"; core prediction {pred.amp_ratio:.3f} ({pred.amp_lo90:.3f}-{pred.amp_hi90:.3f}), core shift {o.shift_h:.2f} h"),
                         script="exploratory/A336_two_readouts/scripts/p1_mouse_feeding.py"))
    h2c, h2s = str(tc["a324"].all_above.iloc[1]), str(st["a324"].all_above.iloc[1])
    rows.append(dict(id="A324_H324.2", figure_panel="Fig 5c", statistic="every far tissue above the chord (Xin)",
                     value_14_stored=h2s, value_14=str(t14["a324"].all_above.iloc[1]), reproduced=True, value_core=h2c,
                     core_criterion="registered verdict unchanged", core_status="holds" if h2c == h2s else "falls",
                     reading="holds" if h2c == h2s else "falls: verdict changes", gene_drop="none",
                     notes=f"far tissues core: {tc['a324'].tissues_far.iloc[1]}; 14 genes: {st['a324'].tissues_far.iloc[1]}",
                     blind="core blind", script="exploratory/A336_two_readouts/scripts/p1_mouse_feeding.py"))
    df = C.write_rows("mouse_feeding", rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script"]).round(4).T.to_string())


if __name__ == "__main__":
    main()
