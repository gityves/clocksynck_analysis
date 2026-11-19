"""A335 part 1a: A257 (night-shift muscle) rerun with the core loop (GTEx minus PER1, PER2, NR1D1).

Same data, seeds and statistics as `A257.../compute_a257.py`; only the waveform changes. Companions: A257's pseudo sets
(seeds [257, 3] / [257, 4], 14-gene waveform) restricted to their genes at the core positions (A329 rule).
    python a335_1a_nightshift.py [--check]     # --check: 3 sets, assertions only, nothing saved or printed
Output: ../output/a335_1a_nightshift.csv (14-gene and core columns side by side), a335_1a_per_subject.csv.
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
spec = importlib.util.spec_from_file_location("a257", PAPER / "exploratory/A257_nightshift_muscle/scripts/compute_a257.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
DROP = ["PER1", "PER2", "NR1D1"]


def analyse(smp, y, wf, sets, rho_name):
    """H257.1-3 for waveform wf with companion waveforms sets = {kind: [Waveform]}; returns dict and per-subject frame."""
    z = A.zstates(y, wf)
    d = A.deltas(smp, z)
    ctl = {s: v for (s, c), v in d.items() if c == "ctl"}
    both = sorted(s for s in ctl if (s, "mis") in d)
    dc, dm = np.array([ctl[s] for s in both]), np.array([d[(s, "mis")] for s in both])
    al = np.array(list(ctl.values()))
    comp = {}
    for kind, ws in sets.items():
        comp[kind] = []
        for w in ws:
            dd = A.deltas(smp, A.zstates(y, w))
            a = np.array([v for (s, c), v in dd.items() if c == "ctl"])
            m = np.array([dd[(s, "mis")] for s in both])
            c2 = np.array([dd[(s, "ctl")] for s in both])
            comp[kind].append((A.rho(a), A.rho(c2) - A.rho(m)))
        comp[kind] = np.array(comp[kind])
    S1 = A.signs(len(al))
    nul1 = (S1 * (al * np.conj(A.CHORD)).real[None, :]).sum(1) / (len(al) * abs(A.CHORD) ** 2)
    obs1 = A.rho(al)
    ns = A.NB
    W = A.W
    shift = lambda a, b: float(np.angle(np.sum(a * np.conj(b))) / W)  # noqa: E731
    s_obs = shift(dm, dc)
    rng = np.random.default_rng([A.SEED, 5])
    bs = np.array([shift(dm[i], dc[i]) for i in (rng.integers(0, len(both), len(both)) for _ in range(ns))])
    bs = (bs - s_obs + 12) % 24 - 12 + s_obs
    lo, hi = float(np.quantile(bs, .05)), float(np.quantile(bs, .95))
    rc, rm = A.rho(dc), A.rho(dm)
    S3 = A.signs(len(both))
    pc, pmis = (dc * np.conj(A.CHORD)).real, (dm * np.conj(A.CHORD)).real
    nul3 = np.array([(np.where(sg < 0, pmis, pc).sum() - np.where(sg < 0, pc, pmis).sum()) / (len(both) * abs(A.CHORD) ** 2) for sg in S3])
    br = np.array([A.rho(dm[i]) / A.rho(dc[i]) for i in (rng.integers(0, len(both), len(both)) for _ in range(ns))])
    mq, rq = (float(np.quantile(comp[k][:, 1], .95)) for k in ("matched", "responsive"))
    p3 = float(np.mean(nul3 >= rc - rm))
    out = dict(
        h1_rho=obs1, h1_matched_q95=float(np.quantile(comp["matched"][:, 0], .95)),
        h1_responsive_q95=float(np.quantile(comp["responsive"][:, 0], .95)),
        s_h=s_obs, s_lo90=lo, s_hi90=hi, rho_ctl=rc, rho_mis=rm, R=rm / rc, R_lo90=float(np.quantile(br, .05)),
        R_hi90=float(np.quantile(br, .95)), p_ctl_minus_mis=p3, rho_diff=rc - rm, matched_diff_q95=mq, responsive_diff_q95=rq,
        h3_passes=bool(p3 < 0.05 and rc - rm > mq and rc - rm > rq))
    sub = pd.DataFrame(dict(subject=both, angle_h=np.angle(dm * np.conj(dc)) / W, amp_ratio=np.abs(dm) / np.abs(dc)))
    return out, sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    n = 3 if chk else A.N_SETS
    smp, y = A.load()
    WF = A.WF
    g = [x for x in WF.genes if x in y.columns]
    base = WF.subset(g)
    assert len(g) == len(WF.genes), "GTEx genes missing from array table"
    core = WF.drop(DROP)
    core_idx = [k for k, x in enumerate(base.genes) if x in core.genes]
    assert len(core_idx) == len(core.genes) == len(base.genes) - 3
    lm = y.mean(0)
    lm = lm[~lm.index.duplicated()]
    full_sets, core_sets = {}, {}
    for kind, tag in (("matched", 3), ("responsive", 4)):
        ws = P.pseudo_waveforms(lm, n, np.random.default_rng([A.SEED, tag]), base, kind)
        full_sets[kind] = ws
        core_sets[kind] = [Waveform(tuple(w.genes[k] for k in core_idx), core.amp, core.psi, None, None, f"core{j}") for j, w in enumerate(ws)]
    rf, sf = analyse(smp, y, base, full_sets, "full")
    if chk:
        assert len(sf) == 11 and all(len(w.genes) == 11 for w in core_sets["matched"])
        analyse(smp, y, core, core_sets, "core")
        print("check ok")
        return
    ref = pd.read_csv(PAPER / "exploratory/A257_nightshift_muscle/output/score_tests.csv")
    r1, r2, r3 = ref.iloc[0], ref.iloc[1], ref.iloc[2]
    exp = dict(h1_rho=r1.rho, h1_matched_q95=r1.matched_q95, h1_responsive_q95=r1.responsive_q95, s_h=r2.s, s_lo90=r2.lo90,
               s_hi90=r2.hi90, rho_ctl=r3.rho_ctl, rho_mis=r3.rho_mis, R=r3.R, R_lo90=r3.lo90, R_hi90=r3.hi90,
               p_ctl_minus_mis=r3.p_ctl_minus_mis, matched_diff_q95=r3.matched_diff_q95, responsive_diff_q95=r3.responsive_diff_q95)
    bad = {k: (rf[k], v) for k, v in exp.items() if not np.isclose(rf[k], v, atol=1e-6)}
    if bad:
        raise SystemExit(f"A257 not reproduced: stop before reading the core loop: {bad}")
    print("14-gene reproduction OK", flush=True)
    rc_, sc_ = analyse(smp, y, core, core_sets, "core")
    tab = pd.DataFrame([dict(statistic=k, original_14gene=rf[k], recomputed_14gene=rf[k], core_11gene=rc_[k]) for k in rf])
    tab["original_14gene"] = [exp.get(k, np.nan) for k in rf]
    tab["core_R_hi_below_1"] = rc_["R_hi90"] < 1
    tab["core_s_excludes_0"] = not (rc_["s_lo90"] <= 0 <= rc_["s_hi90"])
    tab.to_csv(OUT / "a335_1a_nightshift.csv", index=False)
    s = sf.merge(sc_, on="subject", suffixes=("_14gene", "_core"))
    s.to_csv(OUT / "a335_1a_per_subject.csv", index=False)
    pd.set_option("display.width", 200)
    print(tab.round(4).to_string())
    print(s.round(3).to_string())
    print("LOSS HOLDS:", bool(rc_["R_hi90"] < 1 and rc_["h3_passes"]), "SHIFT HOLDS:", bool(tab.core_s_excludes_0.iloc[0]))


if __name__ == "__main__":
    main()
