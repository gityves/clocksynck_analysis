"""A336 part 1, A257 (night-shift muscle, Figure 1d): two readouts, gene-drop null, trio alone, swap-test tie fixed.

Same data, seeds and statistics as A257 / A335 1a. Swap tests compare with a tolerance (null >= obs - 1e-12), so the
unswapped labelling always counts. Statistics: H257.1 aligned rotation rho; H257.2 shift s (h); H257.3 loss, read as
1 - R with R = rho_mis / rho_ctl (direction +1: more loss). Trio contribution: H257.2 / H257.3 on A257's own subject
bootstrap (rng [257, 5], shared by both waveforms), H257.1 on 200 subject bootstraps (rng [336, 1]). New blind layer
statistic (descriptive): per-subject trio-alone shift and amplitude ratio, with the H257.1 gate on the trio.

    python p1_a257.py [--check]
Outputs: ../output/part1/a257.csv, ../output/a257_trio_subjects.csv
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
sys.path.insert(0, str(HERE.parent))
spec = importlib.util.spec_from_file_location("a257", PAPER / "exploratory/A257_nightshift_muscle/scripts/compute_a257.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
import a336_common as C  # noqa: E402
from clocksync import paired as P  # noqa: E402

TOL = 1e-12
OUTD = HERE.parents[1] / "output"


def core_stats(smp, y, wf):
    z = A.zstates(y, wf)
    d = A.deltas(smp, z)
    ctl = {s: v for (s, c), v in d.items() if c == "ctl"}
    both = sorted(s for s in ctl if (s, "mis") in d)
    dc, dm = np.array([ctl[s] for s in both]), np.array([d[(s, "mis")] for s in both])
    al = np.array(list(ctl.values()))
    return al, dc, dm, both


def shift(a, b):
    return float(np.angle(np.sum(a * np.conj(b))) / A.W)


def stats(al, dc, dm):
    return dict(rho=A.rho(al), s=shift(dm, dc), loss=1 - A.rho(dm) / A.rho(dc))


def swap_p(al, dc, dm):
    S1 = A.signs(len(al))
    nul1 = (S1 * (al * np.conj(A.CHORD)).real[None, :]).sum(1) / (len(al) * abs(A.CHORD) ** 2)
    S3 = A.signs(len(dc))
    pc, pm = (dc * np.conj(A.CHORD)).real, (dm * np.conj(A.CHORD)).real
    nul3 = np.array([(np.where(sg < 0, pm, pc).sum() - np.where(sg < 0, pc, pm).sum()) / (len(dc) * abs(A.CHORD) ** 2) for sg in S3])
    o1, o3 = A.rho(al), A.rho(dc) - A.rho(dm)
    return float(np.mean(nul1 >= o1 - TOL)), float(np.mean(nul3 >= o3 - TOL))


def companions(smp, y, sets, both):
    out = []
    for w in sets:
        dd = A.deltas(smp, A.zstates(y, w))
        a = np.array([v for (s, c), v in dd.items() if c == "ctl"])
        m = np.array([dd[(s, "mis")] for s in both])
        c2 = np.array([dd[(s, "ctl")] for s in both])
        out.append((A.rho(a), A.rho(c2) - A.rho(m)))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    smp, y = A.load()
    g = [x for x in A.WF.genes if x in y.columns]
    base = A.WF.subset(g)
    assert len(g) == 14
    core, trio = C.core(base), C.trio(base)
    lm = y.mean(0)
    lm = lm[~lm.index.duplicated()]
    n = 3 if chk else A.N_SETS
    sets = {k: P.pseudo_waveforms(lm, n, np.random.default_rng([A.SEED, tag]), base, k) for k, tag in (("matched", 3), ("responsive", 4))}
    # ---- 14 genes: reproduce A257 -----------------------------------------------------------------------------------
    al, dc, dm, both = core_stats(smp, y, base)
    s14 = stats(al, dc, dm)
    ref = pd.read_csv(PAPER / "exploratory/A257_nightshift_muscle/output/score_tests.csv")
    stored = dict(rho=ref.iloc[0].rho, s=ref.iloc[1].s, loss=1 - ref.iloc[2].R)
    rep = {k: bool(np.isclose(s14[k], v, atol=1e-6)) for k, v in stored.items()}
    p1_14, p3_14 = swap_p(al, dc, dm)
    if not all(rep.values()):
        raise SystemExit(f"A257 not reproduced: {rep}")
    if chk:
        print("14-gene reproduction ok; tolerance p:", p1_14, p3_14)
    # ---- core ------------------------------------------------------------------------------------------------------
    alc, dcc, dmc, _ = core_stats(smp, y, core)
    sc = stats(alc, dcc, dmc)
    p1c, p3c = swap_p(alc, dcc, dmc)
    compc = {k: companions(smp, y, [C.restrict(w, base, core) for w in ws], both) for k, ws in sets.items()}
    # shift bootstrap and R bootstrap exactly as A257 / A335 (rng [257, 5]); shared indices for both waveforms
    nb = 20 if chk else A.NB
    rng = np.random.default_rng([A.SEED, 5])
    idx_s = [rng.integers(0, len(both), len(both)) for _ in range(nb)]
    idx_r = [rng.integers(0, len(both), len(both)) for _ in range(nb)]

    def sb(d_c, d_m, so):
        b = np.array([shift(d_m[i], d_c[i]) for i in idx_s])
        return (b - so + 12) % 24 - 12 + so

    bs14, bsc = sb(dc, dm, s14["s"]), sb(dcc, dmc, sc["s"])
    br14 = np.array([1 - A.rho(dm[i]) / A.rho(dc[i]) for i in idx_r])
    brc = np.array([1 - A.rho(dmc[i]) / A.rho(dcc[i]) for i in idx_r])
    rb = np.random.default_rng([336, 1])
    idx_a = [rb.integers(0, len(al), len(al)) for _ in range(200)]
    ba = np.array([A.rho(al[i]) - A.rho(alc[i]) for i in idx_a])
    # ---- gene-drop null --------------------------------------------------------------------------------------------
    dl, tl = C.drops(base)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:3]
    dv = {k: {} for k in ("rho", "s", "loss")}
    for lab, w in dl:
        st = stats(*core_stats(smp, y, w)[:3])
        for k in dv:
            dv[k][lab] = st[k]
    wfs = dict(dl)
    for k in dv:
        assert np.isclose(dv[k][tl], sc[k])
    # ---- trio alone -------------------------------------------------------------------------------------------------
    alt, dct, dmt, _ = core_stats(smp, y, trio)
    p1t, _ = swap_p(alt, dct, dmt)
    gate = bool(p1t < 0.05)
    st_t = stats(alt, dct, dmt)
    subj = pd.DataFrame(dict(subject=both, trio_angle_h=np.angle(dmt * np.conj(dct)) / A.W,
                             trio_amp_ratio=np.abs(dmt) / np.abs(dct)))
    if chk:
        print("check ok")
        return
    mq = {k: float(np.quantile(compc[k][:, 0], .95)) for k in compc}
    dq = {k: float(np.quantile(compc[k][:, 1], .95)) for k in compc}
    st1 = "holds" if (p1c < 0.05 and sc["rho"] > mq["matched"] and sc["rho"] > mq["responsive"]) else "falls"
    s_lo, s_hi = np.quantile(bsc, [0.05, 0.95])
    st2 = "holds" if not (s_lo <= 0 <= s_hi) else "falls"
    R_hi = 1 - np.quantile(brc, 0.05)
    diffc = A.rho(dcc) - A.rho(dmc)
    h3 = p3c < 0.05 and diffc > dq["matched"] and diffc > dq["responsive"]
    st3 = "holds" if (R_hi < 1 and h3) else "at the line" if (R_hi < 1 or h3) else "falls"
    rows = []
    for k, stat, status, crit, diff, draws, tv, note in (
        ("rho", "H257.1 aligned rotation rho", st1, "swap p < 0.05 and > core matched and responsive q95",
         ba, "200 subject bootstraps (rng [336, 1])", st_t["rho"],
         f"core swap p {p1c:.6f} (tolerance; 14 genes {p1_14:.6f}); core matched q95 {mq['matched']:.3f}, responsive q95 {mq['responsive']:.3f}"),
        ("s", "H257.2 shift of the muscle clock (h)", st2, "90 % bootstrap interval of s excludes 0",
         bs14 - bsc, "A257 subject bootstrap (rng [257, 5]), shared", st_t["s"],
         f"core s 90 % {s_lo:.2f} to {s_hi:.2f} h; 14-vs-core difference NOT BLIND (A335 per-subject table)"),
        ("loss", "H257.3 loss 1 - R (R = rho_mis / rho_ctl)", st3, "R 90 % upper < 1 and H257.3 swap / companion rule",
         br14 - brc, "A257 subject bootstrap (rng [257, 5]), shared", st_t["loss"],
         f"core R {1 - sc['loss']:.3f}, R upper {R_hi:.3f}; core swap p {p3c:.6f} (tolerance; 14 genes {p3_14:.6f}, A335 printed 0.0 from a floating tie); "
         f"rho_ctl - rho_mis {diffc:.3f} vs matched q95 {dq['matched']:.3f}, responsive q95 {dq['responsive']:.3f}")):
        pr = C.drop_percentiles(dv[k], tl, wfs, direction=1)
        rows.append(dict(id=f"A257_{stat.split()[0]}", figure_panel="Fig 1d", statistic=stat, direction=1,
                         value_14_stored=stored[k], value_14=s14[k], reproduced=True, value_core=sc[k], core_criterion=crit,
                         core_status=status, diff_14_minus_core=s14[k] - sc[k], diff_lo90=float(np.quantile(diff, 0.05)),
                         diff_hi90=float(np.quantile(diff, 0.95)), diff_draws=draws, **pr, trio_gate_passes=gate,
                         trio_alone_value=tv if gate else np.nan, reading=C.reading(status, pr["trio_drop_below_q10"]),
                         blind="core not blind (A335 1a)", notes=note, script="exploratory/A336_two_readouts/scripts/p1_a257.py"))
    C.write_rows("a257", rows)
    subj.assign(trio_gate_passes=gate, trio_gate_p=p1t).to_csv(OUTD / "a257_trio_subjects.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).drop(columns=["script", "notes"]).round(4).T.to_string())
    for r in rows:
        print(r["notes"])
    print(f"trio gate p {p1t:.5f} passes {gate}")
    print(subj.round(3).to_string() if gate else "trio gate failed: per-subject trio shifts not read")


if __name__ == "__main__":
    main()
