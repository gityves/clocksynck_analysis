"""A129 / D145.2 gap by simulation (roadmap 4.6, DESCRIPTIVE; no new hypothesis). Why does the clocksync drop (-0.293, Alb-Cre ->
BLKO) exceed the timed oracle's (-0.186)? Registered A158 simulated noise-only coupling loss and found the two estimators
agree there (gap <= 0.014); this adds the mechanisms of the roadmap list at the real design (Manella 2021 AL: real ZTs, genes
after QC, 24 mice per genotype, depth 3.2e6, Liver read on WAT + Lung host).

Control (every replicate, shared across mechanisms): liver population coupling 0.99, anchors 0.98, shared mouse phase jitter
SD 0.23 rad (calibrate_control.py: reproduces observed control clocksync 0.969 / oracle 0.934). Knockout mechanisms, 5 levels each:
  noise      liver psi raised so population coupling = x * 0.99, x in .9 .8 .7 .6 .5 (A158 mechanism)
  amplitude  liver amplitude z2 x f (f = .7 .5 .35 .2 .1), psi fixed at the control value (amplitude damping, phase intact)
  phasedisp  independent per-mouse liver phase noise SD .2 .4 .6 .8 1.0 rad (phase dispersion), amplitude and psi as control
  hostjitter liver stops sharing the mouse's host deviation (weight w = .75 .5 .25 0 and jitter x 2 in the KO mouse), amplitude as control
  hostspread KO-arm host (anchor) amplitude x .8 .6 .45 .3 .2, liver untouched (the host re-estimated per arm at a different spread)
  small_n    noise mechanism read at n = 12 mice (every other ZT), to size the noise floor of the gap
  noise_hj0  noise mechanism (x = .9 .8 .7 .6 .5) together with hostjitter w = 0 (the KO liver follows the ZT, not the mouse)
Estimators: clocksync (Dataset.fit coupling, anchors WAT + Lung; difference of independent cohorts at the same ZT design, so
equals d*) and the A145 timed oracle (two-step on true ZT). 200 replicates, seed 1290.
Outputs: output/mech_raw.csv, output/mech_summary.csv, output/mech_ratio.csv
"""
from __future__ import annotations
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[v] = "1"
import multiprocessing as mp, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from lib_a129sim import ANCH, Compartment, fit_cs, oracle, pop_v, psi_of, real_inputs, sim_table

OUT = Path(__file__).resolve().parent / "output"
N_REP, SEED, JOBS = 200, 1290, int(os.environ.get("SIM_JOBS", 6))
AMP, DEPTH, Z2C, C_LIV, C_ANC, JIT = float(np.sqrt(14.0)), 3.2e6, 1.926, 0.99, 0.98, 0.23
OBS = dict(clocksync=-0.2934, oracle=-0.1859)
LEVELS = dict(noise=(0.9, 0.8, 0.7, 0.6, 0.5), amplitude=(0.7, 0.5, 0.35, 0.2, 0.1), phasedisp=(0.2, 0.4, 0.6, 0.8, 1.0),
              hostjitter=(0.75, 0.5, 0.25, 0.0, -0.25), hostspread=(0.8, 0.6, 0.45, 0.3, 0.2), small_n=(0.9, 0.8, 0.7, 0.6, 0.5),
              noise_hj0=(0.9, 0.8, 0.7, 0.6, 0.5))
_G = {}


def init():
    pb, d, wf = real_inputs()
    ids = d.index[(d.genotype == "AlbCre") & (d.feeding == "AL")]
    _G.update(wf=wf, t=d.zt.reindex(ids).astype(float).to_numpy())
    _G["V1"] = pop_v(_G["t"])


def cohort(rng, mech="control", lev=None, t=None):
    wf, V1 = _G["wf"], _G["V1"]
    t = _G["t"] if t is None else t
    psi_a, psi_c = psi_of(C_ANC, 1.0, AMP, V1), psi_of(C_LIV, Z2C, AMP, V1)
    z2, psi, w, psd, am_a, jit = Z2C, psi_c, 1.0, 0.0, 1.0, JIT
    if mech in ("noise", "small_n", "noise_hj0"):
        psi = psi_of(lev * C_LIV, Z2C, AMP, V1)
        w = 0.0 if mech == "noise_hj0" else 1.0
    elif mech == "amplitude":
        z2 = Z2C * lev
    elif mech == "phasedisp":
        psd = lev
    elif mech == "hostjitter":
        w = lev
    elif mech == "hostspread":
        am_a = lev
    comps = {k: Compartment(1.0, psi_a) for k in ANCH}
    comps["Liver"] = Compartment(np.sqrt(z2), psi)
    tab = sim_table(wf, comps, t, AMP, DEPTH, rng, jitter_sd=jit, w={"Liver": w}, phase_sd={"Liver": psd},
                    amp_mult={"WAT": am_a, "Lung": am_a, "Liver": 1.0 if mech != "hostspread" else 1.0})
    r = fit_cs(tab, wf, ANCH)
    c_true = z2 * AMP**2 * V1 / (z2 * AMP**2 * V1 + 2 * psi)
    return dict(cs=r["c"], z2=r["z2"], oracle=oracle(tab, wf, t), c_true=c_true)


def one(rep):
    rng = np.random.default_rng([SEED, rep])
    ctrl = cohort(rng)
    ctrl12 = cohort(rng, t=_G["t"][::2])
    rows = []
    for mech, levs in LEVELS.items():
        for lev in levs:
            n12 = mech == "small_n"
            ko = cohort(rng, mech, lev, t=_G["t"][::2] if n12 else None)
            c0 = ctrl12 if n12 else ctrl
            rows.append(dict(rep=rep, mech=mech, level=lev, cs_ctrl=c0["cs"], cs_ko=ko["cs"], or_ctrl=c0["oracle"], or_ko=ko["oracle"],
                             c_true_ko=ko["c_true"], c_true_ctrl=c0["c_true"], z2_ko=ko["z2"]))
    return rows


def main():
    OUT.mkdir(exist_ok=True); t0 = time.time(); rows = []
    with mp.get_context("spawn").Pool(JOBS, initializer=init) as pool:
        for i, r in enumerate(pool.imap_unordered(one, range(N_REP))):
            rows += r
            if i % 10 == 0:
                print(f"[{time.time() - t0:.0f} s] {i + 1}/{N_REP}", flush=True)
    R = pd.DataFrame(rows).sort_values(["mech", "level", "rep"])
    R["chg_cs"], R["chg_or"] = R.cs_ko - R.cs_ctrl, R.or_ko - R.or_ctrl
    R["gap"] = R.chg_cs - R.chg_or
    R.to_csv(OUT / "mech_raw.csv", index=False)
    S = R.groupby(["mech", "level"]).agg(cs_ctrl=("cs_ctrl", "mean"), or_ctrl=("or_ctrl", "mean"), cs_ko=("cs_ko", "mean"), or_ko=("or_ko", "mean"),
                                         cs_change=("chg_cs", "mean"), or_change=("chg_or", "mean"), gap_mean=("gap", "mean"), gap_sd=("gap", "std"),
                                         cs_change_sd=("chg_cs", "std"), or_change_sd=("chg_or", "std"), c_true_change=("c_true_ko", "mean"),
                                         z2_ko=("z2_ko", "mean")).reset_index()
    S["ratio"] = S.cs_change / S.or_change
    S.to_csv(OUT / "mech_summary.csv", index=False)
    rat = []
    for mech, d in S.groupby("mech"):
        d = d.sort_values("or_change", ascending=False)
        x, yc = d.or_change.to_numpy(), d.cs_change.to_numpy()
        o = np.argsort(x); x, yc = x[o], yc[o]
        inside = bool(x.min() <= OBS["oracle"] <= x.max())
        cs_at = float(np.interp(OBS["oracle"], x, yc))
        k = np.polyfit(x, yc, 1)
        rat.append(dict(mech=mech, oracle_drop_matched=OBS["oracle"], cs_drop_at_matched_oracle=cs_at, ratio_at_match=cs_at / OBS["oracle"],
                        inside_range=inside, slope_cs_on_oracle=float(k[0]), observed_ratio=OBS["clocksync"] / OBS["oracle"]))
    pd.DataFrame(rat).to_csv(OUT / "mech_ratio.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(S.round(3).to_string()); print(pd.DataFrame(rat).round(3).to_string()); print(f"done {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
