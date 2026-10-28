"""A158: why is the clocksync coupling drop in hepatocyte Bmal1 KO larger than the timed oracle's? Simulation.

Registered in paper/docs/registered_experiments.md (A158). DESIGN USED: the registration names the Jouffe 2022 tables
(`jouffe2022_bmal1ko`: whole-body KO, all samples at ZT12), but the observed values to reproduce (clocksync -0.293, oracle -0.186)
are D145.2 / A129 on Manella 2021 (hepatocyte KO vs Alb-Cre, AL, liver read on WAT + lung, 24 + 24 mice, ZT 0-22 every 2 h,
two per ZT), and a timed oracle does not exist at a single ZT. The simulated design is therefore the Manella one (see README).
Per replicate (seed [158, rep]): one control cohort (liver true coupling c0 = 0.97, z2 = 1.926 = observed Alb-Cre) and, per
multiplier x in {0.5, 0.6, 0.7, 0.8} (and 1.0, descriptive), one KO cohort (liver true coupling x c0), 24 mice each at the real ZTs,
tissues Liver (target), WAT and Lung (anchors, z = 1, coupling 0.9, unchanged). Host amplitude sqrt(14.0) (observed V = 14.0),
depth 3.2e6, simulator defaults otherwise. True change = c_KO - c0 (population coupling at the design's spread).
Arms: amp = "control" (KO liver z2 = control) and "observed" (KO z2 = 0.215 / 1.926 of control, as in the BLKO fit).
Estimators: clocksync (Dataset.fit coupling, anchors WAT + Lung; both genotypes share the ZT design, hence the same spread: the
difference equals d*), timed oracle (A145 `two_step` on the true ZT); clocksync with the liver as a third anchor (leave-one-out host).
Outputs: output/a158_raw.csv (replicate x arm x multiplier), a158_bias.csv, a158_reproduce.csv, run.log
"""
from __future__ import annotations

import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[v] = "1"
import multiprocessing as mp
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from lib158 import ANCH, Compartment, fit_cs, oracle, pop_v, psi_of, real_inputs, sim_table  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
N_REP, SEED, JOBS = 200, 158, 4
MULT = (0.5, 0.6, 0.7, 0.8, 1.0)
C0, Z2_CTRL, Z2_KO_OBS, C_ANCH, AMP, DEPTH = 0.97, 1.926, 0.215, 0.9, float(np.sqrt(14.0)), 3.2e6
OBS = dict(clocksync=-0.2934, oracle=-0.1859)
_G = {}


def init():
    pb, d, wf = real_inputs()
    ids = d.index[(d.genotype == "AlbCre") & (d.feeding == "AL")]
    _G.update(wf=wf, t=d.zt.reindex(ids).astype(float).to_numpy())
    _G["V1"] = pop_v(_G["t"])


def cohort(rng, c, z2):
    wf, t, V1 = _G["wf"], _G["t"], _G["V1"]
    comps = {k: Compartment(1.0, psi_of(C_ANCH, 1.0, AMP, V1)) for k in ANCH}
    comps["Liver"] = Compartment(np.sqrt(z2), psi_of(c, z2, AMP, V1))
    tab = sim_table(wf, comps, t, AMP, DEPTH, rng)
    r = fit_cs(tab, wf, ANCH)
    r["c_anchor"] = fit_cs(tab, wf, ANCH + ["Liver"])["c"]
    r["c_oracle"] = oracle(tab, wf, t)
    return r


def one(rep):
    rng = np.random.default_rng([SEED, rep])
    V1 = _G["V1"]
    ctrl = cohort(rng, C0, Z2_CTRL)
    rows = []
    for arm, z2 in (("control", Z2_CTRL), ("observed", Z2_KO_OBS)):
        for x in MULT:
            ko = cohort(rng, x * C0, z2)
            rows.append(dict(rep=rep, amp_arm=arm, mult=x, c_true_ctrl=C0, c_true_ko=x * C0, true_change=(x - 1) * C0,
                             cs_ctrl=ctrl["c"], cs_ko=ko["c"], or_ctrl=ctrl["c_oracle"], or_ko=ko["c_oracle"],
                             anc_ctrl=ctrl["c_anchor"], anc_ko=ko["c_anchor"], z2_ko=ko["z2"], boundary_ko=ko["boundary"], psi_ko=ko["psi"],
                             z2_ctrl=ctrl["z2"]))
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    rows = []
    with mp.get_context("spawn").Pool(JOBS, initializer=init) as pool:
        for i, r in enumerate(pool.imap_unordered(one, range(N_REP))):
            rows += r
            if i % 10 == 0:
                print(f"[{time.time() - t0:.0f} s] {i + 1}/{N_REP}", flush=True)
    R = pd.DataFrame(rows).sort_values(["amp_arm", "mult", "rep"])
    R["chg_cs"], R["chg_or"], R["chg_anc"] = R.cs_ko - R.cs_ctrl, R.or_ko - R.or_ctrl, R.anc_ko - R.anc_ctrl
    R.to_csv(OUT / "a158_raw.csv", index=False)
    g = R.groupby(["amp_arm", "mult"])
    B = g.agg(true_change=("true_change", "first"), cs_change=("chg_cs", "mean"), cs_sd=("chg_cs", "std"),
              oracle_change=("chg_or", "mean"), oracle_sd=("chg_or", "std"), anchor_change=("chg_anc", "mean"),
              cs_ko_mean=("cs_ko", "mean"), cs_ctrl_mean=("cs_ctrl", "mean"), or_ko_mean=("or_ko", "mean"),
              or_ctrl_mean=("or_ctrl", "mean"), c_true_ko=("c_true_ko", "first"),
              z2_ko_mean=("z2_ko", "mean"), share_boundary_ko=("boundary_ko", "mean"), share_z2_floor_ko=("z2_ko", lambda s: (s <= 1e-6).mean())).reset_index()
    B["cs_bias"] = B.cs_change - B.true_change
    B["oracle_bias"] = B.oracle_change - B.true_change
    B["anchor_minus_nonanchor_max_abs"] = g.apply(lambda d: float((d.chg_anc - d.chg_cs).abs().max())).to_numpy()
    B.to_csv(OUT / "a158_bias.csv", index=False)
    rep = []
    for arm, d in B.groupby("amp_arm"):
        d = d.sort_values("true_change")
        for est, col in (("clocksync", "cs_change"), ("oracle", "oracle_change")):
            # true change at which the mean estimate equals the observed value: linear interpolation (extrapolation flagged)
            xs, ys = d.true_change.to_numpy(), d[col].to_numpy()
            k = np.polyfit(xs, ys, 1)
            tc = float((OBS[est] - k[1]) / k[0])
            rep.append(dict(amp_arm=arm, estimator=est, observed_change=OBS[est], slope=k[0], intercept=k[1], true_change_reproducing=tc,
                            multiplier_reproducing=1 + tc / C0, inside_simulated_range=bool(xs.min() <= tc <= xs.max())))
    pd.DataFrame(rep).to_csv(OUT / "a158_reproduce.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(B.round(3).to_string()); print(pd.DataFrame(rep).round(3).to_string())
    print(f"done {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
