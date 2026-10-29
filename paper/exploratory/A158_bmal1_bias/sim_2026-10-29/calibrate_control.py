"""Control calibration of the mechanism simulator: which shared-jitter SD and liver coupling reproduce the observed Alb-Cre
control values (clocksync 0.969, timed oracle 0.934)? Control cohort only, 40 replicates per cell, 4 workers.
Output: output/calibration_control.csv"""
from __future__ import annotations
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[v] = "1"
import itertools, multiprocessing as mp, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from lib_a129sim import ANCH, Compartment, fit_cs, oracle, pop_v, psi_of, real_inputs, sim_table

OUT = Path(__file__).resolve().parent / "output"
AMP, DEPTH, Z2C = float(np.sqrt(14.0)), 3.2e6, 1.926
_G = {}

def init():
    pb, d, wf = real_inputs()
    ids = d.index[(d.genotype == "AlbCre") & (d.feeding == "AL")]
    _G.update(wf=wf, t=d.zt.reindex(ids).astype(float).to_numpy()); _G["V1"] = pop_v(_G["t"])

def one(a):
    jit, c_liv, c_anc, rep = a
    wf, t, V1 = _G["wf"], _G["t"], _G["V1"]
    rng = np.random.default_rng([1290, rep, int(jit * 100), int(c_liv * 1000), int(c_anc * 1000)])
    comps = {k: Compartment(1.0, psi_of(c_anc, 1.0, AMP, V1)) for k in ANCH}
    comps["Liver"] = Compartment(np.sqrt(Z2C), psi_of(c_liv, Z2C, AMP, V1))
    tab = sim_table(wf, comps, t, AMP, DEPTH, rng, jitter_sd=jit)
    return dict(jitter_sd=jit, c_liver=c_liv, c_anchor=c_anc, rep=rep, cs=fit_cs(tab, wf, ANCH)["c"], oracle=oracle(tab, wf, t))

if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); t0 = time.time()
    grid = [(j, cl, ca, r) for j, cl, ca in itertools.product((0.0, 0.15, 0.3, 0.45), (0.97, 0.99), (0.9, 0.97)) for r in range(40)]
    with mp.get_context("spawn").Pool(4, initializer=init) as p:
        rows = p.map(one, grid, chunksize=8)
    D = pd.DataFrame(rows).groupby(["jitter_sd", "c_liver", "c_anchor"]).agg(cs=("cs", "mean"), cs_sd=("cs", "std"), oracle=("oracle", "mean"), oracle_sd=("oracle", "std")).reset_index()
    D.to_csv(OUT / "calibration_control.csv", index=False)
    print(D.round(3).to_string(), f"\n{time.time()-t0:.0f} s")
