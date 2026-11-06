"""A204: how do orderers produce "uncoupled tumour clocks" under working hours? Simulated fits per replicate.

Registered in paper/docs/registered_experiments.md (A204). Two task families, cached one pickle each under <out>/runs/
(a rerun resumes):
- cofe:  (window, truth, rep). 300 tumours and about 60 matched normals; COFE trained on the tumours, normals projected
         by COFE's test rule; clocksync's A138 H138.4 pairing statistic on the matched patients (1,000 permutations).
- tt:    (window, c_T, rep). TimeTeller trained on 144 timed whole-day healthy samples, tested on 108 tumours.
Scoring: score_a204.py.

    RUN_CAP_GB=10 paper/common/run.sh paper/exploratory/A204_orderer_uncoupling/scripts/compute_a204.py --jobs 8
    options: --reps N (default 30), --out DIR, --only cofe|tt, --mesor-sd (A204b: --only tt --mesor-sd 0)
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import multiprocessing as mp
import shutil
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent))
sys.path.insert(0, str(PAPER))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import paired as PD  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402

OUT = HERE.parent / "output"
HOURS = 24.0 / (2 * np.pi)
WF = cs.gtex()
GENES = list(WF.genes)
SEED = 204
N_PATIENTS, P_NORMAL = 300, 0.2
N_TRAIN, N_TEST = 144, 108
C_NORMAL, AMP_T, OFFSET_T_H = 0.9, 0.7, 3.0
N_PERM = 1000
WINDOWS = {"F": ((0.0, 24.0), 0.0), "W": ((8.0, 18.0), 1.0), "W_narrow": ((10.0, 16.0), 1.0)}
COFE_WINDOWS = ("F", "W")
COFE_TRUTHS = {"c0.9": 0.9, "c0.5": 0.5, "uncoupled": None}
TT_CS = (0.0, 0.25, 0.5, 0.75, 0.9)  # 0 = uncoupled
Z_T = AMP_T * np.exp(1j * OFFSET_T_H / HOURS)


def host_variance(window, jitter) -> float:
    """V = trace Cov(cos, sin) of the sampled host time (A145 lib.host_variance)."""
    r = np.random.default_rng(0)
    t = r.uniform(*window, 400_000) + (r.normal(0, jitter, 400_000) if jitter > 0 else 0.0)
    return float(np.trace(np.cov(np.column_stack([np.cos(t / HOURS), np.sin(t / HOURS)]).T)))


V_FULL = host_variance((0.0, 24.0), 0.0)


def psi_of(c: float, z: float) -> float:
    """psi giving coupling c at whole-day host variance (the biology is fixed; the window changes)."""
    return z * z * V_FULL * (1 - c) / (2 * c)


def tumour(c, mesor_sd: float = 0.5) -> Compartment:
    """Tumour compartment: coupled at c, or uncoupled (c None or 0) with the psi of c = 0.5."""
    if not c:
        return Compartment(z=Z_T, psi=psi_of(0.5, AMP_T), coupled=False, mesor_sd=mesor_sd)
    return Compartment(z=Z_T, psi=psi_of(c, AMP_T), mesor_sd=mesor_sd)


def log2(tab: pd.DataFrame) -> pd.DataFrame:
    from common.benchmark.data import log2_normalised

    return log2_normalised(tab.set_index("resection")[GENES].astype(float))


# ---- tasks ----------------------------------------------------------------------------------------------------------

def cofe_task(window: str, truth: str, rep: int, runs: str) -> dict:
    from clocksync.config import limit_threads

    limit_threads(1)
    f = Path(runs) / f"cofe_{window}_{truth}_{rep:03d}.pkl"
    if f.exists():
        return dict(key=f.stem, runtime_s=pd.read_pickle(f)["runtime_s"], failed=False)
    from common.benchmark.cofe import run_cofe

    (win, jit) = WINDOWS[window]
    comps = {"N": Compartment(z=1.0, psi=psi_of(C_NORMAL, 1.0), presence=P_NORMAL), "T": tumour(COFE_TRUTHS[truth])}
    tab, t = cs.simulate(WF, comps, n=N_PATIENTS, window_h=win, jitter_h=jit,
                         seed=np.random.SeedSequence([SEED, list(WINDOWS).index(window), list(COFE_TRUTHS).index(truth), rep]))
    t = pd.Series(t, index=[f"r{i:03d}" for i in range(N_PATIENTS)])
    T, N = tab[tab.compartment == "T"], tab[tab.compartment == "N"]
    LT, LN = log2(T), log2(N)
    keep = LT.columns[LT.std() > 0]
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_cofe(LT[keep], seed=rep, ncores=1, test=LN[keep])
    # clocksync: A138 H138.4 statistic on the matched patients
    m = N.resection.to_numpy()
    case, ref = T.set_index("resection").loc[m].reset_index(), N.reset_index(drop=True)
    bA, vA, bB, vB = PD.cohort_states(case, ref, counts=True, waveform=WF)
    obs = PD.cross(bA, vA, bB, vB)
    rng = np.random.default_rng([SEED, rep, 99])
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        p = rng.permutation(len(bB))
        null[k] = PD.cross(bA, vA, bB[p], vB[p])["corr"]
    out = dict(t=t, matched=list(m), phase_T=r["phase"], phase_N=r["test_phase"], interior_min=r["interior_min"],
               best_s=r["best_s"], cs_corr=obs["corr"], cs_offset_h=obs["offset_h"],
               cs_p=float((1 + np.sum(null >= obs["corr"])) / (1 + N_PERM)), cs_null_q95=float(np.quantile(null, 0.95)),
               runtime_s=time.time() - t0)
    pd.to_pickle(out, f)
    return dict(key=f.stem, runtime_s=out["runtime_s"], failed=False)


def tt_task(window: str, c: float, rep: int, runs: str, mesor_sd: float = 0.5) -> dict:
    from clocksync.config import limit_threads

    limit_threads(1)
    f = Path(runs) / f"tt_{window}_{c:.2f}_{rep:03d}.pkl"
    if f.exists():
        return dict(key=f.stem, runtime_s=pd.read_pickle(f)["runtime_s"], failed=False)
    from common.benchmark.timeteller import run_timeteller

    tr, ttr = cs.simulate(WF, {"H": Compartment(z=1.0, psi=psi_of(C_NORMAL, 1.0), mesor_sd=mesor_sd)}, n=N_TRAIN, window_h=(0.0, 24.0),
                          jitter_h=0.0, seed=np.random.SeedSequence([SEED, 7, rep]))  # same training set across cells
    (win, jit) = WINDOWS[window]
    te, tte = cs.simulate(WF, {"T": tumour(c, mesor_sd)}, n=N_TEST, window_h=win, jitter_h=jit,
                          seed=np.random.SeedSequence([SEED, 8, list(WINDOWS).index(window), TT_CS.index(c), rep]))
    Ltr, Lte = log2(tr), log2(te)
    Lte.index = [f"t{i}" for i in Lte.index]
    tmp = Path(runs) / "tmp" / f.stem
    t0 = time.time()
    try:
        res = run_timeteller(Ltr, pd.Series(ttr, index=Ltr.index), Lte, tmp, cores=1)
        err = ""
    except Exception as e:  # a failed fit is data
        res, err = None, f"{type(e).__name__}: {e}"[:300]
    shutil.rmtree(tmp, ignore_errors=True)
    out = dict(t=pd.Series(tte, index=Lte.index), pred=res, error=err, runtime_s=time.time() - t0)
    pd.to_pickle(out, f)
    return dict(key=f.stem, runtime_s=out["runtime_s"], failed=res is None)


# ---- main ---------------------------------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--only", choices=("cofe", "tt"))
    ap.add_argument("--mesor-sd", type=float, default=0.5, help="TimeTeller arm: per-gene tissue offset sd (A204b: 0)")
    args = ap.parse_args(argv)
    runs = Path(args.out) / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    todo = []
    if args.only in (None, "cofe"):
        todo += [("cofe", w, tr, r) for r in range(args.reps) for w in COFE_WINDOWS for tr in COFE_TRUTHS]
    if args.only in (None, "tt"):
        todo += [("tt", w, c, r) for r in range(args.reps) for w in WINDOWS for c in TT_CS]
    print(f"workers = {args.jobs} x 1 thread; {len(todo)} tasks ({args.reps} replicates)", flush=True)
    t0, done, last = time.time(), 0, time.time()
    with cf.ProcessPoolExecutor(args.jobs, mp_context=mp.get_context("spawn")) as ex:
        pend = {(ex.submit(cofe_task, w, x, r, str(runs)) if k == "cofe" else ex.submit(tt_task, w, x, r, str(runs), args.mesor_sd)):
                (k, w, x, r) for k, w, x, r in todo}
        while pend:
            fin, _ = cf.wait(pend, timeout=600, return_when=cf.FIRST_COMPLETED)
            for fu in fin:
                r = fu.result()
                pend.pop(fu)
                done += 1
                if done % 20 == 0 or r["failed"]:
                    print(f"[{time.time() - t0:.0f} s] {done}/{len(todo)} {r['key']} {r['runtime_s']:.0f} s failed={r['failed']}",
                          flush=True)
            if not fin or time.time() - last > 600:
                print(f"[{time.time() - t0:.0f} s] heartbeat: {done}/{len(todo)} done, {len(pend)} pending", flush=True)
                last = time.time()
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
