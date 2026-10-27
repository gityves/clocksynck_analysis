"""A140 Q1: simulation of the preprint's design (registered in paper/docs/registered_experiments.md, A140).

Run order (each step by `RUN_CAP_GB=30 paper/common/run.sh <this script> --step <s> --jobs N`, detached):
  timing   2 replicates per arm (A-D), K = 1, seed 500 without core genes: cost of one CYCLOPS fit and one gate
  gate     simulator gate: arm D, K = 1, seed 1,500 genes (with and without core genes), replicates 0-49, base seed 1500
  primary  H140.1/H140.2 (arm B, K = 1, seed 500, replicates 0-49, base seed 500) and H140.3 (clocksync gate in arms A, B, C,
           50 replicates, 08-18 h); scored by q1_score.py
  grid     descriptive grid: arms B, C, A, D (this order) x K in {1, 3} x seed {50, 150, 500, 1500} x with/without core genes,
           50 replicates; one file per (arm, K, replicate) in output/q1/grid/ as it finishes (resumable)

Outputs: output/q1/<step>/<arm>_K<K>_r<rep>.csv (one row per fit: readouts of the published pipeline) and
`gate_<step>.csv` rows (clocksync gate). Aggregates only (simulated data).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import pandas as pd
from joblib import Parallel, delayed

import q1_lib as Q


SCALES = (0.0, 0.1, 0.25, 0.5)


def task(step: str, arm: str, K: int, rep: int, combos: list, base_seed: int, threads: int, do_fit: bool, t0: float,
         scale: float = 1.0):
    from clocksync.config import limit_threads
    limit_threads(1)
    out = Q.OUT / step
    out.mkdir(parents=True, exist_ok=True)
    f = out / (f"{arm}_K{K}_r{rep:03d}.csv" if scale == 1.0 else f"{arm}_K{K}_s{scale:g}_r{rep:03d}.csv")
    if f.exists():
        return str(f)
    # data seed: arms differ, replicate and K identify the world (loadings shared by the arms of one replicate)
    P = Q.load_panel()
    wseed = base_seed + 1000 * K + rep
    L = Q.loadings(P, K, wseed)
    tab, t, S = Q.simulate_cohort(P, arm, K, wseed + 100_000 * (1 + Q.ARM_INDEX[arm]), L, scale=scale)
    row0 = dict(step=step, arm=arm, K=K, rep=rep, base_seed=base_seed, factor_scale=scale)
    t1 = time.time()
    g = Q.clocksync_gate(tab, seed=wseed)
    gate_s = time.time() - t1
    rows = []
    if not do_fit:
        rows.append({**row0, **g, "gate_s": gate_s, "n_seed": 0, "core": False})
    else:
        rk = Q.seed_ranking(P, K, L, wseed, scale)
        for n, core in combos:
            tmp = Q.cell_dir(f"{arm}{K}{rep}_{n}{int(core)}")
            t2 = time.time()
            try:
                ph = Q.run_fit(tab, P, Q.seed_set(P, rk, n, core), tmp, threads, seed=0)
                m = Q.pipeline_metrics(tab, ph, P, t)
                status = "ok"
            except Exception as e:  # a failed fit is a row, never silently dropped
                m, status = {}, f"error: {type(e).__name__}"
            rows.append({**row0, **g, "gate_s": gate_s, "n_seed": n, "core": core, "fit_s": time.time() - t2,
                         "status": status, **m})
            Q.cleanup(tmp)
            print(f"[{time.time() - t0:7.0f}s] {step} {arm} K{K} r{rep} seed{n} core{int(core)} {status} "
                  f"fit {time.time() - t2:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(f.with_suffix(".tmp"), index=False)
    f.with_suffix(".tmp").rename(f)
    return str(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--step", required=True, choices=["timing", "gate", "primary", "grid", "diag"])
    ap.add_argument("--jobs", type=int, default=1, help="parallel replicates (workers = jobs x threads, at most 16)")
    ap.add_argument("--threads", type=int, default=1, help="Julia workers per CYCLOPS fit")
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--arms", default="BCAD", help="grid: arm order")
    args = ap.parse_args()
    t0 = time.time()
    tasks = []  # (step, arm, K, rep, combos, base_seed, do_fit, factor scale)
    if args.step == "timing":
        tasks = [("timing", a, 1, r, [(500, False)], 99, True, 1.0) for a in "ABCD" for r in range(2)]
    elif args.step == "gate":
        tasks = [("gate", "D", 1, r, [(1500, False)], 1500, True, 1.0) for r in range(args.reps)]
    elif args.step == "primary":
        tasks = [("primary", "B", 1, r, [(500, False)], 500, True, 1.0) for r in range(args.reps)]
        tasks += [("primary", a, 1, r, [], 500, False, 1.0) for a in "AC" for r in range(args.reps)]
    elif args.step == "diag":  # post hoc diagnostic: arm D, seed 1,500, factor sd scaled (not registered)
        tasks = [("diag", "D", 1, r, [(1500, False)], 1500, True, sc) for sc in SCALES for r in range(args.reps)]
    else:
        combos = [(n, c) for n in Q.SEED_SIZES for c in (False, True)]
        for a in args.arms:
            for K in (1, 3):
                tasks += [("grid", a, K, r, combos, 2000, True, 1.0) for r in range(args.reps)]
    print(f"{len(tasks)} tasks, jobs {args.jobs} x threads {args.threads}", flush=True)
    done = Parallel(n_jobs=args.jobs, backend="loky", return_as="generator_unordered")(
        delayed(task)(s, a, K, r, c, b, args.threads, fit, t0, sc) for s, a, K, r, c, b, fit, sc in tasks)
    for i, d in enumerate(done, 1):
        print(f"[{time.time() - t0:7.0f}s] done {i}/{len(tasks)} {Path(d).name}", flush=True)
    print("finished", flush=True)


if __name__ == "__main__":
    main()
