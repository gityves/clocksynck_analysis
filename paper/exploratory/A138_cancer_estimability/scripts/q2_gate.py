"""A138 Q2 simulator gate: full arm without heterogeneity (intact clock over 24 h, patient phase SD 0), n = 300,
50 replicates; CoFE (genome-wide as published, and clock genes only) and CYCLOPS 2.0 (default seeds, 80 models).
Only the gate is run; no other arm, size or hypothesis.

Registered gate (registered_experiments.md, A138 Q2): "in the full arm without heterogeneity, CoFE and CYCLOPS recover
the truth to within their published GTEx-like error (CoFE clock genes median acrophase error <= 2 h at n = 300); if not,
Q2 stops and the simulator is fixed first."

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A138_cancer_estimability/scripts/q2_gate.py --jobs 4 [--reps 50]
Outputs (output/q2_gate/, one file per replicate and method, resumable): reps/<method>_r<rep>.csv, genes/<method>_r<rep>.csv.
Score with q2_gate_score.py.

Choices not in the registration (all minor):
- Panel: K = 10 leading PCs of the BRCA panel matrix (782 tumours), scores resampled as whole rows (no replacement);
  residual SD per gene from the BRCA residual; followers random among GTEx genes with p < 0.01 in >= 3 tissues; the
  CYCLOPS default seed genes present in BRCA are forced into the panel; panel fixed (seed 138), data seed per replicate.
- Clock genes: the 15 CoFE genes; waveform from clocksync's GTEx table (13 genes), NFIL3 and BHLHE40 from pooled GTEx
  cosinor. GTEx amplitudes are natural-log, converted to log2.
- Patient phase SD 0: the tumour clock state equals its sampling time (uniform 0-24 h).
- CoFE: as published (`common/benchmark/cofe.py`, 3 x 5-fold CV, 5 restarts, 20 sparsity values), on the simulated log2
  table without further filtering (the panel genes are BRCA genes that passed COFE's filter). CYCLOPS: template settings,
  80 models, default (30 ubiquity-cycler) seed genes, expression rank cutoff 1,500 of 2,000 genes, no covariate, input
  2^Y - 1; one Julia worker per fit.
- Gate statistic: per replicate the median over the 15 clock genes of |acrophase error| after oracle alignment of
  rotation and reflection (benchmark metric, `arrangement_error` against the cosinor of the true times); also reported
  against the true waveform. Pass: median over replicates <= 2 h. CYCLOPS bound: 2 h too (GTEx n = 80 full day: 1.73 h).
  Verdict "CoFE" = CoFE as published (genome-wide); clock-only CoFE reported beside it.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))

import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import q2_lib as Q  # noqa: E402

N = 300
BASE_SEED = 2138
METHODS = ("cyclops", "cofe_genome", "cofe_clock")  # longest first


def task(method: str, rep: int, n: int, t0: float, out: Path):
    from clocksync.config import limit_threads
    limit_threads(1)
    f = out / "reps" / f"{method}_r{rep:03d}.csv"
    if f.exists():
        return f.name, "cached"
    P = Q.load_panel()
    Y, t = Q.simulate(P, n, BASE_SEED + rep, arm="full", sd_h=0.0)
    tmp = Path(tempfile.mkdtemp(prefix=f"a138q2_{method}{rep}_"))
    t1 = time.time()
    try:
        ph, meta = Q.run_method(method, Y, seed=rep, tmp=tmp, threads=1)
        m, genes = Q.score(ph, Y, t, P)
        status = "ok"
    except Exception as e:  # a failed fit is a row, never dropped
        meta, m, genes, status = {}, {}, None, f"error: {type(e).__name__}: {str(e)[:120]}"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    row = dict(method=method, rep=rep, n=n, status=status, fit_s=time.time() - t1, **meta, **m)
    (out / "reps").mkdir(parents=True, exist_ok=True)
    (out / "genes").mkdir(parents=True, exist_ok=True)
    if genes is not None:
        genes.assign(method=method, rep=rep).to_csv(out / "genes" / f.name)
    pd.DataFrame([row]).to_csv(f.with_suffix(".tmp"), index=False)
    f.with_suffix(".tmp").rename(f)
    print(f"[{time.time() - t0:7.0f}s] {method} r{rep} {status} {row['fit_s']:.0f}s "
          f"clock_err_med {m.get('clock_err_med_h', float('nan')):.2f} h", flush=True)
    return f.name, status


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=1, help="parallel tasks (one thread each; at most 4 here)")
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--n", type=int, default=N)
    ap.add_argument("--out", default=str(Q.OUT))
    args = ap.parse_args()
    assert args.jobs <= 4
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    Q.load_panel()  # build once before the workers start
    t0 = time.time()
    tasks = [(m, r) for r in range(args.reps) for m in METHODS]
    print(f"{len(tasks)} tasks (n = {args.n}, {args.reps} replicates), jobs {args.jobs}", flush=True)
    stop = threading.Event()

    def beat():
        while not stop.wait(600):
            k = len(list((out / "reps").glob("*.csv"))) if (out / "reps").exists() else 0
            print(f"[{time.time() - t0:7.0f}s] heartbeat: {k}/{len(tasks)} tasks cached", flush=True)
    threading.Thread(target=beat, daemon=True).start()
    done = Parallel(n_jobs=args.jobs, backend="loky", return_as="generator_unordered")(
        delayed(task)(m, r, args.n, t0, out) for m, r in tasks)
    for i, (name, st) in enumerate(done, 1):
        print(f"[{time.time() - t0:7.0f}s] done {i}/{len(tasks)} {name} {st}", flush=True)
    stop.set()
    print("finished", flush=True)


if __name__ == "__main__":
    main()
