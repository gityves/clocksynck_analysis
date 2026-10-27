"""A139: waveform refit under hospital hours. Cells (window x n) x scenarios (S0, S1, S3, SW) x 100 replicates, cached
per (window, n, scenario, replicate) in output/cache/ so that a rerun resumes. Stages: `primary` (08-18 h, n = 80),
`rest` (the other cells, D139.3), `all`; `--quick` runs 3 replicates of the primary cell into output/cache_quick/.
Score with score_a139.py."""
from __future__ import annotations

import argparse
import pickle
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402


def task(cache: Path, window, n, scen, rep):
    from clocksync.config import limit_threads
    limit_threads(1)
    f = cache / f"{window}_{n}_{scen}_{rep:03d}.pkl"
    if f.exists():
        return "cached", 0.0
    t0 = time.time()
    try:
        res = lib.run_rep(window, n, scen, rep)
    except Exception:
        res = ("FAIL", dict(window=window, n=n, scenario=scen, rep=rep, error=traceback.format_exc(limit=3)))
    pickle.dump(res, open(f.with_suffix(".tmp"), "wb"))
    f.with_suffix(".tmp").rename(f)
    return ("fail" if isinstance(res[0], str) else "ok"), time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["primary", "rest", "all"], default="all")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    cache = lib.OUT / ("cache_quick" if a.quick else "cache")
    cache.mkdir(parents=True, exist_ok=True)
    cells = [(w, n) for w in lib.WINDOWS for n in lib.NS]
    prim = [("08-18h", 80)]
    cells = prim if (a.stage == "primary" or a.quick) else [c for c in cells if c not in prim] if a.stage == "rest" else prim + [c for c in cells if c not in prim]
    reps = 3 if a.quick else a.reps
    # heaviest first (E-b scenarios), replicates interleaved so that partial progress covers all scenarios
    jobs = [(w, n, s, r) for (w, n) in cells for r in range(reps) for s in ("SW", "S0", "S3", "S1")]
    print(f"{len(jobs)} tasks, {a.jobs} workers, cache {cache}", flush=True)
    t0, done = time.time(), 0
    import threading
    state = {"done": 0}

    def beat():  # heartbeat every 10 min (a task takes several minutes)
        while state["done"] < len(jobs):
            time.sleep(600)
            print(f"[{time.strftime('%F %T')}] heartbeat: {state['done']}/{len(jobs)} tasks, {(time.time()-t0)/60:.0f} min", flush=True)
    threading.Thread(target=beat, daemon=True).start()
    gen = Parallel(n_jobs=a.jobs, return_as="generator_unordered")(delayed(task)(cache, *j) for j in jobs)
    for st, dt in gen:
        done += 1
        state['done'] = done
        if done % 4 == 0 or done == len(jobs) or st == 'fail':
            el = time.time() - t0
            print(f"[{time.strftime('%F %T')}] {done}/{len(jobs)} tasks, {el/60:.1f} min, last {st} {dt:.0f}s", flush=True)
    print(f"done in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
