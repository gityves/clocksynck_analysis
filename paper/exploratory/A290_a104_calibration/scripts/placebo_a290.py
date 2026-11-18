"""A290 (registered in `paper/docs/registered_experiments.md`): placebo calibration of A271 Part B (H271.2, responsive
family) in the A285 form.

The A104 statistic, Spearman rho(own amplitude, coupling) over the 27 peripheral GTEx tissues (every tissue an anchor,
leave-one-out host), for the clock and the 200 responsive sets of A271 Part B (gene lists of `a104_sets.csv`, clock
amplitudes and phases). Draw 0 is the real pairing; draws 1..N permute the donor labels within each tissue, one
permutation per draw for every set. Per draw: clock rho, responsive median and q95, pi (share of sets below the clock).
One worker per draw; rows are appended per draw, so a rerun skips the draws already written.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A290_a104_calibration/scripts/placebo_a290.py [--draws 400] [--jobs 8]
    python placebo_a290.py --check   (draw 0 with the clock and 5 sets)
    python placebo_a290.py --score   (tests.csv from the draws written so far; the weekend run is stopped by the clock,
                                      see the registered deviation)
Outputs (`../output/`): draws.csv (one row per draw), tests.csv.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
SETS = PAPER / "exploratory/A271_responsive_companions/output/a104_sets.csv"
SEED = 290
META = ["donor", "compartment", "depth"]
WF = cs.gtex()
S = {}


def waveforms(n_sets=None):
    s = pd.read_csv(SETS).query("block == 'periphery' and kind in ('clock', 'responsive')").sort_values("set")
    s = s if n_sets is None else pd.concat([s[s.kind == "clock"], s[s.kind == "responsive"].head(n_sets)])
    return [(r.kind, Waveform(tuple(r.genes.split(";")), WF.amp, WF.psi, None, None, f"{r.kind}{r.set}")) for r in s.itertuples()]


def rho_of(pb, w):
    try:
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=S["comps"], waveform=w, name=w.name)
        return float(A104.rho(A104.per_tissue(ds)))
    except Exception:  # noqa: BLE001  (an undefined fit counts as below the clock, as registered)
        return np.nan


def draw(d):
    warnings.filterwarnings("ignore")
    pb = S["pb"]
    if d:
        rng = np.random.default_rng([SEED, d])
        pb = pb.copy()
        for idx in pb.groupby("compartment").indices.values():
            pb.iloc[idx, pb.columns.get_loc("donor")] = rng.permutation(pb.donor.to_numpy()[idx])
    r = np.array([rho_of(pb, w) for _, w in S["sets"]])
    clock, resp = r[0], r[1:]
    ok = resp[np.isfinite(resp)]
    pi = float(np.mean(np.where(np.isfinite(resp), resp, -np.inf) < clock)) if np.isfinite(clock) else 0.0
    return dict(draw=d, clock_rho=clock, resp_median=float(np.median(ok)) if len(ok) else np.nan,
                resp_q95=float(np.quantile(ok, 0.95)) if len(ok) else np.nan, resp_valid=len(ok), pi=pi,
                above_q95=bool(np.isfinite(clock) and len(ok) and clock > np.quantile(ok, 0.95)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--draws", type=int, default=400)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--score", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.score:
        return score(t0)
    from common import gtex
    sets = waveforms(5 if a.check else None)
    genes = sorted({g for _, w in sets for g in w.genes})
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    limit_threads(1)
    comps = sorted(c for c in pb.compartment.unique() if not c.startswith("Brain"))
    S.update(sets=sets, comps=comps, pb=pb[pb.compartment.isin(comps)][META + genes].reset_index(drop=True))
    print(f"{len(comps)} tissues, {len(sets)} sets, {len(genes)} genes, {time.time() - t0:.0f}s", flush=True)
    OUT.mkdir(exist_ok=True)
    if a.check:
        print(draw(0), f"{time.time() - t0:.0f}s", flush=True)
        return
    f = OUT / "draws.csv"
    done = set(pd.read_csv(f).draw) if f.exists() else set()
    todo = [d for d in range(a.draws + 1) if d not in done]
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, row in enumerate(p.imap_unordered(draw, todo)):
            pd.DataFrame([row]).to_csv(f, mode="a", header=not f.exists(), index=False)
            print(f"draw {row['draw']} pi {row['pi']:.3f} ({k + 1}/{len(todo)}) {time.time() - t0:.0f}s", flush=True)
    score(t0)


def score(t0):
    D = pd.read_csv(OUT / "draws.csv").sort_values("draw")
    r0, null = D[D.draw == 0].iloc[0], D[D.draw > 0]
    T = pd.DataFrame([dict(clock_rho=r0.clock_rho, resp_median=r0.resp_median, resp_q95=r0.resp_q95, pi_obs=r0.pi,
                           draws=len(null), v290_rejection=float(null.above_q95.mean()),
                           v290_calibrated=bool(null.above_q95.mean() <= 0.0685), null_rate_pi1=float((null.pi >= 1).mean()),
                           p_cal=float((1 + (null.pi >= r0.pi).sum()) / (1 + len(null))))])
    T["decidable"] = len(null) >= 60  # registered deviation: fewer draws, H290 not decidable
    T["h290_pass"] = (T.p_cal < 0.05).where(T.decidable)
    T.to_csv(OUT / "tests.csv", index=False)
    print(T.T.to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
