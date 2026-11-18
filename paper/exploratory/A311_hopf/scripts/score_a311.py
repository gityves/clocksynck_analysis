"""A311 scoring (registered; run once after the power note): orientation of the own variation of the 27 GTEx peripheral
tissues (clock and 100 expression-matched sets) and H311.1 / H311.2. Aggregates only (GTEx controlled access).

Rows, fit and matched sets as A294b (`A294_floor/scripts/score_a294.py`: same loading, seed [294, 1]). Per set and
tissue `clocksync.orientation.orientation` (clock 200 simulations, seed [311, tissue index]; sets 50).

    python score_a311.py [--jobs 4]
Outputs (`../output/`): tissues.csv (clock per tissue), sets.csv (per set and tissue), tests.csv.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync import orientation  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

sys.path.insert(0, str(HERE.parents[2] / "A294_floor" / "scripts"))
import score_a294 as A294  # noqa: E402

OUT = HERE.parents[1] / "output"
FLOOR = HERE.parents[2] / "A294_floor" / "output" / "tissues.csv"
G = {}


def run(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    kind, w = G["sets"][i]
    pb, comps = G["pb"], G["comps"]
    ds = cs.Dataset(pb[A294.META + list(w.genes)], anchors=comps, waveform=w, name=f"s{i}")
    rows = []
    for k, t in enumerate(comps):
        try:
            o = orientation.orientation(ds, t, n_sim=200 if kind == "clock" else 50, seed=[311, k])
        except (ValueError, np.linalg.LinAlgError) as e:
            o = dict(target=t, error=f"{type(e).__name__}: {e}"[:200])
        rows.append(dict(set=i, kind=kind, **o))
    out = pd.DataFrame(rows).rename(columns={"target": "tissue"})
    if kind == "clock":
        pt = A104.per_tissue(ds)
        out["S"] = out.tissue.map(pt.z2 * pt.V)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    pb = A294.load()
    limit_threads(1)
    META = A294.META
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    wf = cs.gtex()
    sets = [("clock", wf)] + [("matched", w) for w in pseudo_waveforms(lm, 100, np.random.default_rng([A294.SEED, 1]), wf, kind="matched")]
    keep = sorted(set(g for _, w in sets for g in w.genes))
    G.update(pb=pb[META + keep].copy(), sets=sets, comps=sorted(pb.compartment.unique()))
    del pb
    print(f"{len(sets) - 1} matched sets ({time.time() - t0:.0f}s)", flush=True)
    clock = run(0)
    floor = pd.read_csv(FLOOR, index_col=0)
    dev = float(np.max(np.abs(clock.set_index("tissue").own - floor.own.reindex(clock.tissue).to_numpy())))
    print(f"clock own vs A294b: max |diff| {dev:.3g} ({time.time() - t0:.0f}s)", flush=True)
    res = [clock]
    with Pool(a.jobs) as pool:
        for n, r in enumerate(pool.imap_unordered(run, range(1, len(sets))), 1):
            res.append(r)
            if n % 10 == 0:
                print(f"{n}/{len(sets) - 1} sets ({time.time() - t0:.0f}s)", flush=True)
    allr = pd.concat(res, ignore_index=True)
    allr.drop(columns=["S"]).to_csv(OUT / "sets.csv", index=False)
    c = clock.set_index("tissue")
    c["E"] = floor.own - floor.own_matched_median
    c["share_departure"] = c.share * c.own / c.E
    c.to_csv(OUT / "tissues.csv")
    rng = np.random.default_rng([311, 1])
    sh = c.share.to_numpy()
    med = float(np.median(sh))
    boot = np.median(sh[rng.integers(0, len(sh), (2000, len(sh)))], axis=1)
    lo, hi = np.quantile(boot, [0.05, 0.95])
    m = allr[(allr.kind == "matched")].dropna(subset=["share"])
    mmed = m.groupby("set").share.median()
    q95 = float(mmed.quantile(0.95))
    rho = spearmanr(c.share, c.S).statistic
    null = np.array([spearmanr(rng.permutation(sh), c.S).statistic for _ in range(2000)])
    p = float((1 + np.sum(null >= rho)) / 2001)
    tests = [dict(test="H311.1", median_share=med, lo90=lo, hi90=hi, matched_median_q95=q95,
                  matched_median_median=float(mmed.median()), matched_sets=int(len(mmed)),
                  percentile_vs_matched=float(np.mean(mmed < med)), passes=bool(lo > 0 and med > q95)),
             dict(test="H311.2", rho=float(rho), p_one_sided=p, passes=bool(rho > 0 and p < 0.05)),
             dict(test="descriptive", median_share_departure=float(c.share_departure.median()),
                  tissues_share_departure_gt_half=int((c.share_departure > 0.5).sum()),
                  tissues_share_gt_0=int((c.share > 0).sum()),
                  errors=int(allr.get("error", pd.Series(dtype=str)).notna().sum()))]
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    print(c[["S", "own", "E", "share", "share_departure"]].round(3).sort_values("S").to_string(), flush=True)
    print(pd.DataFrame(tests).round(4).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
