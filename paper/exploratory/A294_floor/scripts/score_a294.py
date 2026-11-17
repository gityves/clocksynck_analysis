"""A294 (registered in `paper/docs/registered_experiments.md`): is the own-variation floor of peripheral tissues specific
to the clock genes?

A104 rows and fit (GTEx fast deaths, tissues with >= 40 donors, periphery block, every tissue an anchor, leave-one-out
host). For the clock and 100 expression-matched pseudo-clock sets (`pseudo_waveforms(kind="matched")`, matching as A271
Part B): per tissue own variance 2 psi and shared variance S = z2 V. H294.1: median over tissues of
log(2 psi_clock / median_sets 2 psi) > 0 with a tissue-bootstrap 90 % interval above 0. H294.2 (if H294.1 passes):
Theil-Sen slope of log(excess) on log S_clock. Aggregates only.

    python score_a294.py [--n-sets 100] [--jobs 4] [--check]
Outputs (`../output/`): sets.csv (per set and tissue), tissues.csv, tests.csv.
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
from scipy.stats import theilslopes  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 294
META = ["donor", "compartment", "depth"]
S = {}


def load():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    return pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)


def fit(i):
    warnings.filterwarnings("ignore")
    kind, w = S["sets"][i]
    pb = S["pb"]
    comps = sorted(pb.compartment.unique())
    try:
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=f"s{i}")
        pt = A104.per_tissue(ds)
        pt = pt.assign(set=i, kind=kind, S=pt.z2 * pt.V, own=2 * pt.psi).rename_axis("tissue").reset_index()
        return pt
    except Exception as e:  # noqa: BLE001
        return pd.DataFrame([dict(set=i, kind=kind, error=f"{type(e).__name__}: {e}"[:200])])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sets", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    pb = load()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    wf = cs.gtex()
    n = 3 if a.check else a.n_sets
    sets = [("clock", wf)] + [("matched", w) for w in pseudo_waveforms(lm, n, np.random.default_rng([SEED, 1]), wf, kind="matched")]
    keep = sorted(set(g for _, w in sets for g in w.genes))
    S.update(pb=pb[META + keep].copy(), sets=sets)
    del pb
    print(f"{len(sets) - 1} matched sets, {len(keep)} genes {time.time() - t0:.0f}s", flush=True)
    with Pool(a.jobs) as pool:
        res = pd.concat(pool.map(fit, range(len(sets))), ignore_index=True)
    res.to_csv(OUT / "sets.csv", index=False)
    clock = res[res.kind == "clock"].set_index("tissue")
    rep = A104.rho(clock[["own_amplitude", "coupling"]])
    print(f"clock rho {rep:.3f} (A104 0.734); errors {res.get('error', pd.Series(dtype=str)).notna().sum()} {time.time() - t0:.0f}s", flush=True)
    if round(rep, 3) != 0.734:
        raise SystemExit("A104 not reproduced: stop")
    if a.check:
        return
    m = res[res.kind == "matched"].dropna(subset=["own"])
    floor = m.groupby("tissue").own.median()
    q95 = m.groupby("tissue").own.quantile(0.95)
    t = clock[["S", "own", "coupling", "own_amplitude"]].join(floor.rename("own_matched_median")).join(q95.rename("own_matched_q95"))
    t["log_ratio"] = np.log(t.own / t.own_matched_median)
    t["matched_S_median"] = m.groupby("tissue").S.median()
    t.to_csv(OUT / "tissues.csv")
    rng = np.random.default_rng(SEED)
    lr = t.log_ratio.to_numpy()
    boot = np.array([np.median(rng.choice(lr, len(lr))) for _ in range(2000)])
    lo, hi = np.quantile(boot, [0.05, 0.95])
    row = dict(test="H294.1", tissues=len(t), median_log_ratio=float(np.median(lr)), lo90=float(lo), hi90=float(hi),
               tissues_above_matched_q95=int((t.own > t.own_matched_q95).sum()),
               median_own_clock=float(t.own.median()), median_own_matched=float(t.own_matched_median.median()))
    row["passes"] = bool(row["median_log_ratio"] > 0 and lo > 0)
    rows = [row]
    if row["passes"]:
        e = (t.own - t.own_matched_median)
        ok = e > 0
        xs, ys = np.log(t.S[ok].to_numpy()), np.log(e[ok].to_numpy())
        sl = theilslopes(ys, xs)[0]
        bs = []
        for _ in range(2000):
            i = rng.integers(0, ok.sum(), ok.sum())
            bs.append(theilslopes(ys[i], xs[i])[0] if len(np.unique(xs[i])) > 2 else np.nan)
        blo, bhi = np.nanquantile(bs, [0.05, 0.95])
        rows.append(dict(test="H294.2", tissues=int(ok.sum()), slope=float(sl), lo90=float(blo), hi90=float(bhi),
                         signature="phase/amplitude (slope 1)" if blo <= 1 <= bhi and blo > 0 else
                         "additive (slope 0)" if blo <= 0 <= bhi and bhi < 1 else "undecided"))
    # descriptive: the same decomposition slope for the matched floor (does it scale with the clock's S?)
    xs = np.log(t.S.to_numpy())
    rows.append(dict(test="descriptive", slope_log_own_clock_on_log_S=float(theilslopes(np.log(t.own.to_numpy()), xs)[0]),
                     slope_log_floor_on_log_S=float(theilslopes(np.log(t.own_matched_median.to_numpy()), xs)[0])))
    pd.DataFrame(rows).to_csv(OUT / "tests.csv", index=False)
    print(t.round(3).sort_values("S").to_string(), flush=True)
    print(pd.DataFrame(rows).round(3).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
