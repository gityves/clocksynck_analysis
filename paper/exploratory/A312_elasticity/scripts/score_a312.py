"""A312 (registered): amplitude elasticity of the 27 GTEx peripheral tissues to the persons' shared clock amplitude,
and its gradient with gain (H312.1). Aggregates only (GTEx controlled access).

A104 rows (fast deaths, tissues with >= 40 donors, periphery, every tissue an anchor). Per target tissue,
`clocksync.elasticity` on 20 halvings of the other tissues (seed [312, tissue index]); donor bootstrap (200 draws of the
common donors, hosts fixed, seed [312, 1000 + tissue index]).

    python score_a312.py --stage power    # leverages, bootstrap SEs and the power of H312.1; no elasticity written
    python score_a312.py --stage score    # elasticities, H312.1 and the descriptive readouts
Outputs (`../output/`): power_tissues.csv, power.csv | tissues.csv, tests.csv.
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
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
TABLE = HERE.parents[2] / "A294_floor" / "output" / "tissues.csv"  # S = z2 V per tissue (A104 fit)
N_SPLITS, N_BOOT = 20, 200
G = {}


def one(k):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    ds, comps = G["ds"], G["comps"]
    t = comps[k]
    hosts = E.split_hosts(ds, [c for c in comps if c != t], N_SPLITS, seed=[312, k])
    est = E.from_hosts(ds, t, hosts)
    donors = E.common_donors(ds, t, hosts)
    rng = np.random.default_rng([312, 1000 + k])
    boot = np.array([E.from_hosts(ds, t, hosts, donors=rng.choice(donors, len(donors)))["elasticity"] for _ in range(N_BOOT)])
    lo, hi = np.nanquantile(boot, [0.05, 0.95])
    return dict(tissue=t, elasticity=est["elasticity"], leverage=est["leverage"], n=est["n"], n_common=len(donors),
                se=float(np.nanstd(boot, ddof=1)), lo90=float(lo), hi90=float(hi))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    G.update(ds=cs.Dataset(pb, anchors=comps, name="gtex_periphery"), comps=comps)
    with Pool(a.jobs) as pool:
        rows = []
        for r in pool.imap_unordered(one, range(len(comps))):
            rows.append(r)
            print(f"{len(rows)}/{len(comps)} tissues ({time.time() - t0:.0f}s)", flush=True)
    t = pd.DataFrame(rows).set_index("tissue").sort_index()
    t["S"] = pd.read_csv(TABLE, index_col=0).S.reindex(t.index)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        t[["S", "leverage", "n", "n_common", "se"]].to_csv(OUT / "power_tissues.csv")
        rng = np.random.default_rng([312, 9])
        ls = np.log(t.S.to_numpy())
        alt = 2.0 - (ls - ls.min()) / (ls.max() - ls.min())  # 2.0 at the lowest S to 1.0 at the highest
        hits = 0
        for _ in range(2000):
            x = alt + rng.normal(0, t.se.to_numpy())
            r0 = spearmanr(x, ls).statistic
            null = np.array([spearmanr(rng.permutation(x), ls).statistic for _ in range(500)])
            hits += (1 + np.sum(null <= r0)) / 501 < 0.05
        lev = float(t.leverage.median())
        res = dict(median_leverage=lev, min_leverage=float(t.leverage.min()), median_se=float(t.se.median()),
                   power=hits / 2000)
        res["gate"] = bool(lev >= 0.15 and res["power"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(t[["S", "leverage", "n_common", "se"]].round(3).sort_values("S").to_string())
        print(pd.Series(res).to_string(), flush=True)
        return
    t.to_csv(OUT / "tissues.csv")
    rng = np.random.default_rng([312, 2])
    x, ls = t.elasticity.to_numpy(), np.log(t.S.to_numpy())
    rho = spearmanr(x, ls).statistic
    null = np.array([spearmanr(rng.permutation(x), ls).statistic for _ in range(2000)])
    p = float((1 + np.sum(null <= rho)) / 2001)
    o = t.sort_values("S")
    low, high = float(o.elasticity.iloc[:9].median()), float(o.elasticity.iloc[-9:].median())
    tests = [dict(test="H312.1", rho=float(rho), p_one_sided=p, passes=bool(rho < 0 and p < 0.05)),
             dict(test="descriptive", median_elasticity=float(t.elasticity.median()), low9_median=low, high9_median=high,
                  ratio_low_high=low / high if high > 0 else np.nan, median_leverage=float(t.leverage.median()))]
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    print(o[["S", "elasticity", "lo90", "hi90", "leverage"]].round(3).to_string())
    print(pd.DataFrame(tests).round(4).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
