"""A181 V181.1 (with the registered design note): GTEx depth control of the fetal clock gate (H107b.1 design).

Donors: GTEx deaths at 08-18 h (any Hardy class) with >= 2 of {Brain - Frontal Cortex (BA9), Heart - Left Ventricle,
Kidney - Cortex, Lung}. Each subsample draws 23 donors and assigns the 23 H107b.1 fetal organ sets at random (a donor lacking
an organ of its assigned set is redrawn; organs outside the set are masked). Clock-gene counts are thinned binomially,
p = min(1, target_clock_depth / observed_clock_depth) with target = Cao median whole-organ depth (build_cao.md) x
clock-gene count / library size, i.e. p = min(1, Cao depth / library size); the depth is scaled by p (clock share kept).
Gate: clocksync Dataset.gate(n_perm=1000), GTEx global waveform; pass = p_pairing and p_acrophase < 0.05. The unthinned
gate on the same subsamples is the reference. Controlled access: time of death is used in memory to select donors only;
only aggregates are written.

    python v181_1_depth_control.py [--n-sub 200] [--jobs 4] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data, gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

OUT = HERE.parents[1] / "output"
ORGAN = {"Cerebrum": "Brain - Frontal Cortex (BA9)", "Heart": "Heart - Left Ventricle", "Kidney": "Kidney - Cortex",
         "Lung": "Lung"}
CAO_DEPTH = {"Cerebrum": 1.4e8, "Heart": 5.1e6, "Kidney": 9.2e6, "Lung": 1.1e7}  # build_cao.md, median whole-organ
N_DON, N_PERM, SEED = 23, 1000, 181
G = None  # per-process globals: (table by donor, pool, fetal sets, clock genes)


def setup():
    global G
    pb, _ = gtex.load(list(ORGAN.values()), hardy=(0, 1, 2, 3, 4), depth="library_size")  # tod needed: in memory only
    _, d = gtex.load(list(ORGAN.values()), hardy=(0, 1, 2, 3, 4), depth="library_size")
    win = set(d.index[in_window(d.tod.to_numpy())])
    pb = pb[pb.donor.isin(win)].reset_index(drop=True)
    del d
    genes = list(cs.gtex().genes)
    n = pb.groupby("donor").compartment.nunique()
    pool = sorted(n.index[n >= 2])
    cao, _ = data.load("fetal_cao", level="organ")
    el = S.eligible(cao, 1000, 8, compartments=list(ORGAN))
    sets = [tuple(sorted(g.compartment)) for _, g in el.groupby("donor")]
    assert len(sets) == 23
    tabs = {k: g.set_index("compartment") for k, g in pb[pb.donor.isin(pool)].groupby("donor")}
    cao_clock = el[genes].sum(1).groupby(el.compartment).median() if set(genes) <= set(el) else None
    G = (tabs, pool, sets, genes, cao_clock)


def draw(rng):
    tabs, pool, sets, genes, _ = G
    for _try in range(1000):
        order = rng.permutation(len(sets))
        used, assign = set(), []
        for j in order:
            ok = [d for d in pool if d not in used and all(ORGAN[o] in tabs[d].index for o in sets[j])]
            if not ok:
                break
            d = ok[rng.integers(len(ok))]
            used.add(d)
            assign.append((d, sets[j]))
        else:
            return assign, _try
    raise RuntimeError("no valid draw")


def one(i):
    limit_threads(1)
    tabs, pool, sets, genes, _ = G
    rng = np.random.default_rng([SEED, i])
    assign, restarts = draw(rng)
    inv = {v: k for k, v in ORGAN.items()}
    rows = []
    for d, st in assign:
        for o in st:
            r = tabs[d].loc[ORGAN[o]]
            rows.append(dict(donor=d, compartment=o, depth=float(r.depth), **{g: float(r[g]) for g in genes}))
    raw = pd.DataFrame(rows)
    th = raw.copy()
    obs = th[genes].sum(1).to_numpy()
    p = np.minimum(1.0, th.compartment.map(CAO_DEPTH).to_numpy() / th.depth.to_numpy())
    th[genes] = rng.binomial(th[genes].to_numpy(int), p[:, None])
    th["depth"] = th.depth * p
    out = []
    for tag, tab in (("thinned", th), ("unthinned", raw)):
        g = cs.Dataset(tab, anchors=None, name="a181").gate(n_perm=N_PERM, seed=SEED)
        out.append(dict(sub=i, kind=tag, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share,
                        n_pairs=g.n_pairs, passed=bool(g.p_pairing < .05 and g.p_acrophase < .05), restarts=restarts))
    depth = pd.DataFrame(dict(compartment=th.compartment, clock_raw=obs, clock_thinned=th[genes].sum(1).to_numpy(),
                              p=p)).assign(sub=i)
    return pd.DataFrame(out), depth


def ci(k, n):
    r = binomtest(int(k), int(n)).proportion_ci(0.95, method="exact")
    return r.low, r.high


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    limit_threads(1)
    setup()
    n = 4 if a.check else a.n_sub
    print(f"pool donors {len(G[1])}; {len(G[2])} fetal organ sets; {n} subsamples", flush=True)
    res, dep = [], []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (r, d) in enumerate(pool.imap(one, range(n))):
            res.append(r), dep.append(d)
            if (k + 1) % 10 == 0:
                print(f"{k + 1}/{n} {time.time() - t0:.0f}s", flush=True)
    res, dep = pd.concat(res), pd.concat(dep)
    sfx = "_check" if a.check else ""
    res.to_csv(OUT / f"v181_1_subsamples{sfx}.csv", index=False)
    rows = []
    for kind, r in res.groupby("kind"):
        k, m = int(r.passed.sum()), len(r)
        lo, hi = ci(k, m)
        rows.append(dict(kind=kind, n_subsamples=m, n_pass=k, pass_rate=k / m, ci_lo=lo, ci_hi=hi,
                         median_p_pairing=r.p_pairing.median(), median_p_acrophase=r.p_acrophase.median(),
                         frac_p_pairing_lt05=(r.p_pairing < .05).mean(), frac_p_acro_lt05=(r.p_acrophase < .05).mean(),
                         restarts_total=int(r.restarts.sum() / 2)))
    pd.DataFrame(rows).to_csv(OUT / f"v181_1_passrates{sfx}.csv", index=False)
    cao = G[4]
    dd = dep.groupby("compartment").agg(median_clock_depth_unthinned=("clock_raw", "median"),
                                        median_clock_depth_thinned=("clock_thinned", "median"),
                                        median_p=("p", "median"), n_rows=("p", "size"))
    dd["cao_median_clock_depth"] = cao.reindex(dd.index) if cao is not None else np.nan
    dd["cao_pseudobulk_depth_target"] = pd.Series(CAO_DEPTH)
    dd.to_csv(OUT / f"v181_1_depths{sfx}.csv")
    print(pd.DataFrame(rows).round(4).to_string(index=False))
    print(dd.round(4).to_string())
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
