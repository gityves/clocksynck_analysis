"""A271 Part B (registered in `paper/docs/registered_experiments.md`): is "stronger clocks are more coupled" (A104, GTEx
fast deaths, Spearman rho(own amplitude, coupling) over tissues) specific to the clock genes?

Rows, blocks and fit as A104 (`results/figure_2/scripts/backbone/compute_oscillator_gtex.py`: tissues with >= 40 donors,
periphery and brain blocks, every tissue an anchor, leave-one-out host). For the clock and each of 200 expression-matched
and 200 responsive pseudo-clock sets (`clocksync.paired.pseudo_waveforms`, matching variable = mean log fraction over all
rows, genes with mean count >= 1): rho per block (point estimate, no bootstrap). Aggregates only (GTEx controlled access).

    python score_a271_a104.py [--n-sets 200] [--jobs 6] [--check]
Outputs (`../output/`): a104_sets.csv (one row per set and block), a104_tests.csv.
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

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 271
WF = cs.gtex()
META = ["donor", "compartment", "depth"]
S = {}


def load():
    import anndata as ad
    from common import gtex
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    return pb


def task(i):
    warnings.filterwarnings("ignore")
    kind, w = S["sets"][i]
    pb, rows = S["pb"], []
    for block in ("periphery", "brain"):
        comps = sorted(c for c in pb.compartment.unique() if c.startswith("Brain") == (block == "brain"))
        out = dict(set=i, kind=kind, block=block, genes=";".join(w.genes))
        try:
            ds = cs.Dataset(pb.loc[pb.compartment.isin(comps), META + list(w.genes)], anchors=comps, waveform=w, name=f"{block}{i}")
            pt = A104.per_tissue(ds)
            out.update(tissues=int(pt[["own_amplitude", "coupling"]].notna().all(axis=1).sum()), rho=A104.rho(pt),
                       median_coupling=float(pt.coupling.median()))
        except Exception as e:  # noqa: BLE001
            out["error"] = f"{type(e).__name__}: {e}"[:200]
        rows.append(out)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sets", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, n = time.time(), (3 if a.check else a.n_sets)
    pb = load()  # before limit_threads: pyarrow with OMP_NUM_THREADS=1 peaks far higher on wide tables
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    print(f"{pb.donor.nunique()} donors, {pb.compartment.nunique()} tissues, {len(lm)} genes for matching {time.time() - t0:.0f}s",
          flush=True)
    sets = [("clock", WF)]
    for j, kind in enumerate(("matched", "responsive"), 1):
        sets += [(kind, w) for w in pseudo_waveforms(lm, n, np.random.default_rng([SEED, j]), WF, kind)]
    used = sorted({g for _, w in sets for g in w.genes})
    S.update(sets=sets, pb=pb[META + used].copy())
    del pb
    rows = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, range(len(sets)))):
            rows += r
            if (k + 1) % 25 == 0 or k + 1 == len(sets):
                print(f"{k + 1}/{len(sets)} {time.time() - t0:.0f}s", flush=True)
    R = pd.DataFrame(rows).sort_values(["block", "set"])
    sfx = "_check" if a.check else ""
    R.to_csv(OUT / f"a104_sets{sfx}.csv", index=False)
    tests = []
    for block, g in R.groupby("block"):
        c = g[g.kind == "clock"].iloc[0]
        out = dict(block=block, rho_clock=c.rho, tissues=c.tissues)
        spec = True
        for kind in ("matched", "responsive"):
            r = g[g.kind == kind].rho.dropna()
            out[f"{kind}_valid"] = len(r)
            out[f"{kind}_rho_median"] = float(r.median())
            out[f"{kind}_rho_q95"] = float(r.quantile(0.95))
            out[f"{kind}_p"] = float((1 + (r >= c.rho).sum()) / (1 + len(r)))
            spec &= bool(c.rho > out[f"{kind}_rho_q95"])
        out["clock_specific"] = spec
        tests.append(out)
    T = pd.DataFrame(tests)
    T.to_csv(OUT / f"a104_tests{sfx}.csv", index=False)
    print(T.T.to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
