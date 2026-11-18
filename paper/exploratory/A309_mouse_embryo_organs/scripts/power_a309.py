"""A309 adult power (before any embryo statistic): timed adult mouse stand-ins at the design of each readout.

Registered in paper/docs/registered_experiments.md (A309, "Power first"). Per readout, every embryo is stood for by a
distinct adult sample of the organ (Wolff 2023: adrenal, lung, kidney, heart, all ages; Manella 2021: liver, AlbCre ad
libitum), assigned to the readout's batch structure at random; clock-gene counts thinned gene by gene to the cell type's
median count; clock at 1, 0.75, 0.5 (`simulate.scale_clock`) and shuffled; C with 500 permutations within batch; 200
subsamples. Only the embryo design (counts per batch, median clock-gene counts) enters; no embryo statistic.

    python power_a309.py [--n-sub 200] [--jobs 4]
Outputs (../output/): power_design.csv, power_subsamples.csv, power_passrates.csv.
"""
import argparse
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

import lib_a309 as L
from clocksync.config import limit_threads
from clocksync.simulate import scale_clock
from common import data
from common import development as dev

SCALES = (1.0, 0.75, 0.5)
N_PERM = 500
S = {}
ADULT = {"adrenal": ("mouse_ageing_wolff", "adrenal"), "lung": ("mouse_ageing_wolff", "lung"),
         "kidney": ("mouse_ageing_wolff", "kidney"), "heart": ("mouse_ageing_wolff", "heart"),
         "liver": ("manella2021", "Liver")}


def pools(rng):
    out = {}
    for r, (name, comp) in ADULT.items():
        pb, don = data.load(name)
        x = pb[pb.compartment == comp]
        if name == "manella2021":
            d = don.reindex(x.donor)
            x = x[((d.genotype == "AlbCre") & (d.feeding == "AL")).to_numpy()]
        x = x[["donor", "compartment", "depth", *L.GENES]].reset_index(drop=True)
        out[r] = dict(n=len(x), depth=x.depth.to_numpy(float), med=np.median(x[L.GENES].to_numpy(float), 0),
                      x={a: scale_clock(x, L.WF, a, depth=x.depth, rng=rng)[L.GENES].to_numpy(float) for a in SCALES})
    return out


def one(i):
    rng, out = np.random.default_rng([L.SEED, 1, i]), []
    for j, (r, f) in enumerate(S["design"].items()):
        P = S["pools"][r]
        pick = rng.choice(P["n"], len(f), replace=False)
        p = np.minimum(1.0, S["target"][r] / np.maximum(P["med"], 1e-9))
        base = P["x"][1.0][pick]
        frac = base / base.sum(1, keepdims=True)
        null = np.column_stack([rng.permutation(frac[:, k]) for k in range(frac.shape[1])]) * base.sum(1, keepdims=True)
        t0 = f.reset_index(drop=True)[["donor", "compartment", "block"]].copy()
        t0["compartment"] = "standin"
        t0["block"] = rng.permutation(t0.block.to_numpy())
        for setting, x in [*((f"a{a:g}", P["x"][a][pick]) for a in SCALES), ("shuffled", null)]:
            t = t0.copy()
            t[L.GENES] = rng.binomial(np.round(x).astype(np.int64), p).astype(float)
            t["depth"] = P["depth"][pick]  # thinning by a gene factor is a constant log offset per gene
            tab = dev.adjust(t, L.GENES, L.cov)[0]
            c, pv, _ = dev.coherence(tab, L.WF, N_PERM, np.random.default_rng([L.SEED, 2, i, j]), "block")
            out.append(dict(readout=r, setting=setting, sub=i, C=c, p=pv))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    limit_threads(1)
    L.OUT.mkdir(parents=True, exist_ok=True)
    rows = L.late_rows()
    design = {r: rows[rows.readout == r] for r in L.READOUTS}
    target = {r: np.median(f[L.GENES].to_numpy(float), 0) for r, f in design.items()}
    S.update(design=design, target=target, pools=pools(np.random.default_rng([L.SEED, 0])))
    drow = pd.DataFrame([dict(readout=r, embryos=len(f), batches=f.block.nunique(),
                              per_batch=";".join(map(str, f.block.value_counts().sort_index())),
                              median_clock_count=float(target[r].sum()), adult_pool=S["pools"][r]["n"],
                              adult_median_clock_count=float(S["pools"][r]["med"].sum()),
                              genes_thinned=int((target[r] < S["pools"][r]["med"]).sum())) for r, f in design.items()])
    drow.to_csv(L.OUT / "power_design.csv", index=False)
    pd.set_option("display.width", 250)
    print(drow.to_string(index=False), flush=True)
    res = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, range(a.n_sub))):
            res.extend(r)
            if (k + 1) % 20 == 0:
                print(f"{k + 1}/{a.n_sub} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res)
    res.to_csv(L.OUT / "power_subsamples.csv", index=False)
    rates = res.groupby(["readout", "setting"], sort=False).agg(pass_rate=("p", lambda s: float((s < .05).mean())),
                                                                 median_C=("C", "median"), n=("p", "size")).reset_index()
    rates.to_csv(L.OUT / "power_passrates.csv", index=False)
    print(rates.round(3).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
