"""A154: gate calibration under shared non-clock structure (registered in paper/docs/registered_experiments.md, A154).

Simulator `clocksync.simulate` (current code, GTEx waveform of the 14 clock genes, gate n_perm 1,000). Designs: donors n in
{30, 100} x compartments K in {4, 8} x depth in {single cell, bulk} (depth values: `A155.../scripts/depth_values.py`).
Arms and readouts (one task kind each; `--kinds`):
  N0     no clock, no factor: gate pass rate (V154.1).
  N1     no clock + shared non-clock donor factor (rank 1 or 3, loading on 10 or 30 % of 2,000 genes, size 0.25-2 x the A140 Q1
         calibration): gate pass rate (descriptive map).
  P      one coupled clock, c in {0.3, 0.6, 0.9}, without the factor and with it (scale 1, rank 1, 30 %): gate pass rate.
  GNULL  global null: one cohort of 2n donors split in two arms of n, clock c = 0.6 in both, `Dataset.contrast`: p of every non-anchor target (V154.2).
  COVER  arm P without factor: `Dataset.fit(gate=False)` 95 % bootstrap interval of c of every non-anchor target
         (V154.3) and the simulation lower bound (`Dataset.lower_bound`) when the fit is at the boundary.
The first K / 2 compartments are anchors (targets of fit / contrast: the others); the gate uses all compartments (free host).
Replicates 200 per cell, seed 154 (SeedSequence([154, design, cell, replicate])). Outputs under `--out`: raw/<kind>/<cell>_<block>.csv
(resumable), then `--score` writes the tables of `score_a154.py`.

    RUN_CAP_GB=8 paper/common/run.sh paper/exploratory/A154_gate_calibration/scripts/compute_a154.py --kinds N0,GNULL,COVER --jobs 6
    options: --reps N (200), --block B (replicates per task), --designs i,j (indices 0-7), --out DIR
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
A155 = PAPER / "exploratory" / "A155_anchor_self_inclusion" / "scripts"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(A155))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import Compartment, thin_clock  # noqa: E402

import depth_values  # noqa: E402
from compute_a155 import psi_for  # noqa: E402

OUT = HERE.parent / "output"
SEED, WINDOW, JITTER = 154, (8.0, 18.0), 1.0
CLOCK = list(cs.gtex().genes)
G_PANEL = 2000  # genes over which the factor loads (A140 panel size); the clock genes are the first 14
LAM = np.array([411.4, 243.6, 195.2])  # A140 Q1 calibration (gate-passing tumour immune genome-wide eigenvalues)
LAM_FILE = PAPER / "exploratory" / "A140_dip_cyclops" / "output" / "q1" / "tumour_immune_genome" / "calibration_lambda.csv"
C_GRID, C_NULL = (0.3, 0.6, 0.9), 0.6
DESIGNS = [dict(n=n, K=K, depth=d) for n in (30, 100) for K in (4, 8) for d in ("single_cell", "bulk")]
KINDS = ("N0", "N1", "P", "GNULL", "COVER")


def lam() -> np.ndarray:
    return pd.read_csv(LAM_FILE).lam.to_numpy() if LAM_FILE.exists() else LAM


def design_name(i: int) -> str:
    d = DESIGNS[i]
    return f"n{d['n']}_K{d['K']}_{d['depth']}"


def cell_specs() -> list[dict]:
    """Every cell of an arm: kind, id and the simulation parameters (factor: rank, load, scale; c)."""
    out = [dict(kind="N0", cell="N0", c=None, fac=None)]
    for r in (1, 3):
        for ld in (0.1, 0.3):
            for s in (0.25, 0.5, 1.0, 2.0):
                out.append(dict(kind="N1", cell=f"N1_r{r}_l{int(ld * 100)}_s{s}", c=None, fac=(r, ld, s)))
    for c in C_GRID:
        out.append(dict(kind="P", cell=f"P_c{c}_nofac", c=c, fac=None))
        out.append(dict(kind="P", cell=f"P_c{c}_fac", c=c, fac=(1, 0.3, 1.0)))
        out.append(dict(kind="COVER", cell=f"COVER_c{c}", c=c, fac=None))
    out.append(dict(kind="GNULL", cell="GNULL", c=C_NULL, fac=None))
    return out


def loadings(rank: int, load: float, rng: np.random.Generator) -> np.ndarray:
    """G_PANEL x rank loadings: each axis loads (normal, unit norm) on a random `load` share of the genes, clock genes eligible."""
    L = np.zeros((G_PANEL, rank))
    for k in range(rank):
        idx = rng.choice(G_PANEL, int(round(load * G_PANEL)), replace=False)
        v = rng.normal(size=len(idx))
        L[idx, k] = v / np.linalg.norm(v)
    return L


def simulate(i: int, spec: dict, rep: int, arm: int = 0, n: int | None = None) -> pd.DataFrame:
    """One cohort: long table (donor, compartment, depth, 14 clock genes) at design `i`, cell `spec`."""
    d = DESIGNS[i]
    key = [SEED, i, [s["cell"] for s in cell_specs()].index(spec["cell"]), rep]
    if spec["kind"] == "COVER":  # same data as the P arm without factor
        key[2] = [s["cell"] for s in cell_specs()].index(f"P_c{spec['c']}_nofac")
    ss = np.random.SeedSequence(key + [arm])
    seed = int(ss.generate_state(1)[0])
    rng = np.random.default_rng(ss.spawn(1)[0])
    clock = spec["c"] is not None
    psi = psi_for(spec["c"] if clock else 0.6, 1.0)
    comps = {f"c{k}": Compartment(z=1.0, psi=psi) for k in range(d["K"])}
    dp = depth_values.load()[d["depth"]]
    n = d["n"] if n is None else n
    tab, _ = cs.simulate(cs.gtex(), comps, n=n, window_h=WINDOW, jitter_h=JITTER, amplitude=1.0 if clock else 0.0,
                         depth=depth_values.depth_for(dp), seed=seed)
    if spec["fac"] is not None:
        rank, load, scale = spec["fac"]
        L = loadings(rank, load, rng)[:len(CLOCK)]
        S = rng.normal(size=(d["n"], rank)) * np.sqrt(lam()[:rank]) * scale
        pid = tab.resection.str[1:].astype(int).to_numpy()
        tab = tab.copy()
        tab[CLOCK] = thin_clock(tab[CLOCK].to_numpy(float), L, S[pid], rng)
    return tab.rename(columns={"resection": "donor"})


def dataset(tab: pd.DataFrame, K: int, anchors: bool) -> cs.Dataset:
    names = [f"c{k}" for k in range(K)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cs.Dataset(tab[["donor", "compartment", "depth"] + CLOCK], anchors=names[:max(K // 2, 2)] if anchors else None)


# ---- one replicate --------------------------------------------------------------------------------------------------

def run(i: int, spec: dict, rep: int) -> list[dict]:
    d, kind = DESIGNS[i], spec["kind"]
    base = dict(design=design_name(i), n=d["n"], K=d["K"], depth=d["depth"], cell=spec["cell"], kind=kind, rep=rep)
    t0 = time.time()
    if kind in ("N0", "N1", "P"):
        ds = dataset(simulate(i, spec, rep), d["K"], anchors=False)
        g = ds.gate(n_perm=1000, seed=rep)
        return [dict(base, passes=bool(g.passes), p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                     clock_share=g.clock_share, share_null=g.share_null_mean, s=time.time() - t0)]
    K, targets = d["K"], [f"c{k}" for k in range(d["K"] // 2 if d["K"] > 4 else 2, d["K"])]
    if kind == "GNULL":
        # one cohort of 2n donors (same compartment offsets in both arms), the first n donors are arm 1: a global null
        tab = simulate(i, spec, rep, 0, n=2 * d["n"])
        arm1 = set(sorted(tab.donor.unique())[:d["n"]])
        ds = dataset(tab, K, anchors=True)
        lab = pd.Series(["control" if x in arm1 else "case" for x in sorted(tab.donor.unique())],
                        index=sorted(tab.donor.unique()))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = ds.contrast(lab, levels=("control", "case"), targets=targets, n_perm=1000, n_boot=0, seed=rep).table
        return [dict(base, compartment=k, p=r.loc[k, "p"], p_star=r.loc[k, "p_star"], d=r.loc[k, "d"],
                     status=r.loc[k, "status"], s=time.time() - t0) for k in targets]
    ds = dataset(simulate(i, spec, rep), K, anchors=True)  # COVER
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = ds.fit(gate=False, targets=targets, n_boot=200, seed=rep)
    rows = []
    for k in targets:
        r = f.table.loc[k]
        row = dict(base, compartment=k, true_c=spec["c"], coupling=r.coupling, lo=r.coupling_lo, hi=r.coupling_hi,
                   boundary=bool(r.boundary), c_lo=np.nan)
        if row["boundary"]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                row["c_lo"] = ds.lower_bound(k, n_sim=200, seed=rep)["c_lo"]
        rows.append(row)
    for r in rows:
        r["s"] = time.time() - t0
    return rows


def block(i: int, cell: str, reps: list, out: str) -> str:
    limit_threads(1)
    f = Path(out) / "raw" / cell.split("_")[0] / f"{design_name(i)}__{cell}__{reps[0]:03d}.csv"
    if f.exists():
        return str(f)
    spec = next(s for s in cell_specs() if s["cell"] == cell)
    rows = [r for rep in reps for r in run(i, spec, rep)]
    f.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(f.with_suffix(".tmp"), index=False)
    f.with_suffix(".tmp").rename(f)
    return str(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kinds", default="N0,GNULL,COVER")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--first", type=int, default=0, help="first replicate (for a timing run)")
    ap.add_argument("--block", type=int, default=5)
    ap.add_argument("--designs", default=None)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    kinds = args.kinds.split(",")
    ds = [int(x) for x in args.designs.split(",")] if args.designs else range(len(DESIGNS))
    cost = lambda i, s: DESIGNS[i]["n"] * DESIGNS[i]["K"] * (2 if s["kind"] == "GNULL" else 1)  # tasks run in kind order, longest first within a kind
    tasks = [(i, s["cell"], list(range(r, min(r + args.block, args.first + args.reps))))
             for i in ds for s in cell_specs() if s["kind"] in kinds
             for r in range(args.first, args.first + args.reps, args.block)]
    prio = {k: j for j, k in enumerate(("N0", "GNULL", "COVER", "N1", "P"))}
    tasks.sort(key=lambda t: (prio[next(s for s in cell_specs() if s["cell"] == t[1])["kind"]],
                              -cost(t[0], next(s for s in cell_specs() if s["cell"] == t[1]))))
    pending = [t for t in tasks if not (args.out / "raw" / t[1].split("_")[0] /
                                       f"{design_name(t[0])}__{t[1]}__{t[2][0]:03d}.csv").exists()]
    print(f"{len(tasks)} tasks ({len(pending)} to run), {args.jobs} jobs, kinds {kinds}", flush=True)
    t0, done = time.time(), 0
    with cf.ProcessPoolExecutor(args.jobs) as ex:
        futs = [ex.submit(block, i, c, r, str(args.out)) for i, c, r in pending]
        for f in cf.as_completed(futs):
            f.result()
            done += 1
            if done % 10 == 0 or done == len(futs):
                print(f"{time.strftime('%H:%M:%S')} {done}/{len(futs)} tasks ({time.time() - t0:.0f} s)", flush=True)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
