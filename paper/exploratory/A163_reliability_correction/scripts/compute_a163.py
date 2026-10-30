"""A163: host-reliability correction of the coupling (registered in paper/docs/registered_experiments.md, A163).

The A162 simulator (A154 designs: n in {30, 100} x K in {4, 8} x depth in {single cell, bulk}; anchors = the first K/2
compartments at c_anchor, targets the others at c_target; no factor). Every compartment is fitted on its anchored host
(anchors on the leave-one-out host where it exists: K = 8) with a 200-replicate donor bootstrap.

  --mode dev       A162 data (`compute_a154.simulate`, seeds SeedSequence([162, ...])): every candidate reliability
                   kind of `clocksync._core.host.reliability` applied to the same point fit and bootstrap replicates.
  --mode validate  fresh seeds SeedSequence([163, design, cell, rep]); the frozen readout through the public API,
                   `Dataset.fit(attenuation=<--kind>)`. Also the A154 equal-coupling cells (c_anchor = c_target,
                   descriptive, `--eq-reps` replicates).
Raw rows: output/<mode>/raw/<design>__<cell>__<first rep>.csv (resumable); `score_a163.py` scores them.

    RUN_CAP_GB=8 paper/common/run.sh paper/exploratory/A163_reliability_correction/scripts/compute_a163.py --mode dev --reps 40 --jobs 4
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
sys.path.insert(0, str(PAPER / "exploratory" / "A154_gate_calibration" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402

import compute_a154 as C  # noqa: E402  (designs, depths, psi_for, dataset)

OUT = HERE.parent / "output"
CA, CT, C_EQ = (0.6, 0.8), (0.1, 0.3, 0.6, 0.9), (0.3, 0.6, 0.9)
N_BOOT = 200


def cells(mode: str) -> list[dict]:
    out = [dict(cell=f"ca{ca}_ct{ct}", ca=ca, ct=ct) for ca in CA for ct in CT]
    if mode == "validate":
        out += [dict(cell=f"eq_c{c}", ca=c, ct=c) for c in C_EQ]
    return out


def simulate(i: int, spec: dict, rep: int, mode: str) -> pd.DataFrame:
    if mode == "dev":  # the A162 data sets themselves
        return C.simulate(i, next(s for s in C.cell_specs() if s["cell"] == f"P162_ca{spec['ca']}_ct{spec['ct']}"), rep)
    d = C.DESIGNS[i]
    ss = np.random.SeedSequence([163, i, [c["cell"] for c in cells(mode)].index(spec["cell"]), rep])
    na = max(d["K"] // 2, 2)
    comps = {f"c{k}": Compartment(z=1.0, psi=C.psi_for(spec["ca"] if k < na else spec["ct"], 1.0)) for k in range(d["K"])}
    dp = C.depth_values.load()[d["depth"]]
    tab, _ = cs.simulate(cs.gtex(), comps, n=d["n"], window_h=C.WINDOW, jitter_h=C.JITTER, amplitude=1.0,
                         depth=C.depth_values.depth_for(dp), seed=int(ss.generate_state(1)[0]))
    return tab.rename(columns={"resection": "donor"})


def fit_all_kinds(ds: cs.Dataset, names: list, seed: int) -> dict:
    """Point fit and bootstrap of every compartment, with R of every reliability kind (as `Dataset.fit`)."""
    d = ds.data
    H = A._hosts(d, ds.anchors, names)
    point = {k: (A.on_host(d, k, H[k])["coupling"], {q: A.reliability(d, k, H[k], q) for q in A.RELIABILITY}) for k in names}
    rng = np.random.default_rng(seed)
    boots = {k: {q: [] for q in ("raw",) + A.RELIABILITY} for k in names}
    for _ in range(N_BOOT):
        b, _ = d.resample(rng)
        try:
            Hb = A._hosts(b, ds.anchors, names, max_iter=A.BOOT_MAX_ITER)
        except (ValueError, np.linalg.LinAlgError):
            continue
        for k in names:
            c = A.on_host(b, k, Hb[k])["coupling"]
            boots[k]["raw"].append(c)
            for q in A.RELIABILITY:
                boots[k][q].append(A.corrected(c, A.reliability(b, k, Hb[k], q)))
    out = {}
    for k in names:
        c, R = point[k]
        row = dict(coupling=c)
        for q in ("raw",) + A.RELIABILITY:
            x = np.asarray(boots[k][q], float)
            x = x[np.isfinite(x)]
            est = c if q == "raw" else A.corrected(c, R[q])
            row.update({f"{q}_est": est, f"{q}_R": np.nan if q == "raw" else R[q],
                        f"{q}_lo": np.quantile(x, 0.025) if len(x) else np.nan,
                        f"{q}_hi": np.quantile(x, 0.975) if len(x) else np.nan})
        out[k] = row
    return out


def run(i: int, spec: dict, rep: int, mode: str, kind: str) -> list[dict]:
    K = C.DESIGNS[i]["K"]
    na = max(K // 2, 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = C.dataset(simulate(i, spec, rep, mode), K, anchors=True)
        names = [f"c{k}" for k in range(K)] if K >= 8 else [f"c{k}" for k in range(na, K)]
        if mode == "dev":
            res = fit_all_kinds(ds, names, rep)
        else:
            t = ds.fit(gate=False, targets=names, n_boot=N_BOOT, seed=rep, attenuation=kind).table
            res = {k: {"coupling": t.loc[k, "coupling"], "raw_est": t.loc[k, "coupling"], "raw_lo": t.loc[k, "coupling_lo"],
                       "raw_hi": t.loc[k, "coupling_hi"], f"{kind}_est": t.loc[k, "coupling_corrected"],
                       f"{kind}_R": t.loc[k, "host_reliability"], f"{kind}_lo": t.loc[k, "coupling_corrected_lo"],
                       f"{kind}_hi": t.loc[k, "coupling_corrected_hi"], "boundary": bool(t.loc[k, "boundary"])} for k in names}
    base = dict(design=C.design_name(i), n=C.DESIGNS[i]["n"], K=K, depth=C.DESIGNS[i]["depth"], cell=spec["cell"],
                c_anchor=spec["ca"], c_target=spec["ct"], rep=rep)
    return [dict(base, compartment=k, role="anchor" if int(k[1:]) < na else "target",
                 true_c=spec["ca"] if int(k[1:]) < na else spec["ct"], **r) for k, r in res.items()]


def block(i: int, cell: str, reps: list, mode: str, kind: str, out: str) -> str:
    limit_threads(1)
    f = Path(out) / mode / "raw" / f"{C.design_name(i)}__{cell}__{reps[0]:03d}.csv"
    if f.exists():
        return str(f)
    spec = next(c for c in cells(mode) if c["cell"] == cell)
    rows = [r for rep in reps for r in run(i, spec, rep, mode, kind)]
    f.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(f.with_suffix(".tmp"), index=False)
    f.with_suffix(".tmp").rename(f)
    return str(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", choices=("dev", "validate"), required=True)
    ap.add_argument("--kind", default="df", help="validate: the frozen reliability kind")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--eq-reps", type=int, default=100)
    ap.add_argument("--block", type=int, default=5)
    ap.add_argument("--designs", default=None)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    ds = [int(x) for x in args.designs.split(",")] if args.designs else range(len(C.DESIGNS))
    tasks = [(i, c["cell"], list(range(r, min(r + args.block, n))))
             for i in ds for c in cells(args.mode)
             for n in [args.eq_reps if c["cell"].startswith("eq") else args.reps] for r in range(0, n, args.block)]
    tasks.sort(key=lambda t: -C.DESIGNS[t[0]]["n"] * C.DESIGNS[t[0]]["K"])  # longest first
    pending = [t for t in tasks if not (args.out / args.mode / "raw" / f"{C.design_name(t[0])}__{t[1]}__{t[2][0]:03d}.csv").exists()]
    print(f"{len(tasks)} tasks ({len(pending)} to run), {args.jobs} jobs, mode {args.mode}", flush=True)
    t0, done = time.time(), 0
    with cf.ProcessPoolExecutor(args.jobs) as ex:
        futs = [ex.submit(block, i, c, r, args.mode, args.kind, str(args.out)) for i, c, r in pending]
        for f in cf.as_completed(futs):
            f.result()
            done += 1
            if done % 10 == 0 or done == len(futs):
                print(f"{time.strftime('%H:%M:%S')} {done}/{len(futs)} tasks ({time.time() - t0:.0f} s)", flush=True)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
