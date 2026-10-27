"""A145: what is specific to clocksync against two-step orderer pipelines? Fits and readouts per (arm, replicate, method).

Registered in paper/docs/registered_experiments.md (A145). Simulated design (D145.1, H145.1-3): `lib.py`, arms F, W,
S_small, S_large, N, 30 replicates (seed 145). Real known answer (D145.2): A129, Manella 2021 hepatocyte Bmal1 knockout (BLKO)
vs Alb-Cre, ad libitum, liver read against the WAT + lung host, per genotype.

Tasks (cached one pickle each under <out>/runs/, a rerun resumes): the orderers CHIRAL, CYCLOPS 2.0 (clock-gene seed) and
CoFE on the donor x (gene, anchor) log2 matrix, and clocksync's Dataset.fit. Then, in the main process, step 2 (regression of
clocksync's clock-plane state of each target on the inferred phase) and the timed oracle (the same on the true times).

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A145_two_step/scripts/compute_a145.py --jobs 4
    options: --reps N (default 30), --out DIR (default exploratory/A145_two_step/output), --arms ..., --no-real

Workers: --jobs processes, 1 thread each (CYCLOPS: 1 Julia worker; CoFE: 1 core). Outputs: results_long.csv (simulation),
runtime.csv, real_long.csv (D145.2), real_clocksync.csv.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import importlib.util
import multiprocessing as mp
import shutil
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import lib as A  # noqa: E402

COST = {"cyclops": 0, "cofe": 1, "chiral": 2, "clocksync": 3}
REAL = ("AlbCre", "BLKO")


# ---- tasks ----------------------------------------------------------------------------------------------------------

def sim_task(arm: str, rep: int, method: str, runs: str, threads: int = 1) -> dict:
    from clocksync.config import limit_threads

    limit_threads(1)
    runs = Path(runs)
    f = runs / f"{arm}_{rep:03d}_{method}.pkl"
    if f.exists():
        r = pd.read_pickle(f)
        return dict(key=(arm, rep, method), runtime_s=r["runtime_s"], failed=r["failed"], cached=True)
    tab, _, _ = A.cohort(arm, rep)
    t0 = time.time()
    if method == "clocksync":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = A.clocksync_fit(tab, seed=rep)
        r.update(runtime_s=time.time() - t0, failed=False, error="")
    else:
        tmp = runs / "tmp" / f.stem
        ph, err, dt = A.run_orderer(method, A.orderer_input(tab), tmp, threads=threads, seed=rep)
        shutil.rmtree(tmp, ignore_errors=True)
        r = dict(phase=ph, runtime_s=dt, failed=ph is None, error=err)
    pd.to_pickle(r, f)
    return dict(key=(arm, rep, method), runtime_s=r["runtime_s"], failed=r["failed"], cached=False)


def real_inputs():
    """Manella 2021, AL mice per genotype: pseudobulk, donor table, waveform without the genes failing QC, genes."""
    spec = importlib.util.spec_from_file_location("a111c", A.PAPER / "exploratory/A111_feeding_controls/scripts/_common.py")
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
    from common import data

    pb, donors = data.load("manella2021")
    qc = C.gene_qc(pb, "manella2021")
    wf = C.waveform_without(sorted(qc[~qc.kept].gene.unique()))
    return pb, donors, wf


def real_task(geno: str, method: str, runs: str, threads: int = 1) -> dict:
    from clocksync.config import limit_threads

    limit_threads(1)
    runs = Path(runs)
    f = runs / f"real_{geno}_{method}.pkl"
    if f.exists():
        r = pd.read_pickle(f)
        return dict(key=("real", geno, method), runtime_s=r["runtime_s"], failed=r["failed"], cached=True)
    import clocksync as cs

    pb, donors, wf = real_inputs()
    ids = donors.index[(donors.genotype == geno) & (donors.feeding == "AL")]
    x = pb[pb.donor.isin(ids)]
    t0 = time.time()
    if method == "clocksync":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds = cs.Dataset(x, anchors=["WAT", "Lung"], donors=donors, waveform=wf, name=geno)
            fit = ds.fit(n_boot=0, gate=False, targets=["Liver"])
            g = ds.gate(n_perm=1000, targets=["Liver"], hosts=["WAT", "Lung"])
        r = dict(table=fit.table.loc[["Liver"]], gate=dict(passes=bool(g.passes), p_pairing=g.p_pairing, p_acrophase=g.p_acrophase),
                 runtime_s=time.time() - t0, failed=False, error="")
    else:
        tmp = runs / "tmp" / f.stem
        ph, err, dt = A.run_orderer(method, A.orderer_input(x, ["WAT", "Lung"], list(wf.genes)), tmp, threads=threads, seed=0)
        shutil.rmtree(tmp, ignore_errors=True)
        r = dict(phase=ph, runtime_s=dt, failed=ph is None, error=err)
    pd.to_pickle(r, f)
    return dict(key=("real", geno, method), runtime_s=r["runtime_s"], failed=r["failed"], cached=False)


# ---- assembly -----------------------------------------------------------------------------------------------------------

def sim_rows(arm: str, rep: int, runs: Path) -> list[dict]:
    tab, d, t = A.cohort(arm, rep)
    dtruth = d.assign(offset_h=d.offset_h)
    rows = []
    base = lambda m, k: dict(arm=arm, rep=rep, method=m, target=k, c_true=d.c[k], offset_true_h=d.offset_h[k], amp_true=d.amp[k])
    # clocksync
    r = pd.read_pickle(runs / f"{arm}_{rep:03d}_clocksync.pkl")
    rel = A.clocksync_relative(r["table"])
    for k in A.TARGETS:
        rows.append(dict(**base("clocksync", k), failed=False, gate_pass=r["gate"]["passes"], c_hat=rel.coupling[k],
                         offset_hat_h=rel.offset_h[k], log_amp_hat=rel.log_amp_rel[k]))
    # two-step: orderers and the timed oracle
    phis = {"oracle": pd.Series(t / A.HOURS, index=sorted(tab.donor.unique()))}
    for m in A.METHODS:
        ph = pd.read_pickle(runs / f"{arm}_{rep:03d}_{m}.pkl")["phase"]
        phis[m] = ph
    for m, ph in phis.items():
        if ph is None:
            res = None
        else:
            res = A.two_step(tab, ph)
        for k in A.TARGETS:
            e = dict(**base(m, k), failed=res is None, gate_pass=np.nan)
            if res is not None:
                e.update(c_hat=res.coupling[k], offset_hat_h=res.offset_h[k], log_amp_hat=res.log_amp_rel[k])
            rows.append(e)
    return rows


def real_rows(runs: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    pb, donors, wf = real_inputs()
    genes = list(wf.genes)
    out, cl = [], []
    for geno in REAL:
        ids = donors.index[(donors.genotype == geno) & (donors.feeding == "AL")]
        x = pb[pb.donor.isin(ids)]
        r = pd.read_pickle(runs / f"real_{geno}_clocksync.pkl")
        cl.append(dict(genotype=geno, coupling=r["table"].coupling["Liver"], offset_h=r["table"].offset_h["Liver"], **r["gate"]))
        phis = {"oracle": (donors.zt.reindex(ids).astype(float) / 24 * 2 * np.pi)}
        for m in A.METHODS:
            phis[m] = pd.read_pickle(runs / f"real_{geno}_{m}.pkl")["phase"]
        for m, ph in phis.items():
            e = dict(genotype=geno, method=m, failed=ph is None)
            if ph is not None:
                res = A.two_step(x, ph, targets=["Liver"], anchors=["WAT", "Lung"], genes=genes)
                e.update(c_hat=res.coupling["Liver"], offset_h=res.offset_h["Liver"], log_amp_rel=res.log_amp_rel["Liver"])
            out.append(e)
    return pd.DataFrame(out), pd.DataFrame(cl)


# ---- main ---------------------------------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", default=str(A.OUT))
    ap.add_argument("--arms", nargs="*", default=list(A.ARMS))
    ap.add_argument("--no-real", action="store_true")
    ap.add_argument("--no-sim", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.out)
    runs = out / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    todo = []
    if not args.no_sim:
        for arm in args.arms:
            for rep in range(args.reps):
                for m in ("clocksync",) + A.METHODS:
                    if not (runs / f"{arm}_{rep:03d}_{m}.pkl").exists():
                        todo.append(("sim", arm, rep, m))
    if not args.no_real:
        for g in REAL:
            for m in ("clocksync",) + A.METHODS:
                if not (runs / f"real_{g}_{m}.pkl").exists():
                    todo.append(("real", g, 0, m))
    todo.sort(key=lambda k: (COST[k[3]], k[0] != "real", k[1], k[2]))
    print(f"workers = {args.jobs} x 1 thread; {len(todo)} tasks to run (arms {args.arms}, {args.reps} replicates)", flush=True)
    t0, done, last = time.time(), 0, time.time()
    if todo:
        with cf.ProcessPoolExecutor(args.jobs, mp_context=mp.get_context("spawn")) as ex:
            pend = {}
            for kind, a, rep, m in todo:
                fu = ex.submit(sim_task, a, rep, m, str(runs)) if kind == "sim" else ex.submit(real_task, a, m, str(runs))
                pend[fu] = (kind, a, rep, m)
            while pend:
                fin, _ = cf.wait(pend, timeout=600, return_when=cf.FIRST_COMPLETED)
                for fu in fin:
                    r = fu.result()
                    pend.pop(fu)
                    done += 1
                    print(f"[{time.time() - t0:.0f} s] {done}/{len(todo)} {r['key']} {r['runtime_s']:.0f} s failed={r['failed']}", flush=True)
                if not fin or time.time() - last > 600:
                    print(f"[{time.time() - t0:.0f} s] heartbeat: {done}/{len(todo)} done, {len(pend)} pending", flush=True)
                last = time.time()
    # assembly (main process; step 2 is fast)
    if not args.no_sim:
        rows = [r for arm in args.arms for rep in range(args.reps) for r in sim_rows(arm, rep, runs)]
        pd.DataFrame(rows).to_csv(out / "results_long.csv", index=False)
        rt = []
        for arm in args.arms:
            for rep in range(args.reps):
                for m in ("clocksync",) + A.METHODS:
                    r = pd.read_pickle(runs / f"{arm}_{rep:03d}_{m}.pkl")
                    rt.append(dict(arm=arm, rep=rep, method=m, runtime_s=round(r["runtime_s"], 1), failed=r["failed"], error=r["error"]))
        pd.DataFrame(rt).to_csv(out / "runtime.csv", index=False)
        print(pd.DataFrame(rt).groupby("method").runtime_s.describe().round(1).to_string(), flush=True)
    if not args.no_real:
        lg, cl = real_rows(runs)
        lg.to_csv(out / "real_long.csv", index=False)
        cl.to_csv(out / "real_clocksync.csv", index=False)
        print(lg.round(3).to_string(), flush=True)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
