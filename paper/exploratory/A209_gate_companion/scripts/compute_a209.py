"""A209: false-claim rate of the gate plus the expression-matched companion, the rule real-data claims use.

Registered in paper/docs/registered_experiments.md (A209). Simulator and design = A145 (`A145_two_step/scripts/lib.py`:
anchors, targets, `design`, `psi_of`, `host_variance`, `shared_sd`), with the 14-gene GTEx waveform extended by 200
non-clock genes (amplitude 0, expression-matched mesors). One pickle per (arm, rep) under <out>/runs/ (a rerun resumes).
Per replicate: the clock genes fitted with the gate (forced, as A145), then 100 matched pseudo sets (and 100 responsive
sets in N_resp) renamed to the clock genes and fitted on the same anchors with the gate off.

    RUN_CAP_GB=10 paper/common/run.sh paper/exploratory/A209_gate_companion/scripts/compute_a209.py --jobs 8

Outputs (<out>/): per_target.csv, rates.csv, tests.csv.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import multiprocessing as mp
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent))
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory" / "A145_two_step" / "scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import lib as A145  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parent / "output"
SEED, N_EXTRA, N_RESP_POOL, N_SETS = 209, 200, 30, 100
WF = cs.gtex()
GENES = list(WF.genes)
RESP_CLOCK = ["PER1", "PER2", "NR1D1", "DBP", "TEF"]
EXTRA = [f"PSG{i:03d}" for i in range(N_EXTRA)]
RESP_POOL, PLAIN = EXTRA[:N_RESP_POOL], EXTRA[N_RESP_POOL:]
ARMS = {  # window, jitter, amplitude, iid shared program, responsive factor
    "N": dict(window=(0.0, 24.0), jitter=0.0, amplitude=0.0, shared=True, resp=False),
    "N_resp": dict(window=(0.0, 24.0), jitter=0.0, amplitude=0.0, shared=True, resp=True),
    "W": dict(window=(8.0, 18.0), jitter=1.0, amplitude=1.0, shared=False, resp=False),
}
C_MIN = 0.2


def extended_waveform() -> Waveform:
    """GTEx waveform plus N_EXTRA non-clock genes: amplitude 0, mesors resampled from the clock genes + N(0, 0.3)."""
    r = np.random.default_rng([SEED, 0])
    mes = r.choice(np.asarray(WF.mesor, float), N_EXTRA) + r.normal(0.0, 0.3, N_EXTRA)
    return Waveform(tuple(GENES + EXTRA), np.r_[WF.amp, np.zeros(N_EXTRA)], np.r_[WF.psi, np.zeros(N_EXTRA)],
                    np.r_[WF.mesor, mes], None, "a209")


BIG = extended_waveform()


def cohort(arm: str, rep: int) -> pd.DataFrame:
    """A145 `lib.cohort` on the extended waveform (seed 209), with the arm's shared programs."""
    a, d = ARMS[arm], A145.design(rep)
    V = A145.host_variance(a["window"], a["jitter"])
    comps = {k: Compartment(z=1.0, psi=A145.psi_of(A145.C_ANCHOR, 1.0, V)) for k in A145.ANCHORS}
    for k, r in d.iterrows():
        comps[k] = Compartment(z=r.amp * np.exp(1j * r.offset_h / A145.HOURS), psi=A145.psi_of(r.c, r.amp, V), coupled=r.c > 0)
    sd = A145.shared_sd("large") if a["shared"] else 0.0
    load = None
    if a["resp"]:
        load = np.where(np.isin(BIG.genes, RESP_CLOCK + RESP_POOL), 2 * A145.shared_sd("large"), 0.0)
    tab, _ = cs.simulate(BIG, comps, n=A145.N_DONORS, window_h=a["window"], jitter_h=a["jitter"], amplitude=a["amplitude"],
                         shared_sd=sd, shared_load=load, seed=np.random.SeedSequence([SEED, list(ARMS).index(arm), rep]))
    return tab.rename(columns={"resection": "donor"})


def couplings(tab: pd.DataFrame, genes: list[str], gate: bool):
    """Target couplings of a fit on `genes` renamed to the clock genes (position i stands for clock gene i)."""
    sub = tab[["donor", "compartment", "depth", *genes]].set_axis(["donor", "compartment", "depth", *GENES], axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = cs.Dataset(sub, anchors=A145.ANCHORS, name="a209").fit(n_boot=0, gate=gate, force=True)
    c = f.table.loc[A145.TARGETS, "coupling"].to_numpy(float)
    return (c, f.gate) if gate else c


def pseudo_sets(rng, responsive: bool) -> list[list[str]]:
    sets = []
    for _ in range(N_SETS):
        s = list(rng.choice(PLAIN, len(GENES), replace=False))
        if responsive:
            rp = iter(rng.choice(RESP_POOL, len(RESP_CLOCK), replace=False))
            s = [next(rp) if g in RESP_CLOCK else x for g, x in zip(GENES, s)]
        sets.append(s)
    return sets


def task(arm: str, rep: int, runs: str) -> dict:
    from clocksync.config import limit_threads

    limit_threads(1)
    f = Path(runs) / f"{arm}_{rep:03d}.pkl"
    if f.exists():
        return pd.read_pickle(f)
    t0 = time.time()
    tab = cohort(arm, rep)
    c, g = couplings(tab, GENES, gate=True)
    rng = np.random.default_rng([SEED, list(ARMS).index(arm), rep, 1])
    matched = np.array([couplings(tab, s, gate=False) for s in pseudo_sets(rng, False)])
    resp = np.array([couplings(tab, s, gate=False) for s in pseudo_sets(rng, True)]) if ARMS[arm]["resp"] else None
    out = dict(arm=arm, rep=rep, c=c, true_c=A145.design(rep).c.to_numpy(), gate_pass=bool(g.passes),
               p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share, matched=matched, resp=resp,
               runtime_s=time.time() - t0)
    pd.to_pickle(out, f)
    return out


def summarise(res: list[dict], out: Path):
    from scipy import stats

    rows = []
    for r in res:
        q_m = np.quantile(r["matched"], 0.95, axis=0)
        q_r = np.quantile(r["resp"], 0.95, axis=0) if r["resp"] is not None else np.full(len(r["c"]), -np.inf)
        pct_m = (r["matched"] < r["c"]).mean(0)
        for j, t in enumerate(A145.TARGETS):
            base = r["gate_pass"] and r["c"][j] > C_MIN
            rows.append(dict(arm=r["arm"], rep=r["rep"], target=t, true_c=r["true_c"][j], c=r["c"][j], gate_pass=r["gate_pass"],
                             matched_q95=q_m[j], resp_q95=q_r[j], matched_pct=pct_m[j], claim_gate=base,
                             claim_matched=base and r["c"][j] > q_m[j],
                             claim_full=base and r["c"][j] > q_m[j] and r["c"][j] > q_r[j]))
    d = pd.DataFrame(rows).sort_values(["arm", "rep", "target"])
    d.to_csv(out / "per_target.csv", index=False)
    rates = d.groupby(["arm", "true_c"])[["gate_pass", "claim_gate", "claim_matched", "claim_full"]].mean().reset_index()
    rates.to_csv(out / "rates.csv", index=False)

    def false_rate(arm, col, h):
        x = d[d.arm == arm]
        k, n = int(x[col].sum()), len(x)
        p = float(stats.binomtest(k, n, 0.127, alternative="less").pvalue)
        return dict(h=h, arm=arm, rule=col, k=k, n=n, rate=k / n, p_vs_0127=p, gate_only_rate=float(x.claim_gate.mean()),
                    passes_rule=(k / n <= 0.05) and p < 0.05)

    w = d[(d.arm == "W") & (d.true_c >= 0.5)]
    r_gate, r_comp = float(w.claim_gate.mean()), float(w.claim_matched.mean())
    tests = pd.DataFrame([
        false_rate("N", "claim_matched", "H209.1"),
        false_rate("N_resp", "claim_full", "H209.2"),
        dict(h="H209.3", arm="W", rule="claim_matched / claim_gate (true c >= 0.5)", k=int(w.claim_matched.sum()), n=len(w),
             rate=r_comp, gate_only_rate=r_gate, ratio=r_comp / r_gate if r_gate > 0 else np.nan,
             passes_rule=r_gate > 0 and r_comp >= 0.8 * r_gate),
    ])
    p = tests.p_vs_0127.fillna(0.0).to_numpy()
    order, adj, run = np.argsort(p), np.empty(len(p)), 0.0
    for i, j in enumerate(order):
        run = max(run, min(1.0, (len(p) - i) * p[j]))
        adj[j] = run
    tests["p_holm"] = np.where(tests.p_vs_0127.notna(), adj, np.nan)
    tests["verdict"] = np.where(tests.passes_rule & (tests.p_holm.fillna(0) < 0.05), "PASS", "FAIL")
    tests.to_csv(out / "tests.csv", index=False)
    print(tests.to_string(), flush=True)
    print(rates.round(3).to_string(), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    out = Path(args.out)
    runs = out / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    todo = [(a, r) for r in range(args.reps) for a in ARMS]
    print(f"workers = {args.jobs}; {len(todo)} tasks", flush=True)
    t0, res, last = time.time(), [], time.time()
    with cf.ProcessPoolExecutor(args.jobs, mp_context=mp.get_context("spawn")) as ex:
        for i, fu in enumerate(cf.as_completed([ex.submit(task, a, r, str(runs)) for a, r in todo])):
            res.append(fu.result())
            if time.time() - last > 300 or i + 1 == len(todo):
                print(f"[{time.time() - t0:.0f} s] {i + 1}/{len(todo)}", flush=True)
                last = time.time()
    summarise(res, out)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
