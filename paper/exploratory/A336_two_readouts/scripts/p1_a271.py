"""A336 part 1, A271 H271.2 (GTEx periphery rho between tissue clock amplitude and coupling, A104 rows, 27 tissues).

Original: `A271_responsive_companions/scripts/score_a271_a104.py` (imported, not edited; `A104.per_tissue`, `A104.rho`).
Periphery block only. Same rows, fit (every tissue an anchor, leave-one-out host), seeds [271, 1] matched / [271, 2]
responsive (200 each, `pseudo_waveforms` on the genome-wide log mean fraction); only the waveform changes. Companions:
the original 14-gene draws restricted to the core positions. 14 genes: rho_clock recomputed and the original set draws
verified gene by gene against `a104_sets.csv` (the 400 companion rhos are not refitted at 14 genes: stored, see report).
Criterion: rho_core above the q95 of the matched and of the responsive core companions (holds: both; at the line: one).
Trio contribution: 200 tissue bootstraps (the statistic's unit; rng [336, 1]) of rho, shared. Gene drops: 100 random + trio
(seed 336). Trio alone: `Dataset.gate(n_perm=1000, seed=336)` on the periphery rows, p_pairing and p_acrophase < 0.05.
GTEx: aggregates only.

    python p1_a271.py [--check] [--jobs 4]
Output: ../output/part1/a271.csv
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/backbone"))
sys.path.insert(0, str(PAPER / "exploratory/A271_responsive_companions/scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
import score_a271_a104 as A  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

ENTRY = "a271"
BASE = C.WF
META = A.META
OUTA = PAPER / "exploratory/A271_responsive_companions/output"
STORED = pd.read_csv(OUTA / "a104_tests.csv").set_index("block").loc["periphery"]
SETS = pd.read_csv(OUTA / "a104_sets.csv")
G = {}
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a271.py"


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def fit(w):
    """Periphery A104 per-tissue table for waveform w (every tissue an anchor)."""
    pb = G["pb"]
    comps = sorted(c for c in pb.compartment.unique() if not c.startswith("Brain"))
    ds = cs.Dataset(pb.loc[pb.compartment.isin(comps), META + list(w.genes)], anchors=comps, waveform=w, name="periphery")
    return A104.per_tissue(ds), ds


def task(i):
    warnings.filterwarnings("ignore")
    w = G["wfs"][i]
    try:
        pt, _ = fit(w)
        return i, pt[["own_amplitude", "coupling"]]
    except Exception:  # noqa: BLE001
        return i, None


def run(wfs, jobs):
    G["wfs"] = wfs
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        out = dict(p.imap_unordered(task, range(len(wfs))))
    return [out[i] for i in range(len(wfs))]


def rho_of(pt):
    return A104.rho(pt) if pt is not None else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    n = 3 if chk else 200
    pb = A.load()
    limit_threads(1)
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    core, trio = C.core(BASE), C.trio(BASE)
    sets = {k: pseudo_waveforms(lm, n, np.random.default_rng([271, j]), BASE, k) for j, k in ((1, "matched"), (2, "responsive"))}
    # ---- 14 genes: draws verified, rho_clock reproduced -------------------------------------------------------------
    st = SETS[SETS.block == "periphery"]
    for k, ws in sets.items():
        s = st[st.kind == k].sort_values("set")
        want = [g for g in s.genes.tolist()]
        got = [";".join(w.genes) for w in ws]
        assert want[:n] == got[:n], f"{k} draws differ from the stored sets"
    used = sorted({g for ws in sets.values() for w in ws for g in w.genes} | set(BASE.genes))
    G["pb"] = pb[META + used].copy()
    del pb
    pt14, _ = fit(BASE)
    r14 = A104.rho(pt14)
    if not chk and not np.isclose(r14, STORED.rho_clock, atol=1e-6):
        C.write_rows(ENTRY, [dict(id="A271_H271.2", figure_panel="Fig 6/GTEx", statistic="periphery rho(A, c)", direction=1,
                                  value_14_stored=STORED.rho_clock, value_14=r14, reproduced=False,
                                  notes="rho_clock not reproduced", script=SCRIPT)])
        raise SystemExit("A271 not reproduced")
    log("14-gene stage done")
    # ---- core, companions, drops, trio (one pool pass) --------------------------------------------------------------
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    wfs = [core] + [C.restrict(w, BASE, core) for k in ("matched", "responsive") for w in sets[k]] + [w for _, w in dl]
    t0 = time.time()
    res = run(wfs, a.jobs)
    log(f"fits done {(time.time() - t0) / 60:.1f} min")
    ptc = res[0]
    rc = rho_of(ptc)
    comp = {"matched": [rho_of(p) for p in res[1:1 + n]], "responsive": [rho_of(p) for p in res[1 + n:1 + 2 * n]]}
    dv = {lab: rho_of(p) for (lab, _), p in zip(dl, res[1 + 2 * n:])}
    assert np.isclose(dv[tl], rc, equal_nan=True)
    dwf = dict(dl)
    # gate on the trio
    _, dst = fit(trio)
    g = dst.gate(n_perm=50 if chk else 1000, seed=336)
    gate = bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05)
    tv = rho_of(fit(trio)[0]) if gate else np.nan
    # tissue bootstrap, shared
    rb = np.random.default_rng([336, 1])
    nt = len(pt14)
    bidx = [rb.integers(0, nt, nt) for _ in range(5 if chk else 200)]
    assert list(pt14.index) == list(ptc.index)
    diff = np.array([A104.rho(pt14.iloc[i]) - A104.rho(ptc.iloc[i]) for i in bidx])
    bc = np.array([A104.rho(ptc.iloc[i]) for i in bidx])
    if chk:
        print("check ok")
        return
    q95 = {k: float(np.nanquantile(v, 0.95)) for k, v in comp.items()}
    conds = [bool(rc > q95["matched"]), bool(rc > q95["responsive"])]
    status = "holds" if all(conds) else "at the line" if any(conds) else "falls"
    pr = C.drop_percentiles(dv, tl, dwf, direction=1)
    row = dict(id="A271_H271.2", figure_panel="Fig 6/GTEx (stronger clocks more coupled)",
               statistic="periphery rho between tissue clock amplitude and coupling (27 tissues)", direction=1,
               value_14_stored=STORED.rho_clock, value_14=r14, reproduced=True, value_core=rc,
               core_criterion="rho above q95 of the matched and of the responsive companions (200 each, restricted to core positions)",
               core_status=status, diff_14_minus_core=r14 - rc, diff_lo90=float(np.quantile(diff, .05)),
               diff_hi90=float(np.quantile(diff, .95)), diff_draws="200 tissue bootstraps (rng [336, 1]), shared", **pr,
               trio_gate_passes=gate, trio_alone_value=tv, reading=C.reading(status, pr["trio_drop_below_q10"]), blind="core blind",
               notes=(f"core companions: matched median {np.nanmedian(comp['matched']):.3f} q95 {q95['matched']:.3f} (valid {int(np.isfinite(comp['matched']).sum())}), "
                      f"responsive median {np.nanmedian(comp['responsive']):.3f} q95 {q95['responsive']:.3f} (valid {int(np.isfinite(comp['responsive']).sum())}); "
                      f"core 90 % tissue bootstrap {np.nanquantile(bc, .05):.3f} to {np.nanquantile(bc, .95):.3f}; "
                      f"trio gate p_pairing {g.p_pairing:.4f}, p_acrophase {g.p_acrophase:.4f}; 14-gene companion rhos not refitted (draws verified)"),
               script=SCRIPT)
    df = C.write_rows(ENTRY, [row])
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws"]).round(4).T.to_string())
    print(row["notes"])


if __name__ == "__main__":
    main()
