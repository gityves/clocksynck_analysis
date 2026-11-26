"""A385d (descriptive check of A385b, registered 2026-11-26): does the A385b design protect the sheep person leverage
against a waveform mismatch under working-hours sampling, and does the sheep value depend on the waveform?

Part 1 (simulation of the sheep design, every animal; D1 reads the big projects): clock-gene counts redrawn with true phases in windows of 4-24 h per project,
truth T0 (human waveform) or T1 (A400c mouse pattern), with or without a person factor; read with the human waveform and
the A385b D1 statistic. Part 2 (real sheep): D1 leverage with the sheep read by the human (H), baboon (B, fitted on
Mure 2018 with ZT) and mouse (Mo, A400c Mf) waveforms; sheep-permuted null; 50 + 50 companions per waveform.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A385_sheep_person/scripts/a385d.py --jobs 8
    python a385d.py --check      # 1 replicate, 2 sets, 20 permutations; nothing written
Outputs (`../output/`): d_simulation.csv, d_waveforms.csv, d_leverages.csv, d_summary.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
for p in (PAPER, HERE.parent, EXP / "A394_mouse_organisation" / "scripts", EXP / "A395_baboon_timed_gain" / "scripts"):
    sys.path.insert(0, str(p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a385_a386 as M  # noqa: E402
import a385b as R  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform, fit_waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
GENES = list(M.WF.genes)
WIDTHS, TRUTHS, FACTORS = (4, 8, 12, 24), ("T0", "T1"), (0.0, 0.25)
N_REP, N_SETS, N_PERM = 10, 50, 1000
S = {}


# ---- waveforms -------------------------------------------------------------------------------------------------------
def log_frac(d):
    return np.log((d[GENES].to_numpy(float) + 0.5) / d.depth.to_numpy(float)[:, None])


def mouse_waveform():
    """A400c Mf: all mice, peripheral tissues kept by the A394 QC, group tissue x arm."""
    import a394
    pb, per, _ = a394.mouse_table()
    d = pb[pb.compartment.isin(per)]
    zt = d.donor.astype(str).str.extract(r"_ZT(\d+)_", expand=False).astype(float).to_numpy()
    return fit_waveform(log_frac(d), zt, GENES, group=d.compartment.astype(str) + "|" + d.stratum.astype(str), name="mouse")


def baboon_waveform():
    import a395
    pb, _, zt = a395.load()
    return fit_waveform(log_frac(pb), zt.reindex(pb.donor).to_numpy(float), GENES, group=pb.compartment.astype(str),
                        name="baboon")


def as_set(wf, genes, rd, name):
    w = Waveform(tuple(genes), np.asarray(wf.amp), np.asarray(wf.psi), None, None, name)
    return w if rd == "14 genes" else Waveform(tuple(w.genes[k] for k in M.CORE_IDX), w.amp[M.CORE_IDX],
                                                w.psi[M.CORE_IDX], None, None, name)


# ---- D1 leverage -----------------------------------------------------------------------------------------------------
def states(pb, w):
    ds = cs.Dataset(pb[M.META + list(w.genes)], anchors=S["per"], waveform=w, mixed_strata=True, name=w.name,
                    min_donors=5)
    return R.states(ds, M.halvings(S["per"]))


def d1(st, rng=None):
    """A385b D1: within project (>= 30 animals), Fisher-combined, median over halvings; rng: permute B within project."""
    L = []
    for j, a, b, cov in st:
        pr = S["project"].reindex(j).to_numpy()
        z = w = 0.0
        for p in np.unique(pr):
            m = np.flatnonzero(pr == p)
            if len(m) < M.MIN_PROJECT:
                continue
            bm = b[m] if rng is None else b[m][rng.permutation(len(m))]
            z += (len(m) - 3) * np.arctanh(M.leverage(a[m], bm, cov[m]))
            w += len(m) - 3
        L.append(np.tanh(z / w) if w > 0 else np.nan)
    return float(np.nanmedian(L))


# ---- simulation ------------------------------------------------------------------------------------------------------
def simulate(truth, width, factor, seed):
    rng = np.random.default_rng(seed)
    base = S["sim"]
    don = base.donor.to_numpy()
    animals = np.unique(don)
    proj = S["project"].reindex(animals).to_numpy()
    centre = {p: rng.uniform(0, 24) for p in np.unique(proj)}
    t = pd.Series([centre[p] for p in proj], index=animals) + rng.uniform(-width / 2, width / 2, len(animals))
    lam = pd.Series(np.exp(rng.normal(0, factor, len(animals))) if factor > 0 else 1.0, index=animals)
    wf = S["truth"][truth]
    x = 2 * np.pi * t.reindex(don).to_numpy()[:, None] / 24
    osc = lam.reindex(don).to_numpy()[:, None] * wf.amp[None, :] * np.cos(x - wf.psi[None, :])
    sd = np.sqrt(np.maximum(S["s"] ** 2 - wf.amp[None, :] ** 2 / 2, S["s"] ** 2 / 4))
    mu = S["m"] + osc + rng.normal(size=osc.shape) * sd
    out = base[M.META].copy()
    out[GENES] = rng.poisson(base.depth.to_numpy(float)[:, None] * np.exp(mu))
    return out


def sim_task(job):
    warnings.filterwarnings("ignore")
    truth, width, factor, r = job
    try:
        pb = simulate(truth, width, factor, [385, 7, TRUTHS.index(truth), width, int(factor * 100), r])
        return [dict(truth=truth, width_h=width, person_sd=factor, rep=r, readout=rd, leverage=d1(states(pb, w)))
                for rd, w in (("14 genes", M.WF), ("core", M.CORE))]
    except Exception as e:  # counted, not dropped silently
        return [dict(truth=truth, width_h=width, person_sd=factor, rep=r, readout=None, error=f"{type(e).__name__}: {e}"[:200])]


# ---- real sheep ------------------------------------------------------------------------------------------------------
def swap_task(job):
    warnings.filterwarnings("ignore")
    kind, i, rd, wname, genes = job
    try:
        st = states(S["sheep"], as_set(S["wf"][wname], genes, rd, f"{wname}{kind}{i}"))
        o = dict(kind=kind, set=i, readout=rd, waveform=wname, leverage=d1(st))
        if kind == "clock":
            rng = np.random.default_rng([385, 8])
            null = np.array([d1(st, rng) for _ in range(S["n_perm"])])
            o.update(p_perm=float((1 + np.sum(null >= o["leverage"])) / (1 + len(null))), null_q95=float(np.quantile(null, 0.95)))
        return [o]
    except Exception as e:  # counted, not dropped silently
        return [dict(kind=kind, set=i, readout=rd, waveform=wname, error=f"{type(e).__name__}: {e}"[:200])]


def task(job):
    return sim_task(job[1:]) if job[0] == "sim" else swap_task(job[1:])


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    sheep, per, _ = M.load()
    sheep["donor"] = sheep.donor.astype(str)
    project = sheep.groupby("donor").stratum.first()
    lm = M.log_mean(sheep)
    n = 2 if a.check else N_SETS
    sets = {k: pseudo_waveforms(lm, N_SETS, np.random.default_rng([385, t]), M.WF, kind=k)[:n]
            for k, t in (("matched", 2), ("responsive", 3))}
    genes = sorted({g for ws in sets.values() for w in ws for g in w.genes} | set(GENES))
    sim = sheep[M.META + GENES].reset_index(drop=True)  # every animal: hosts as in A385b; D1 reads the big projects
    lf = pd.DataFrame(log_frac(sim), columns=GENES)
    key = sim.compartment.astype(str) + "|" + sim.stratum.astype(str)
    wf = {"H": M.WF, "B": baboon_waveform(), "Mo": mouse_waveform()}
    t1 = Waveform(M.WF.genes, wf["Mo"].amp * np.mean(M.WF.amp) / np.mean(wf["Mo"].amp), wf["Mo"].psi, None, None, "T1")
    S.update(sheep=sheep[M.META + genes].copy(), per=per, project=project, sim=sim, wf=wf,
             m=lf.groupby(key.to_numpy()).transform("mean").to_numpy(),
             s=lf.groupby(key.to_numpy()).transform("std").fillna(lf.groupby(sim.compartment.to_numpy())
                                                                  .transform("std")).to_numpy(),
             truth={"T0": M.WF, "T1": t1}, n_perm=20 if a.check else N_PERM)
    del sheep
    W = pd.concat([w.table().add_prefix(f"{k}_") for k, w in wf.items()], axis=1)
    print(f"{project.size} animals, {len(per)} tissues; simulation on {sim.donor.nunique()} animals "
          f"({time.time() - t0:.0f}s)", flush=True)
    reps = 1 if a.check else N_REP
    jobs = [("sim", tr, wd, fc, r) for tr in TRUTHS for wd in WIDTHS for fc in FACTORS for r in range(reps)]
    if a.check:
        jobs = jobs[:2]
    jobs += [("swap", "clock", 0, rd, wn, GENES) for wn in ("H", "B", "Mo") for rd in ("14 genes", "core")]
    jobs += [("swap", k, j, rd, wn, w.genes) for wn in ("B", "Mo") for k, ws in sets.items() for j, w in enumerate(ws)
             for rd in ("14 genes", "core")]
    rows = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            rows += r
            if (k + 1) % 20 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    C = pd.DataFrame(rows)
    sim_rows = C[C.truth.notna()] if "truth" in C else C.iloc[:0]
    sw = C[C.waveform.notna()]
    h = sw[(sw.kind == "clock") & (sw.waveform == "H")].set_index("readout").leverage
    print(f"A385b reproduction (H): 14 genes {h['14 genes']:.4f} vs 0.5577; core {h["core"]:.4f} vs 0.3159", flush=True)
    failed = int(C.error.notna().sum()) if "error" in C else 0
    if a.check:
        print(f"check ok; failed {failed}", flush=True)
        return
    prev = pd.read_csv(OUT / "c_leverages.csv")
    prev = prev[(prev.kind != "clock") & (prev.variant == "unadjusted")].assign(waveform="H")
    comp = pd.concat([sw[sw.kind != "clock"], prev[["kind", "set", "readout", "waveform", "leverage"]]])
    summ = []
    for wn in ("H", "B", "Mo"):
        for rd in ("14 genes", "core"):
            clk = sw[(sw.kind == "clock") & (sw.waveform == wn) & (sw.readout == rd)].iloc[0]
            c = comp[(comp.waveform == wn) & (comp.readout == rd)]
            q = {k: float(np.nanquantile(c[c.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}
            summ.append(dict(waveform=wn, readout=rd, clock=clk.leverage, p_perm=clk.p_perm, null_q95=clk.null_q95,
                             matched_q95=q["matched"], responsive_q95=q["responsive"],
                             above_both=bool(clk.leverage > max(q.values()))))
    simsum = (sim_rows.groupby(["truth", "person_sd", "width_h", "readout"]).leverage
              .describe(percentiles=[0.1, 0.5, 0.9])[["count", "10%", "50%", "90%"]].reset_index())
    OUT.mkdir(parents=True, exist_ok=True)
    C.to_csv(OUT / "d_leverages.csv", index=False)
    simsum.to_csv(OUT / "d_simulation.csv", index=False)
    W.to_csv(OUT / "d_waveforms.csv")
    Sm = pd.DataFrame(summ)
    Sm.to_csv(OUT / "d_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(simsum.round(3).to_string(index=False))
    print(Sm.round(4).to_string(index=False), f"\nfailed: {failed} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
