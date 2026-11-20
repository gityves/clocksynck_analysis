"""A336 part 1, A235 H235.1 (T2D and obesity do not shift the fast-death metabolic module, equivalence within +-1 h).

Original: `A235_metabolic_offset/scripts/compute_a235.py` (imported, not edited) on the A231 machinery (`matched_label`):
A164 matching, fast cohort, module = liver, adiposes, muscle, A164 host. Offset = case minus control (> 0 advanced), tissue
weights from the permutation null, interval est +- 1.645 max(bootstrap SE, null SD) (`matched_label.pooled_tests`, set
"fast", readout offset_h), equivalent_1h if inside +-1 h. Only the waveform changes (`cs.Dataset(..., waveform=...)`).

Order: (1) 14 genes: observed, 200 matched-set bootstraps (seeds `seed_of("A231A233", design, "", "fast", 0) + 10**6 + i`)
and the first 8 null draws recomputed and compared (1e-6) with the stored `score_tests.csv` and the committed
`output/draws`; the remaining 14-gene null draws are the committed ones; stop if not reproduced. (2) core: 1,000 permutations
(same seeds) and 200 bootstraps on exactly the same resamples. (3) trio contribution = per-draw pooled offset 14 genes minus
core (tissue weights of each waveform's null), 90 % interval of the 200 shared draws. (4) gene drops (100 random + trio,
clock statistic only: observed pooled offset with the core's tissue weights fixed, so the trio drop equals the core value);
direction +1 (advance is the registered prediction). (5) trio alone: `Dataset.gate(n_perm=1000, seed=336)`, p_pairing and
p_acrophase < 0.05; value = trio pooled offset with the core's weights. Core "holds" if its 90 % interval lies inside +-1 h.
One id per condition (A235_H235.1_t2d, A235_H235.1_obese). The pseudo-clock specificity companion and power are not rerun.
GTEx: aggregates only.

    python p1_a235.py [--check] [--jobs 4]       (long: use paper/common/run.sh)
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
sys.path.insert(0, str(PAPER / "exploratory/A235_metabolic_offset/scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_a235 as M5  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ML, A, B = M5.ML, M5.A, M5.B
ENTRY = "a235"
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a235.py"
OUTA = PAPER / "exploratory/A235_metabolic_offset/output"
STORED = pd.read_csv(OUTA / "score_tests.csv")
STORED = STORED[STORED.host == "A164"].set_index("condition")
BASE = C.WF
G = {}


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def ds_of(pb, anchors, wf):
    return cs.Dataset(pb, anchors=anchors, waveform=wf)


def boot_chunk(pb, anchors, plan, wf, seeds):
    """`matched_label.boot_chunk` with an explicit waveform: (replicates, tissues, 2)."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    cases = sorted({s[0] for v in plan.values() for s in v})
    by_case = {t: {s[0]: s for s in plan[t]} for t in plan}
    out = []
    for sd in seeds:
        draw = np.random.default_rng(sd).choice(len(cases), len(cases))
        sets, mp = {t: [] for t in plan}, set()
        for j, i in enumerate(draw):
            for t in plan:
                s = by_case[t].get(cases[i])
                if s is not None:
                    sets[t].append((s, j))
                    mp.update((d, f"{d}@{j}") for d in s)
        m = pd.DataFrame(sorted(mp), columns=["orig", "copy"])
        rep = pb.merge(m, left_on="donor", right_on="orig")
        rep["donor"] = rep["copy"]
        rep = rep.drop(columns=["orig", "copy"]).reset_index(drop=True)
        ds = ds_of(rep, [a for a in anchors if a in set(rep.compartment)], wf)
        out.append([B.stat2(ds, pd.concat([A.labels([s], copy=j) for s, j in sets[t]]), t) if sets[t] else (np.nan, np.nan)
                    for t in plan])
    return np.array(out)


def perm_core(ds, plan, n, jobs, seed0, label):
    return A.run_chunks(B.perm_chunk, n, jobs, seed0, label, ds, plan)


def obs_offsets(ds, plan, msel):
    """Observed per-tissue offsets of the module tissues (array over msel, in plan order)."""
    tl = list(plan)
    return np.array([B.stat2(ds, A.labels(plan[tl[i]]), tl[i])[1] for i in msel], float)


def summarise(name, plans, dss, nulls, boots, msel):
    """Row of `pooled_tests` (fast, offset_h) and per-draw pooled offsets, tissue weights of the null."""
    obs, _ = ML.observed(dss, plans)
    P = ML.pooled_tests(obs, nulls, boots, msel, "module", name)
    r = P[(P.set == "fast") & (P.readout == "offset_h")].iloc[0]
    tw = B.tissue_weights(nulls["fast"][:, msel["fast"], 1])
    draws = B.wmean(boots["fast"][:, msel["fast"], 1], tw)
    return r, tw, draws


def drop_task(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    w = G["drops"][i][1]
    res = {}
    for name, x in G["designs"].items():
        ds = ds_of(x["pb"], x["anchors"], w)
        res[name] = obs_offsets(ds, x["plan"], x["msel"])
    return i, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    warnings.filterwarnings("ignore", message="below the paper settings")
    t0 = time.time()
    core, trio = C.core(BASE), C.trio(BASE)
    cohorts = ML.load_cohorts(500 if chk else None, 0)
    cohorts = {"fast": cohorts["fast"]}
    n_boot, n_perm = (4, 12) if chk else (200, 1000)
    D = {}
    for name in M5.CONDITIONS:
        d = M5.design_of(name)
        plans, dss, anchors, pbs, Mt, _ = ML.build(d, cohorts, quick=chk, tag="gtex_a235")
        if "fast" not in plans or not ML.selection(plans, M5.MET):
            print("no module tissues", name, "(check)" if chk else "")
            continue
        D[name] = dict(plans=plans, anchors=anchors, pb=pbs["fast"], plan=plans["fast"], msel=ML.selection(plans, M5.MET),
                       s0=A.seed_of("A231A233", name, "", "fast", 0))
    log("matching done", list(D))
    # ---- 1. 14 genes -----------------------------------------------------------------------------------------------
    rep = {}
    for name, x in D.items():
        x["ds14"] = ds_of(x["pb"], x["anchors"]["fast"], BASE)
        boot14 = A.run_chunks(boot_chunk, n_boot, a.jobs, x["s0"] + 10 ** 6, f"boot14 {name}", x["pb"], x["anchors"]["fast"],
                              x["plan"], BASE)
        x["boot14"] = boot14
        if chk:
            x["null14"] = A.run_chunks(B.perm_chunk, n_perm, a.jobs, x["s0"], f"null14 {name}", x["ds14"], x["plan"])
        else:
            x["null14"] = np.load(OUTA / "draws" / f"null_{name}_fast.npy")
            head = A.run_chunks(B.perm_chunk, 8, a.jobs, x["s0"], f"null14 head {name}", x["ds14"], x["plan"])
        x["r14"], x["tw14"], x["d14"] = summarise(name, {"fast": x["plan"]}, {"fast": x["ds14"]}, {"fast": x["null14"]},
                                                  {"fast": boot14}, x["msel"])
        if not chk:
            st = STORED.loc[name]
            rep[name] = bool(np.isclose(x["r14"].estimate, st.fast_offset_h, atol=1e-6, rtol=0)
                             and np.isclose(x["r14"].lo90, st.fast_lo90, atol=1e-6, rtol=0)
                             and np.isclose(x["r14"].hi90, st.fast_hi90, atol=1e-6, rtol=0)
                             and np.allclose(boot14, np.load(OUTA / "draws" / f"boot_{name}_fast.npy"), atol=1e-6, equal_nan=True)
                             and np.allclose(head, x["null14"][:8], atol=1e-6, equal_nan=True))
    if not chk:
        log("reproduction", rep)
        if not all(rep.values()):
            C.write_rows(ENTRY, [dict(id=f"A235_H235.1_{n}", figure_panel="Fig 6", statistic="fast module offset (h)", direction=1,
                                      value_14_stored=STORED.loc[n].fast_offset_h, value_14=D[n]["r14"].estimate, reproduced=False,
                                      notes=f"14-gene estimate / interval / bootstrap / null draws not reproduced: {rep}", script=SCRIPT)
                                 for n in D])
            raise SystemExit(f"A235 not reproduced: {rep}")
    # ---- 2. core -----------------------------------------------------------------------------------------------------
    for name, x in D.items():
        x["dsc"] = ds_of(x["pb"], x["anchors"]["fast"], core)
        nullc = perm_core(x["dsc"], x["plan"], n_perm, a.jobs, x["s0"], f"null core {name}")
        bootc = A.run_chunks(boot_chunk, n_boot, a.jobs, x["s0"] + 10 ** 6, f"boot core {name}", x["pb"], x["anchors"]["fast"],
                             x["plan"], core)
        x["nullc"], x["bootc"] = nullc, bootc
        x["rc"], x["twc"], x["dc"] = summarise(name, {"fast": x["plan"]}, {"fast": x["dsc"]}, {"fast": nullc}, {"fast": bootc},
                                               x["msel"])
        log(name, "core done")
    # ---- 4. gene drops -----------------------------------------------------------------------------------------------
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [q for q in dl if q[0] == tl] + dl[:2]
    G.update(drops=dl, designs={n: dict(pb=x["pb"], anchors=x["anchors"]["fast"], plan=x["plan"], msel=x["msel"]["fast"])
                                for n, x in D.items()})
    with Pool(a.jobs) as p:
        dres = dict(p.imap_unordered(drop_task, range(len(dl))))
    dv = {n: {lab: float(B.wmean(dres[i][n][None], D[n]["twc"])[0]) for i, (lab, _) in enumerate(dl)} for n in D}
    for n in D:
        assert np.isclose(dv[n][tl], float(D[n]["rc"].estimate), equal_nan=True)
    log("drops done")
    wfs = dict(dl)
    # ---- 5. trio alone -----------------------------------------------------------------------------------------------
    for name, x in D.items():
        g = ds_of(x["pb"], x["anchors"]["fast"], trio).gate(n_perm=50 if chk else 1000, seed=336)
        x["gp"] = (float(g.p_pairing), float(g.p_acrophase))
        x["gate"] = bool(x["gp"][0] < 0.05 and x["gp"][1] < 0.05)
        x["trio_v"] = float(B.wmean(obs_offsets(ds_of(x["pb"], x["anchors"]["fast"], trio), x["plan"], x["msel"]["fast"])[None],
                                    x["twc"])[0]) if x["gate"] else np.nan
    if chk:
        print("check ok")
        return
    rows = []
    for name, x in D.items():
        r14, rc = x["r14"], x["rc"]
        status = "holds" if (rc.lo90 > -1 and rc.hi90 < 1) else "falls"
        pr = C.drop_percentiles(dv[name], tl, wfs, direction=1)
        diff = x["d14"] - x["dc"]
        rows.append(dict(
            id=f"A235_H235.1_{name}", figure_panel="Fig 6 (chronic disease: no metabolic phase shift)",
            statistic=f"fast-death module offset (h), case minus control, {name}; equivalence", direction=1,
            value_14_stored=STORED.loc[name].fast_offset_h, value_14=float(r14.estimate), reproduced=True,
            value_core=float(rc.estimate),
            core_criterion="null (equivalence): core 90 % interval (est +- 1.645 max(bootstrap SE, null SD)) inside +-1 h",
            core_status=status, diff_14_minus_core=float(r14.estimate - rc.estimate), diff_lo90=float(np.nanquantile(diff, .05)),
            diff_hi90=float(np.nanquantile(diff, .95)),
            diff_draws="200 matched-set bootstraps (A235 seeds), shared resamples, tissue weights of each waveform's own null",
            **pr, trio_gate_passes=x["gate"], trio_alone_value=x["trio_v"], reading=C.reading(status, pr["trio_drop_below_q10"]),
            blind="core not blind (A335 / A329)",
            notes=(f"core 90% interval {rc.lo90:.4f} to {rc.hi90:.4f} (bootstrap SE {rc.boot_se:.4f}, null SD {rc.null_sd:.4f}, "
                   f"own MDE_80 {rc.mde_80:.3f}); 14 genes {r14.lo90:.4f} to {r14.hi90:.4f}; core one-sided-free two-sided p {rc.p_perm:.4f} "
                   f"(14 genes {r14.p_perm:.4f}); drops: clock statistic only, tissue weights fixed at the core null's; trio gate "
                   f"p_pairing/p_acrophase {x['gp'][0]:.4f}/{x['gp'][1]:.4f} (fast cohort); pseudo-clock companion and power not rerun"),
            script=SCRIPT))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws", "figure_panel", "statistic"]).round(4).T.to_string())
    for r in rows:
        print(r["id"], r["notes"])
    log(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
