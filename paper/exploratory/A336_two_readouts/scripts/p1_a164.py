"""A336 part 1, A164 H164.1 (T2D) and H164.2 (obesity): pooled metabolic-module coupling d*, equivalence within +-0.15.

Original: `A164_gtex_metabolic_matched/scripts/compute_a164.py` (imported, not edited; matching, bootstrap seeds
`seed_of("A164", cohort, design, 0) + 10**6 + i`, 500 matched-set bootstraps, inverse-variance pooling over vent and
fast, interval = D +- 1.645 SE). Only the waveform changes (`cs.Dataset(..., waveform=...)`). The equivalence needs no
permutation null (the 14-gene nulls feed only H164.3 and the power table; power stages not rerun), so the core run is the
observed contrast and the 500 bootstraps on the core waveform, on exactly the same resamples as the 14-gene ones.

Order: (1) 14 genes recomputed (observed + 500 bootstraps) and compared with the stored tests table and committed
`output/draws/boot_*.npy` (1e-6); stop if not reproduced; (2) core bootstraps; (3) trio contribution = per-draw pooled D
(cohort weights fixed at each waveform's own bootstrap-SE weights) 14 genes minus core, 90 % interval of the 500 shared
draws; (4) gene drops (100 random + trio, clock statistic only, observed pooled D with the core's cohort weights fixed, so the
trio drop equals the core D); direction -1 (a loss of coupling is the registered concern of the null); (5) trio alone (gate
`Dataset.gate(n_perm=1000, seed=336)` in vent and fast, both p < 0.05).
Core "holds" if its 90 % interval lies inside +-0.15, else "falls". GTEx: aggregates only.

    python p1_a164.py [--check] [--jobs 4]       (long: use paper/common/run.sh)
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
sys.path.insert(0, str(PAPER / "exploratory/A164_gtex_metabolic_matched/scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_a164 as M  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a164"
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a164.py"
OUTA = PAPER / "exploratory/A164_gtex_metabolic_matched/output"
IDS = {"t2d": "A164_H164.1", "obese": "A164_H164.2"}
BASE = C.WF
G = {}


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def ds_of(pb, anchors, wf):
    return cs.Dataset(pb, anchors=anchors, waveform=wf)


def boot_chunk(pb, anchors, plan, wf, seeds):
    """`compute_a164.boot_chunk` with an explicit waveform."""
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
        out.append([M.dstar(ds, pd.concat([M.labels([s], copy=j) for s, j in sets[t]]), t) for t in plan])
    return np.array(out)


def obs_module(ds, plan):
    """Observed module mean d* of a cohort (finite MET tissues)."""
    return M.stats_of({t: M.dstar(ds, M.labels(plan[t]), t) for t in plan if t in M.MET})[0]


def drop_task(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    w = G["drops"][i][1]
    res = {}
    for (c, dsg), plan in G["plans"].items():
        if any(t in M.MET for t in plan):
            res[(c, dsg)] = obs_module(ds_of(G["pbs"][c], G["anchors"][c], w), plan)
    return i, res


def pooled(est, bt_cols):
    """(D, se_p, weights) from cohort estimates and the (n,) module draws of each cohort."""
    se = {c: float(np.nanstd(v, ddof=1)) for c, v in bt_cols.items()}
    return M.pool(est, se)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    warnings.filterwarnings("ignore", message="below the paper settings")
    t0 = time.time()
    sub = 40 if chk else None
    n_boot, min_cases = (4, 3) if chk else (500, 10)
    core, trio = C.core(BASE), C.trio(BASE)
    ph = M.K.protected()
    pbs, plans, anchors = {}, {}, {}
    for c, h in M.COH.items():
        pb, don = M.K.cohort(h, subsample=sub, seed=0)
        pb = pb[pb.compartment.isin(set(M.S.TISSUES) | set(M.TISS))].reset_index(drop=True)
        q = don[["hardy", "male", "ischaemic_h"]].join(ph[["age", "bmi", "t2d", "steroid"]]).dropna(subset=["age", "ischaemic_h"])
        td = {t: set(pb.donor[pb.compartment == t]) for t in M.TISS}
        pbs[c] = pb
        anchors[c] = [x for x in M.S.TISSUES if x in set(pb.compartment)]
        for dsg in M.PRIMARY:
            plans[c, dsg] = M.plan_of(q, td, dsg, min_cases)[0]
    keys = [k for k, p in plans.items() if any(t in M.MET for t in p)]
    log("matching done", sorted(keys))
    seed0 = {k: M.seed_of("A164", k[0], k[1], 0) + 10 ** 6 for k in keys}

    def boots(wf, label):
        return {k: M.run_chunks(boot_chunk, n_boot, a.jobs, seed0[k], f"{label} {k}", pbs[k[0]], anchors[k[0]], plans[k], wf)
                for k in keys}

    def summary(dss, raw):
        """per design: (D, se_p, weights, per-draw pooled draws) for pooled module d*."""
        out = {}
        for dsg in M.PRIMARY:
            ks = [k for k in keys if k[1] == dsg]
            est = {c: obs_module(dss[c], plans[c, dsg]) for c, _ in ks}
            bt = {c: M.modules(raw[c, dsg], list(plans[c, dsg]))[:, 0] for c, _ in ks}
            D, se, w = pooled(est, bt)
            out[dsg] = dict(D=D, se=se, w=w, est=est, draws=sum(w[c] * bt[c] for c in w), lo=D - M.Z90 * se, hi=D + M.Z90 * se)
        return out

    # ---- 1. 14 genes, reproduction -----------------------------------------------------------------------------------
    dss14 = {c: ds_of(pbs[c], anchors[c], BASE) for c in M.COH}
    raw14 = boots(BASE, "boot14")
    S14 = summary(dss14, raw14)
    if not chk:
        hyp = pd.read_csv(OUTA / "a164_tests_hypotheses.csv").set_index("design")
        rep = {}
        for dsg in M.PRIMARY:
            r = hyp[hyp.hypothesis == IDS[dsg].split("_")[1]].loc[dsg]
            rep[dsg] = bool(np.isclose(S14[dsg]["D"], r.D, atol=1e-6, rtol=0) and np.isclose(S14[dsg]["lo"], r.lo90, atol=1e-6, rtol=0)
                            and np.isclose(S14[dsg]["hi"], r.hi90, atol=1e-6, rtol=0)
                            and all(np.allclose(raw14[k], np.load(OUTA / "draws" / f"boot_{k[0]}_{k[1]}.npy"), atol=1e-6,
                                                equal_nan=True) for k in keys if k[1] == dsg))
        log("reproduction", rep)
        if not all(rep.values()):
            C.write_rows(ENTRY, [dict(id=IDS[d], figure_panel="Fig 6", statistic="pooled module d*", direction=-1,
                                      value_14_stored=hyp[hyp.hypothesis == IDS[d].split("_")[1]].loc[d].D, value_14=S14[d]["D"],
                                      reproduced=False, notes=f"14-gene value / interval / bootstrap draws not reproduced: {rep}",
                                      script=SCRIPT) for d in M.PRIMARY])
            raise SystemExit("A164 not reproduced")
    # ---- 2. core -----------------------------------------------------------------------------------------------------
    dssc = {c: ds_of(pbs[c], anchors[c], core) for c in M.COH}
    rawc = boots(core, "bootcore")
    SC = summary(dssc, rawc)
    log("core bootstraps done")
    # ---- 4. gene drops -----------------------------------------------------------------------------------------------
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [q for q in dl if q[0] == tl] + dl[:2]
    G.update(drops=dl, plans=plans, pbs=pbs, anchors=anchors)
    with Pool(a.jobs) as p:
        dres = dict(p.imap_unordered(drop_task, range(len(dl))))
    dv = {dsg: {} for dsg in M.PRIMARY}
    for i, (lab, _) in enumerate(dl):
        for dsg in M.PRIMARY:
            w = SC[dsg]["w"]
            vals = {c: dres[i][(c, dsg)] for c in w}
            dv[dsg][lab] = sum(w[c] * vals[c] for c in w) if all(np.isfinite(list(vals.values()))) else np.nan
    for dsg in M.PRIMARY:
        assert not SC[dsg]["w"] or np.isclose(dv[dsg][tl], SC[dsg]["D"], equal_nan=True)
    log("drops done")
    wfs = dict(dl)
    # ---- 5. trio alone -----------------------------------------------------------------------------------------------
    gates = {c: ds_of(pbs[c], anchors[c], trio).gate(n_perm=50 if chk else 1000, seed=336) for c in M.COH}
    gp = {c: (float(g.p_pairing), float(g.p_acrophase)) for c, g in gates.items()}
    gate_ok = all(p1 < 0.05 and p2 < 0.05 for p1, p2 in gp.values())
    dsst = {c: ds_of(pbs[c], anchors[c], trio) for c in M.COH} if gate_ok else None
    if chk:
        print("check ok")
        return
    rows = []
    for dsg in M.PRIMARY:
        s14, sc = S14[dsg], SC[dsg]
        st0 = hyp[hyp.hypothesis == IDS[dsg].split("_")[1]].loc[dsg]
        holds = bool(sc["lo"] > -M.BOUND and sc["hi"] < M.BOUND)
        status = "holds" if holds else "falls"
        pr = C.drop_percentiles(dv[dsg], tl, wfs, direction=-1)
        diff = s14["draws"] - sc["draws"]
        tv = np.nan
        if gate_ok:
            tv = sum(sc["w"][c] * obs_module(dsst[c], plans[c, dsg]) for c in sc["w"])
        rows.append(dict(
            id=IDS[dsg], figure_panel="Fig 6 (chronic disease: " + ("T2D" if dsg == "t2d" else "obesity") + ")",
            statistic="pooled module d* (case minus matched control), equivalence", direction=-1, value_14_stored=st0.D,
            value_14=s14["D"], reproduced=True, value_core=sc["D"],
            core_criterion=f"null (equivalence): core 90 % interval (D +- 1.645 SE, 500 matched-set bootstraps) inside +-{M.BOUND}",
            core_status=status, diff_14_minus_core=s14["D"] - sc["D"], diff_lo90=float(np.nanquantile(diff, .05)),
            diff_hi90=float(np.nanquantile(diff, .95)),
            diff_draws="500 matched-set bootstraps (A164 seeds), shared resamples, pooled D with each waveform's own weights",
            **pr, trio_gate_passes=gate_ok, trio_alone_value=tv, reading=C.reading(status, pr["trio_drop_below_q10"]),
            blind="core not blind (A335 / A329)",
            notes=(f"core 90% interval {sc['lo']:.4f} to {sc['hi']:.4f} (SE {sc['se']:.4f}, expected half-width {M.Z90 * sc['se']:.4f}; "
                   f"powered vs bound: {M.Z90 * sc['se'] <= M.BOUND}); 14 genes {s14['lo']:.4f} to {s14['hi']:.4f} (SE {s14['se']:.4f}); "
                   f"cohort weights core {({c: round(x, 3) for c, x in sc['w'].items()})}; drops: clock statistic only, observed pooled D "
                   f"with the core's cohort weights, direction -1 (loss of coupling); trio gate p_pairing/p_acrophase "
                   f"{ {c: tuple(round(x, 4) for x in v) for c, v in gp.items()} } in vent and fast"), script=SCRIPT))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws", "figure_panel", "statistic"]).round(4).T.to_string())
    for r in rows:
        print(r["id"], r["notes"])
    log(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
