"""A336 part 1, A306 (chronic disease: organ gain and departure unchanged; GTEx T2D, hypertension, smoking).

Original: `A306_disease_gain_departure/scripts/compute_a306.py` (imported, not edited) and the A233 matched procedure
(`A231_age_coupling/scripts/matched_label.py`). Same cohorts, matched sets, 1,000 shared within-set permutations (seeds
`A.seed_of("A306", condition, cohort, 0)`), the same Delta (module minus companion of per-tissue log z2 ratio = gain, log psi
ratio = departure, weights from the nulls), two-sided p, Holm over the readouts scored in the registration (the 14-gene
`scored` flags, kept). Only the waveform changes: core `cs.Dataset(..., waveform=core)`; the permutation nulls are
recomputed for the core (the 14-gene ones are the committed `output/draws`, their reproduction checks the pipeline). The
A306 companion is in-design (nerve, skins); no pseudo-clock sets.

Rows: one per condition x readout (6; smoking gain was descriptive, not powered). Nulls ("not detected" kept): core "holds"
if Holm p >= 0.05 and the 90 % interval of the core Delta (200 matched-set bootstraps: cases resampled jointly over
tissues with their sets, host refitted on the replicate, as `matched_label.boot_chunk`; seeds from rng [336, 1]) lies
inside +-MDE_80 (registered, 14-gene `score_tests.csv`); else "falls". Direction -1 (loss). The bootstrap draws are shared by
both waveforms; Delta on a bootstrap uses the weights of its waveform's null. Gene-drop null: Delta on 100 random drops (+ trio)
from the observed contrasts only, with the tissue / cohort weights FIXED at the core null's (no nulls per drop; the trio drop
then equals the core Delta). Trio alone: `Dataset.gate(n_perm=1000, seed=336)` p_pairing and p_acrophase < 0.05 in every
cohort of the condition; value = trio Delta with the core weights. GTEx: aggregates only; per-condition phenotypes stay in memory.

    python p1_a306.py [--check] [--jobs 4]
Output: ../output/part1/a306.csv (long run, about 1.5 h at 4 workers: use run.sh)
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
sys.path.insert(0, str(PAPER / "exploratory/A306_disease_gain_departure/scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_a306 as M6  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ML, A, B = M6.ML, M6.A, M6.B
ENTRY = "a306"
BASE = C.WF
OUTA = PAPER / "exploratory/A306_disease_gain_departure/output"
STORED = pd.read_csv(OUTA / "score_tests.csv").set_index(["condition", "readout"])
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a306.py"
READ = ("gain", "departure")
G = {}


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def ds_of(pb, anchors, wf, name):
    return cs.Dataset(pb, anchors=anchors, waveform=wf, name=name)


def obs_arrays(dss, plans):
    """cohort -> array (1, tissues, 2) of the observed log z2 / log psi ratios."""
    return {c: np.array([[M6.stat4(dss[c], A.labels(plans[c][t]), t) for t in plans[c]]]) for c in plans}


def deltas(arrs, mod, cmp_, nulls):
    """(gain Delta, departure Delta) arrays from replicate arrays and the weights of `nulls`."""
    return tuple(np.asarray(M6.delta(arrs, mod, cmp_, j, nulls)[0], float) for j in range(2))


def boot_chunk(pbs, anchors, plan, wfs, seeds):
    """`matched_label.boot_chunk` computing the two ratios for each waveform in wfs on the same resample:
    (replicates, tissues, waveforms, 2)."""
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
        rep = pbs.merge(m, left_on="donor", right_on="orig")
        rep["donor"] = rep["copy"]
        rep = rep.drop(columns=["orig", "copy"]).reset_index(drop=True)
        anc = [a for a in anchors if a in set(rep.compartment)]
        dsr = [cs.Dataset(rep, anchors=anc, waveform=w) for w in wfs]
        out.append([[M6.stat4(d, pd.concat([A.labels([s], copy=j) for s, j in sets[t]]), t) if sets[t] else [np.nan, np.nan]
                     for d in dsr] for t in plan])
    return np.array(out)


def drop_task(i):
    """Observed arrays of drop i for every (condition, cohort): {(cond, cohort): (1, tissues, 2)}."""
    warnings.filterwarnings("ignore")
    limit_threads(1)
    w = G["drops"][i][1]
    res = {}
    for name in G["conds"]:
        for c in G["plans"][name]:
            ds = ds_of(G["pbs"][name][c], G["anchors"][name][c], w, "drop")
            res[(name, c)] = np.array([[M6.stat4(ds, A.labels(G["plans"][name][c][t]), t) for t in G["plans"][name][c]]])
    return i, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--conditions", nargs="*", default=list(M6.MODULES))
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    warnings.filterwarnings("ignore", message="below the paper settings")
    OUT14 = OUTA
    cohorts = ML.load_cohorts(None, 0)
    core, trio = C.core(BASE), C.trio(BASE)
    conds = a.conditions
    D = {}
    # ---- stage A: 14 genes, cached nulls (reproduction) --------------------------------------------------------------
    for name in conds:
        d = M6.design_of(name)
        tissues = d.modules["module"] + d.modules["companion"]
        plans, dss, anchors, pbs, _, _ = ML.build(d, cohorts, quick=chk, tissues=tissues, tag="gtex_a306")
        mod, cmp_ = ML.selection(plans, d.modules["module"]), ML.selection(plans, d.modules["companion"])
        nulls14 = {c: np.load(OUT14 / ML.DRAW_SUB / f"null_{name}_{c}.npy") for c in plans}
        D[name] = dict(d=d, plans=plans, dss=dss, anchors=anchors, pbs=pbs, mod=mod, cmp=cmp_, nulls14=nulls14)
        obs = obs_arrays(dss, plans)
        D[name]["obs14"] = obs
        D[name]["delta14"] = [float(np.asarray(M6.delta(obs, mod, cmp_, j, nulls14)[0]).item()) for j in range(2)]
        D[name]["nd14"] = [M6.delta(nulls14, mod, cmp_, j, nulls14)[0] for j in range(2)]
        log(name, "14-gene observed done")
    if not chk:
        rep = {}
        for name in conds:
            for j, ro in enumerate(READ):
                st = STORED.loc[(name, ro)]
                rep[(name, ro)] = bool(np.isclose(D[name]["delta14"][j], st.delta, atol=1e-6) and
                                       np.isclose(A.p2(D[name]["delta14"][j], D[name]["nd14"][j]), st.p, atol=1e-9))
        log("reproduction", rep)
        if not all(rep.values()):
            C.write_rows(ENTRY, [dict(id=f"A306_{n}_{r}", figure_panel="Fig 5/6", statistic=f"Delta {r}, {n}", direction=-1,
                                      value_14_stored=STORED.loc[(n, r), "delta"], value_14=D[n]["delta14"][READ.index(r)],
                                      reproduced=rep[(n, r)], notes=f"not reproduced: {rep}", script=SCRIPT)
                                 for n in conds for r in READ])
            raise SystemExit(f"A306 not reproduced: {rep}")
    # ---- stage B: core nulls ----------------------------------------------------------------------------------------
    n_perm = 20 if chk else 1000
    for name in conds:
        x = D[name]
        x["dssc"] = {c: ds_of(x["pbs"][c], x["anchors"][c], core, f"core_{name}_{c}") for c in x["plans"]}
        x["nullsc"] = {c: A.run_chunks(M6.perm_chunk, n_perm, a.jobs, A.seed_of("A306", name, c, 0), f"core null {name} {c}",
                                      x["dssc"][c], x["plans"][c]) for c in x["plans"]}
        x["obsc"] = obs_arrays(x["dssc"], x["plans"])
        x["deltac"] = [float(np.asarray(M6.delta(x["obsc"], x["mod"], x["cmp"], j, x["nullsc"])[0]).item()) for j in range(2)]
        x["ndc"] = [M6.delta(x["nullsc"], x["mod"], x["cmp"], j, x["nullsc"])[0] for j in range(2)]
        log(name, "core null done")
    # ---- bootstraps (shared) ---------------------------------------------------------------------------------------
    nb = 4 if chk else 200
    master = np.random.default_rng([336, 1])
    for name in conds:
        x = D[name]
        x["boot"] = {}
        for c in x["plans"]:
            seeds = [int(s) for s in master.integers(0, 2 ** 31 - 1, nb)]
            size = max(1, int(np.ceil(nb / (a.jobs * 3))))
            chunks = [seeds[i:i + size] for i in range(0, nb, size)]
            with Pool(a.jobs) as p:
                res = p.starmap(boot_chunk, [(x["pbs"][c], x["anchors"][c], x["plans"][c], [BASE, core], ch) for ch in chunks])
            x["boot"][c] = np.vstack(res)  # (nb, tissues, 2 waveforms, 2 ratios)
        b14 = {c: v[:, :, 0, :] for c, v in x["boot"].items()}
        bc = {c: v[:, :, 1, :] for c, v in x["boot"].items()}
        x["b14"] = deltas(b14, x["mod"], x["cmp"], x["nulls14"])
        x["bc"] = deltas(bc, x["mod"], x["cmp"], x["nullsc"])
        log(name, "bootstrap done")
    # ---- gene drops --------------------------------------------------------------------------------------------------
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [q for q in dl if q[0] == tl] + dl[:2]
    G.update(drops=dl, conds=conds, plans={n: D[n]["plans"] for n in conds}, pbs={n: D[n]["pbs"] for n in conds},
             anchors={n: D[n]["anchors"] for n in conds})
    with Pool(a.jobs) as p:
        dres = dict(p.imap_unordered(drop_task, range(len(dl))))
    dv = {(n, r): {} for n in conds for r in READ}
    for i, (lab, _) in enumerate(dl):
        for name in conds:
            x = D[name]
            arrs = {c: dres[i][(name, c)] for c in x["plans"]}
            dd = deltas(arrs, x["mod"], x["cmp"], x["nullsc"])
            for j, r in enumerate(READ):
                dv[(name, r)][lab] = float(dd[j][0])
    for name in conds:
        for j, r in enumerate(READ):
            assert np.isclose(dv[(name, r)][tl], D[name]["deltac"][j], equal_nan=True)
    log("drops done")
    wfs = dict(dl)
    # ---- trio alone --------------------------------------------------------------------------------------------------
    for name in conds:
        x = D[name]
        gp = []
        for c in x["plans"]:
            g = ds_of(x["pbs"][c], x["anchors"][c], trio, "trio").gate(n_perm=50 if chk else 1000, seed=336)
            gp.append((float(g.p_pairing), float(g.p_acrophase)))
        x["gate_p"] = gp
        x["gate"] = bool(all(p1 < 0.05 and p2 < 0.05 for p1, p2 in gp))
        if x["gate"]:
            dst = {c: ds_of(x["pbs"][c], x["anchors"][c], trio, "trio") for c in x["plans"]}
            x["trio_delta"] = [float(v[0]) for v in deltas(obs_arrays(dst, x["plans"]), x["mod"], x["cmp"], x["nullsc"])]
    if chk:
        print("check ok")
        return
    # ---- rows --------------------------------------------------------------------------------------------------------
    scored14 = {k: bool(v.scored) for k, v in STORED.iterrows()}
    pvals = {(n, r): A.p2(D[n]["deltac"][j], D[n]["ndc"][j]) for n in conds for j, r in enumerate(READ)}
    fam = [k for k in pvals if scored14[k]]
    hp = dict(zip(fam, ML.holm(np.array([pvals[k] for k in fam]))))
    rows = []
    for name in conds:
        x = D[name]
        for j, r in enumerate(READ):
            k = (name, r)
            st = STORED.loc[k]
            dc, d14 = x["deltac"][j], x["delta14"][j]
            lo, hi = np.nanquantile(x["bc"][j], [0.05, 0.95])
            diff = x["b14"][j] - x["bc"][j]
            mde = float(st.mde_80)
            detected = bool(k in hp and hp[k] < 0.05)
            holds = (not detected) and (-mde <= lo) and (hi <= mde)
            status = "holds" if holds else "falls"
            pr = C.drop_percentiles(dv[k], tl, wfs, direction=-1)
            sd = float(np.nanstd(x["ndc"][j], ddof=1))
            rows.append(dict(
                id=f"A306_{name}_{r}", figure_panel="Fig 5/6 (chronic disease: gain and departure unchanged)",
                statistic=f"Delta = module minus companion of per-tissue log {'z2' if r == 'gain' else 'psi'} ratio, {name}",
                direction=-1, value_14_stored=st.delta, value_14=d14, reproduced=True, value_core=dc,
                core_criterion=(f"null (not detected): Holm p >= 0.05 and 90 % bootstrap interval inside +-MDE_80 {mde:.3f} (registered, 14 genes)"
                                + ("; DESCRIPTIVE in the registration (not powered), same rule applied" if not scored14[k] else "")),
                core_status=status, diff_14_minus_core=d14 - dc, diff_lo90=float(np.nanquantile(diff, .05)),
                diff_hi90=float(np.nanquantile(diff, .95)), diff_draws="200 matched-set bootstraps (seeds from rng [336, 1]), shared",
                **pr, trio_gate_passes=x["gate"], trio_alone_value=x["trio_delta"][j] if x["gate"] else np.nan,
                reading=C.reading(status, pr["trio_drop_below_q10"]), blind="core blind",
                notes=(f"core p {pvals[k]:.4f}" + (f", Holm {hp[k]:.4f}" if k in hp else ", not in Holm family") +
                       f" (14 genes p {st.p:.4f}); core null SD {sd:.4f}, own MDE_80 {2.49 * sd:.3f}; core 90 % bootstrap {lo:.3f} to {hi:.3f}; "
                       f"trio gate p_pairing / p_acrophase per cohort {[(round(a_, 4), round(b_, 4)) for a_, b_ in x['gate_p']]}; "
                       f"drop weights fixed at the core null's"), script=SCRIPT))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws", "figure_panel", "statistic"]).round(4).to_string())
    for r in rows:
        print(r["id"], r["notes"])


if __name__ == "__main__":
    main()
