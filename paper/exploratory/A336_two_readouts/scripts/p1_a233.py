"""A336 part 1 wrapper: A233 H233.c.1 (body-wide pooled d* equivalence) and H233.c.2-3 (body-wide reordering Q; module d* and
module offset) for hypertension, IHD, COPD, depression and smoking, on the 14 genes, the core loop, the trio contribution,
the gene-drop null and the trio alone. Rows per condition: `A233_H233.c.1:<c>`, `A233_H233.c.2-3:<c>:Q`,
`:module_d_star`, `:module_offset`.

Originals (imported, not edited): A233 `compute_a233.py` (`design_of`, modules) and the A231 machinery `matched_label.py`
(matching, `perm_chunk` / `observed` / `pooled_tests` / `q_test`, A208 pooling). Matching and labels are waveform-free (same
sets). 14 genes: recomputed from the committed null and bootstrap draws (`A233/output/draws`, 1,000 permutations, 200
bootstraps) with the observed contrasts recomputed, and compared with score_tests.csv (1e-6) before any core value.
Core: `cs.Dataset(..., waveform=core)` for the observed contrasts, the permutation null and the matched-set bootstrap
(`boot_chunk` copied with a waveform argument, same seeds `seed_of("A231A233", name, "", cohort, 0)` (+ 10**6 for bootstraps)).

Deviations (registered rules otherwise unchanged): (1) the core permutation null uses the FIRST 200 of the 1,000 permutation
seeds (the same label permutations as the first 200 14-gene replicates), not 1,000, for run time (the full core null is the
3 h power stage); p-values therefore have resolution 1/201 and Holm is over the 5 conditions; (2) the pseudo-clock companion
sets of the original score stage are NOT recomputed on the core: none of these registered criteria uses them.
Criteria on the core: c.1: the 90 % interval (estimate +- 1.645 x max(bootstrap SE, null SD)) inside +-0.15 (equivalent;
reported with the core MDE_80 <= 0.15); c.2: Q one-sided permutation p, Holm over the conditions, >= 0.05 (the null is kept);
c.3: module d* and module offset two-sided p, Holm over the conditions with a module, >= 0.05. A statistic whose null is kept
"holds"; one that crosses reads "falls" (reading() then uses the gene-drop flag).
Trio contribution: the statistic's own 200 matched-set bootstrap draws, pooled with each waveform's own weights (Q: Q of the
per-tissue bootstrap offsets with each waveform's null weights), 90 % interval of 14 minus core. Gene drops: 100 random +
trio, observed statistic only (core tissue weights fixed, the trio drop equals the core), magnitude |estimate| with
direction -1 (a larger departure from the null is the stronger effect). Trio alone: Dataset gate (n_perm 1000, seed 336) in
the vent and fast cohorts on the first condition's tissue set, both p < 0.05, then the observed statistic.
GTEx protected phenotypes: aggregates only.

    python p1_a233.py [--check [--check-donors 60]] [--jobs 4]     (long: paper/common/run.sh, RUN_CAP_GB=40)
Output: ../output/part1/a233.csv; ../output/part1_a233/draws/*.npy (core draws, resumable)
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
sys.path.insert(0, str(PAPER / "exploratory/A233_chronic_conditions/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A231_age_coupling/scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_a233 as M233  # noqa: E402
import matched_label as ML  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

A, B = ML.A, ML.B
ENTRY = "a233"
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a233.py"
OUTA = PAPER / "exploratory/A233_chronic_conditions/output"
OUT = HERE.parents[1] / "output" / "part1_a233"
BASE = C.WF
CONDS = list(M233.COLUMN)
N_PERM_CORE = 200
G = {}


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def ds_of(pb, anchors, wf, name=""):
    return cs.Dataset(pb, anchors=anchors, waveform=wf, name=name)


def boot_chunk(pb, anchors, plan, wf, seeds):
    """`matched_label.boot_chunk` with an explicit waveform."""
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


def pooled_parts(obs, nulls, boots, sel, j):
    """(null, tw, cw, observed pooled, bootstrap pooled draws) of readout j on the tissues `sel`."""
    sn, sb = ML.sub(nulls, sel), ML.sub(boots, sel)
    null, tw, cw, _ = B.pooled_null(sn, j)
    D, _ = B.pool_estimate({c: obs[c][sel[c], j] for c in cw}, tw, cw)
    bt = sum(cw[c] * B.wmean(sb[c][:, :, j], tw[c]) for c in cw)
    return dict(null=null, tw=tw, cw=cw, D=D, bt=bt)


def score(name, d, plans, obs, nulls, boots):
    """Everything the criteria and the trio contribution need for one condition and one waveform."""
    all_sel = ML.selection(plans, lambda t: True)
    mod_sel = ML.selection(plans, d.modules["module"])
    out = dict(all_sel=all_sel, mod_sel=mod_sel, bw=ML.pooled_tests(obs, nulls, boots, all_sel, "bodywide", name).set_index(["readout", "set"]),
               Q=ML.q_test(obs, nulls, all_sel))
    out["mod"] = ML.pooled_tests(obs, nulls, boots, mod_sel, "module", name).set_index(["readout", "set"]) if mod_sel else None
    out["parts"] = {("bw", j): pooled_parts(obs, nulls, boots, all_sel, j) for j in (0, 1)}
    if mod_sel:
        out["parts"].update({("mod", j): pooled_parts(obs, nulls, boots, mod_sel, j) for j in (0, 1)})
    Xb = np.hstack([boots[c][:, all_sel[c], 1] for c in all_sel])
    out["Qboot"] = ML.q_of(Xb, out["Q"]["w"])
    return out


def drop_task(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    w = G["drops"][i][1]
    res = {}
    for name, (plans, pbs, anchors) in G["conds"].items():
        dss = {c: ds_of(pbs[c], anchors[c], w) for c in plans}
        obs, _ = ML.observed(dss, plans)
        res[name] = obs
    return i, res


def pooled_obs(obs, sc, kind):
    """Drop statistics with the core's weights: body-wide d*, Q, module d*, module offset."""
    out = {}
    pb = sc["parts"]
    cw0, tw0 = pb[("bw", 0)]["cw"], pb[("bw", 0)]["tw"]
    out["bw_d"] = B.pool_estimate({c: obs[c][sc["all_sel"][c], 0] for c in cw0}, tw0, cw0)[0]
    xo = np.concatenate([obs[c][sc["all_sel"][c], 1] for c in sc["all_sel"]])[None]
    out["Q"] = float(ML.q_of(xo, sc["Q"]["w"])[0])
    if sc["mod_sel"]:
        for k, j in (("mod_d", 0), ("mod_off", 1)):
            p = pb[("mod", j)]
            out[k] = B.pool_estimate({c: obs[c][sc["mod_sel"][c], j] for c in p["cw"]}, p["tw"], p["cw"])[0]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check-donors", type=int, default=60, help="donors per cohort in --check")
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    warnings.filterwarnings("ignore", message="below the paper settings")
    t0 = time.time()
    core, trio = C.core(BASE), C.trio(BASE)
    cohorts = ML.load_cohorts(a.check_donors if chk else None, 0)
    tissues = ["Lung", "Artery - Aorta", "Heart - Left Ventricle", "Whole Blood", "Brain - Cortex"] if chk else None
    n_perm, n_boot = (4, 4) if chk else (N_PERM_CORE, 200)
    (OUT / "draws").mkdir(parents=True, exist_ok=True)
    stored = pd.read_csv(OUTA / "score_tests.csv").set_index("condition")
    built, S14, bad = {}, {}, []
    names = CONDS[:2] if chk else CONDS
    # ---- 14 genes from the committed draws, reproduced first -------------------------------------------------------
    for name in names:
        d = M233.design_of(name)
        plans, dss, anchors, pbs, _, _ = ML.build(d, cohorts, quick=chk, tissues=tissues, tag="gtex_cond")
        if not plans:
            continue
        built[name] = (d, plans, anchors, pbs, dss)
        if chk:
            nulls = {c: np.random.default_rng(0).standard_normal((n_perm, len(plans[c]), 2)) for c in plans}
            boots = {c: np.random.default_rng(1).standard_normal((n_boot, len(plans[c]), 2)) for c in plans}
        else:
            nulls = ML.nulls_of(name, plans, dss, OUTA, False, None, a.jobs, 0)
            boots = ML.boots_of(name, pbs, anchors, plans, OUTA, False, None, a.jobs, 0)
        obs, _ = ML.observed(dss, plans)
        S14[name] = dict(sc=score(name, d, plans, obs, nulls, boots), nulls=nulls, boots=boots)
    if not chk:
        names = list(built)
        hq = ML.holm([S14[n]["sc"]["Q"]["p_perm"] for n in names])
        mods = [n for n in names if S14[n]["sc"]["mod_sel"]]
        hd = dict(zip(mods, ML.holm([S14[n]["sc"]["mod"].loc[("d_star", "pooled"), "p_perm"] for n in mods])))
        ho = dict(zip(mods, ML.holm([S14[n]["sc"]["mod"].loc[("offset_h", "pooled"), "p_perm"] for n in mods])))
        for i, n in enumerate(names):
            s, st = S14[n]["sc"], stored.loc[n]
            bw = s["bw"].loc[("d_star", "pooled")]
            pairs = [(bw.estimate, st.bodywide_d_star), (bw.lo90, st.bodywide_lo90), (bw.hi90, st.bodywide_hi90), (bw.mde_80, st.bodywide_mde_80),
                     (s["Q"]["Q"], st.Q), (s["Q"]["p_perm"], st.Q_p), (hq[i], st.Q_p_holm)]
            if s["mod_sel"]:
                md, mo = s["mod"].loc[("d_star", "pooled")], s["mod"].loc[("offset_h", "pooled")]
                pairs += [(md.estimate, st.module_d_star), (md.p_perm, st.module_d_star_p), (hd[n], st.module_d_star_p_holm),
                          (mo.estimate, st.module_offset_h), (mo.p_perm, st.module_offset_p), (ho[n], st.module_offset_p_holm)]
            bad += [(n, float(x), float(y)) for x, y in pairs if not np.isclose(x, y, atol=1e-6, rtol=0)]
        log("14-gene reproduction:", "OK" if not bad else bad)
        if bad:
            rows = [dict(id="A233_H233.c.1:" + n, figure_panel="Fig 6", reproduced=False, notes=f"14-gene values not reproduced: {[b for b in bad if b[0] == n]}", script=SCRIPT)
                    for n in names]
            C.write_rows(ENTRY, rows)
            raise SystemExit("A233 not reproduced")
    # ---- core: observed, null (first 200 permutations), bootstraps ---------------------------------------------------
    SC = {}
    for name, (d, plans, anchors, pbs, dss) in built.items():
        dssc = {c: ds_of(pbs[c], anchors[c], core, f"core_{name}_{c}") for c in plans}
        nulls = {c: A.cached(OUT / "draws" / f"null_{name}_{c}.npy", chk, B.perm_chunk, n_perm, a.jobs,
                             A.seed_of("A231A233", name, "", c, 0), f"core null {name} {c}", dssc[c], plans[c]) for c in plans}
        boots = {c: A.cached(OUT / "draws" / f"boot_{name}_{c}.npy", chk, boot_chunk, n_boot, a.jobs,
                             A.seed_of("A231A233", name, "", c, 0) + 10 ** 6, f"core boot {name} {c}", pbs[c], anchors[c], plans[c], core)
                 for c in plans}
        obs, _ = ML.observed(dssc, plans)
        SC[name] = score(name, d, plans, obs, nulls, boots)
        log(f"{name} core done ({(time.time() - t0) / 60:.0f} min)")
    # ---- gene drops -----------------------------------------------------------------------------------------------------
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:1]
    G.update(drops=dl, conds={n: (b[1], b[3], b[2]) for n, b in built.items()})
    with Pool(a.jobs) as p:
        dres = dict(p.imap_unordered(drop_task, range(len(dl))))
    dv = {n: {k: {} for k in ("bw_d", "Q", "mod_d", "mod_off")} for n in built}
    for i, (lab, _) in enumerate(dl):
        for n in built:
            for k, v in pooled_obs(dres[i][n], SC[n], None).items():
                dv[n][k][lab] = v
    for n in built:
        sc = SC[n]
        assert np.isclose(dv[n]["bw_d"][tl], sc["parts"][("bw", 0)]["D"], equal_nan=True)
        assert np.isclose(dv[n]["Q"][tl], sc["Q"]["Q"], equal_nan=True)
    log("drops done")
    # ---- trio alone -----------------------------------------------------------------------------------------------------
    first = next(iter(built))
    _, plans1, anchors1, pbs1, _ = built[first]
    gates, gate_ok = {}, True
    for c in plans1:
        try:
            g = ds_of(pbs1[c], anchors1[c], trio, f"trio_{c}").gate(n_perm=50 if chk else 1000, seed=336)
            gates[c] = (float(g.p_pairing), float(g.p_acrophase))
        except Exception as e:  # noqa: BLE001
            gates[c] = (np.nan, np.nan)
            log("gate failed", c, type(e).__name__)
    gate_ok = all(np.isfinite(p1) and p1 < 0.05 and p2 < 0.05 for p1, p2 in gates.values())
    trio_vals = {}
    if gate_ok:
        for n, (d, plans, anchors, pbs, _) in built.items():
            dsst = {c: ds_of(pbs[c], anchors[c], trio) for c in plans}
            o, _ = ML.observed(dsst, plans)
            trio_vals[n] = pooled_obs(o, SC[n], None)
    if chk:
        print("check ok")
        return
    # ---- criteria and rows -------------------------------------------------------------------------------------------
    names = list(built)
    hq = dict(zip(names, ML.holm([SC[n]["Q"]["p_perm"] for n in names])))
    mods = [n for n in names if SC[n]["mod_sel"]]
    hd = dict(zip(mods, ML.holm([SC[n]["mod"].loc[("d_star", "pooled"), "p_perm"] for n in mods])))
    ho = dict(zip(mods, ML.holm([SC[n]["mod"].loc[("offset_h", "pooled"), "p_perm"] for n in mods])))
    dwf = dict(dl)
    gtxt = f"trio gate (vent, fast) p_pairing/p_acrophase { {c: tuple(round(x, 4) for x in v) for c, v in gates.items()} }"
    rows = []

    def row(sid, panel, stat, key, name, v14, vstored, vc, crit, status, d14, dcore, notes, tv):
        pr = C.drop_percentiles({k: abs(v) if np.isfinite(v) else np.nan for k, v in dv[name][key].items()}, tl, dwf, direction=-1)
        diff = d14 - dcore
        return dict(id=sid, figure_panel=panel, statistic=stat, direction=-1, value_14_stored=vstored, value_14=v14, reproduced=True, value_core=vc,
                    core_criterion=crit, core_status=status, diff_14_minus_core=v14 - vc, diff_lo90=float(np.nanquantile(diff, .05)),
                    diff_hi90=float(np.nanquantile(diff, .95)),
                    diff_draws="200 matched-set bootstraps (A233 seeds), pooled with each waveform's own weights", **pr,
                    trio_gate_passes=gate_ok, trio_alone_value=(abs(tv) if gate_ok and tv is not None else np.nan),
                    reading=C.reading(status, pr["trio_drop_below_q10"]), blind="core blind",
                    notes=f"{notes}; core permutation null 200 (first 200 shared seeds), pseudo-clock companions not recomputed; drops ranked by |estimate|, direction -1; {gtxt}",
                    script=SCRIPT)

    for n in names:
        s14, sc, st = S14[n]["sc"], SC[n], stored.loc[n]
        p14, pc = s14["parts"], sc["parts"]
        bw = sc["bw"].loc[("d_star", "pooled")]
        bw14 = s14["bw"].loc[("d_star", "pooled")]
        eq = bool(bw.equivalent)
        rows.append(row(f"A233_H233.c.1:{n}", "Fig 6 (chronic disease: body-wide coupling)", "body-wide pooled d* equivalence within +-0.15", "bw_d", n,
                        float(bw14.estimate), st.bodywide_d_star, float(bw.estimate), "core 90% interval inside +-0.15 (equivalent), reported with MDE_80 <= 0.15",
                        "holds" if eq else "falls", p14[("bw", 0)]["bt"], pc[("bw", 0)]["bt"],
                        f"core 90% [{bw.lo90:.4f}, {bw.hi90:.4f}] (boot SE {bw.boot_se:.4f}, null SD {bw.null_sd:.4f}, MDE_80 {bw.mde_80:.4f}); 14 genes [{bw14.lo90:.4f}, {bw14.hi90:.4f}]",
                        trio_vals.get(n, {}).get("bw_d")))
        qc = sc["Q"]
        rows.append(row(f"A233_H233.c.2-3:{n}:Q", "Fig 6 (chronic disease: reordering)", "body-wide reordering Q (organ offsets), one-sided permutation p", "Q", n,
                        float(s14["Q"]["Q"]), st.Q, float(qc["Q"]), "core Q p (Holm over 5) >= 0.05 (the null is kept)",
                        "holds" if hq[n] >= 0.05 else "falls", s14["Qboot"], sc["Qboot"],
                        f"core Q p {qc['p_perm']:.4f}, Holm {hq[n]:.4f} (stored 14-gene p {st.Q_p:.4f}, Holm {st.Q_p_holm:.4f}); Q difference draws are Q of bootstrapped offsets",
                        trio_vals.get(n, {}).get("Q")))
        if sc["mod_sel"]:
            for key, rd, sid, h, lab, vst in (("mod_d", "d_star", "module_d_star", hd, "module d*", st.module_d_star), ("mod_off", "offset_h", "module_offset", ho, "module offset (h)", st.module_offset_h)):
                m, m14 = sc["mod"].loc[(rd, "pooled")], s14["mod"].loc[(rd, "pooled")]
                j = 0 if key == "mod_d" else 1
                rows.append(row(f"A233_H233.c.2-3:{n}:{sid}", "Fig 6 (chronic disease: local shift)", lab + " two-sided permutation p", key, n, float(m14.estimate), vst,
                                float(m.estimate), "core p (Holm over conditions with a module) >= 0.05 (the null is kept)",
                                "holds" if h[n] >= 0.05 else "falls", p14[("mod", j)]["bt"], pc[("mod", j)]["bt"],
                                f"core p {m.p_perm:.4f}, Holm {h[n]:.4f}; MDE_80 core {m.mde_80:.4f}", trio_vals.get(n, {}).get(key)))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws", "figure_panel", "statistic", "blind", "direction", "reproduced", "n_undefined"]).round(4).to_string())
    for r in rows:
        print(r["id"], r["notes"][:260])
    log(f"done {(time.time() - t0) / 60:.0f} min")


if __name__ == "__main__":
    main()
