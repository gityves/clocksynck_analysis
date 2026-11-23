"""A336 part 1 wrapper: tumour-normal immune pairing. A227 H227.1 (myeloid), H227.2 (T/NK) and A271 H271.1 (myeloid, T/NK and
CD8-negative T of A228 H228.4 against the responsive family) on the 14 genes, the core loop, the trio contribution, the
gene-drop null and the trio alone. Entries written: a227 (H227.1-2) and a271_h1 (H271.1: one row per component plus an
aggregate row).

Originals (imported, not edited): A227 `compute_a227.py` (tables, `units_of`, `raw_stat`, `perm_null_raw`, `pooled_pseudo`,
`holm`) and A228 `compute_a228.py` (CD8-negative T units and genome-wide subset tables). Same cohorts, patients, seeds
(`seed_of`, permutation null seed + [1], 10,000 within-cohort permutations; pseudo-clock sets `pseudo_waveforms(lm, 200,
default_rng(seed + [44]), cs.gtex(), family)` per cohort and kind); only the waveform changes. Companions: the 200 matched
and 200 responsive 14-gene draws restricted to the core positions. 14 genes: observed pooled statistic, null mean / SD and
p_pairing recomputed against the stored tests tables (1e-6), and the first 10 matched and responsive companion statistics per
cohort against the stored pseudo.csv; the remaining 14-gene companions are NOT refitted (stored q95 used for reference).
Criteria on the core: H227.x: pairing p < 0.05 (Holm over myeloid, T/NK) and statistic above the matched core q95; H271.1:
pooled statistic above the responsive core q95 and responsive p Holm (over the 3 components) < 0.05, every component.
Trio contribution: 200 patient bootstraps within cohort (rng [336, 1]), same resamples for both waveforms, pooled with the
original patient-count weights. Gene drops: 100 random + trio, pooled statistic only. Trio alone: Dataset gate (cohorts as
strata, the two compartments of the pair) and the pooled statistic if it passes. Memory heavy: 2 workers.

    python p1_a227.py [--check] [--jobs 2]        (long: paper/common/run.sh, RUN_CAP_GB=40)
Output: ../output/part1/a227.csv, a271_h1.csv
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
for p in ("A227_tme_compartment_pairing", "A228_immune_body_time", "A149_paired_tumour_normal", "A151_fibroblast_pairing",
          "A226_stroma_pairing_confirm"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a227 as A  # noqa: E402
import compute_a228 as A8  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
import clocksync as cs  # noqa: E402

SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a227.py"
BASE = C.WF
CORE, TRIO = C.core(BASE), C.trio(BASE)
EXP = PAPER / "exploratory"
KINDS = {"myeloid": ("A227", A, "myeloid"), "T_NK": ("A227", A, "T_NK"), "CD8neg_T": ("A228", A8, "CD8neg_T")}
FAMS = ("matched", "responsive")
N_CHECK = 10


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def ab(kind):
    return A.PAIRS[kind] if KINDS[kind][0] == "A227" else A8.TESTS[kind][3:]


def study_task(study, kinds, n_perm, n_pseudo, drop_list):
    """Everything one cohort contributes for each kind it has: states on every waveform, nulls, companions, drops."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = {}
    for kind in kinds:
        src, M, _ = KINDS[kind]
        pb = M.table(study) if src == "A227" else M.pb_of(study, kind)
        if pb is None:
            continue
        a, b = ab(kind)
        if a not in set(pb.compartment) or b not in set(pb.compartment):
            continue
        u = L.pair_states(pb, a, b)
        if not u:
            continue
        pts = u["pts"]
        r = dict(n=u["n"], z14=(u["za"], u["zb"]))
        uc = L.pair_states(pb, a, b, wf=CORE, pts=pts)
        r["zc"] = (uc["za"], uc["zb"])
        r["null14"] = A.perm_null_raw(u, n_perm, np.random.default_rng(M.seed_of(kind, study) + [1]))
        r["nullc"] = A.perm_null_raw(uc, n_perm, np.random.default_rng(M.seed_of(kind, study) + [1]))
        # gene drops and trio (pooled statistic later)
        drops = {}
        for lab, w in drop_list + [("TRIO_ALONE", TRIO)]:
            p = L.pair_states(pb, a, b, wf=w, pts=pts)
            drops[lab] = A.raw_stat(p["za"], p["zb"]) if p else np.nan
        r["drops"] = drops
        # gate rows (clock genes only, stratum = cohort)
        rows = []
        for comp in (a, b):
            x = L._rows(pb, comp).loc[pts].reset_index()
            x = x[["donor", "compartment", "stratum", "depth"] + [g for g in BASE.genes]].copy()
            x["donor"] = study + "|" + x["donor"].astype(str)
            x["stratum"] = study
            x["n_cells"] = 100
            rows.append(x)
        r["gate_rows"] = pd.concat(rows, ignore_index=True)
        # companions
        if src == "A227":
            G = A.genome_pair(study, (a, b), pts)
        else:
            G = M.subset_table(study)
            G = G[G.donor.isin(pts) & G.compartment.isin([a, b])].reset_index(drop=True)
        ref = G[G.compartment == b]
        genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
        frac = ref[genes].div(ref.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        pse = {}
        for fam in FAMS:
            ws = pseudo_waveforms(lm, n_pseudo, np.random.default_rng(M.seed_of(kind, study) + [44]), BASE, fam)
            st14 = []
            for w in ws[:N_CHECK]:
                p = L.pair_states(G, a, b, wf=w, pts=pts)
                st14.append(A.raw_stat(p["za"], p["zb"]) if p else np.nan)
            stc = []
            for w in ws:
                p = L.pair_states(G, a, b, wf=C.restrict(w, BASE, CORE), pts=pts)
                stc.append(A.raw_stat(p["za"], p["zb"]) if p else np.nan)
            pse[fam] = (np.array(st14), np.array(stc))
        r["pse"] = pse
        out[kind] = r
    return study, out


def stored(kind):
    """Stored 14-gene tests row and pseudo columns (matched and responsive) of a kind."""
    if KINDS[kind][0] == "A227":
        d = EXP / "A227_tme_compartment_pairing/output"
    else:
        d = EXP / "A228_immune_body_time/output"
    t = {f: pd.read_csv((d if f == "matched" else d / "responsive") / "tests.csv").set_index("kind").loc[kind] for f in FAMS}
    ps = {f: pd.read_csv((d if f == "matched" else d / "responsive") / "pseudo.csv") for f in FAMS}
    return t, ps


def pool(res, kind, which):
    """Pooled (weights = patients) statistic of waveform `which` (z14 / zc) over the cohorts with the kind."""
    sel = [s for s in res if kind in res[s]]
    w = np.array([res[s][kind]["n"] for s in sel], float)
    st = [A.raw_stat(*res[s][kind][which]) for s in sel]
    return sel, w, float(np.sum(w * st) / w.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    n_perm, n_pseudo = (200, 12) if chk else (10_000, 200)
    cohorts = A.cohort_list(chk)
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    t0 = time.time()
    kinds = list(KINDS)
    res = dict(Parallel(n_jobs=a.jobs, return_as="generator_unordered")(
        delayed(study_task)(s, kinds, n_perm, n_pseudo, dl) for s in cohorts))
    log(f"cohorts done {(time.time() - t0) / 60:.1f} min")
    # ---- 14 genes, reproduced first -------------------------------------------------------------------------------
    st, bad, rep = {}, [], {}
    for kind in kinds:
        t, ps = stored(kind)
        sel, w, obs14 = pool(res, kind, "z14")
        nul = A.pool_null(w, [res[s][kind]["null14"] for s in sel])
        p = float((1 + np.sum(nul >= obs14)) / (1 + len(nul)))
        tm = t["matched"]
        rep[kind] = dict(obs=obs14, null_mean=float(nul.mean()), null_sd=float(nul.std()), p=p, sel=sel, w=w, nul14=nul)
        if not chk:
            checks = [("statistic", obs14, tm.statistic), ("null_mean", nul.mean(), tm.null_mean), ("null_sd", nul.std(), tm.null_sd), ("p_pairing", p, tm.p_pairing)]
            for fam in FAMS:
                got = np.concatenate([res[s][kind]["pse"][fam][0] for s in sel])
                want = np.concatenate([ps[fam][f"{s}|{kind}"].to_numpy()[:N_CHECK] for s in sel])
                if not np.allclose(got, want, atol=1e-6, rtol=0, equal_nan=True):
                    checks.append((f"companions_{fam}", float(np.nanmax(np.abs(got - want))), 0.0))
            bad += [(kind, n, float(x), float(y)) for n, x, y in checks if not np.isclose(x, y, atol=1e-6, rtol=0)]
    log("14-gene reproduction:", "OK" if not bad else bad)
    if bad:
        rows = [dict(id=i, figure_panel="Fig 6", value_14_stored=stored(k)[0]["matched"].statistic, value_14=rep[k]["obs"], reproduced=False,
                     notes=f"14-gene values not reproduced (tolerance 1e-6): {[b for b in bad if b[0] == k]}", script=SCRIPT)
                for i, k in (("A227_H227.1", "myeloid"), ("A227_H227.2", "T_NK"))]
        C.write_rows("a227", rows)
        C.write_rows("a271_h1", [dict(id="A271_H271.1", figure_panel="Fig 6", reproduced=False, notes=f"14-gene values not reproduced: {bad}", script=SCRIPT)])
        raise SystemExit("not reproduced")
    # ---- core: pooled statistic, null, p, companions ---------------------------------------------------------------
    core = {}
    for kind in kinds:
        sel, w = rep[kind]["sel"], rep[kind]["w"]
        _, _, obsc = pool(res, kind, "zc")
        nulc = A.pool_null(w, [res[s][kind]["nullc"] for s in sel])
        d = dict(obs=obsc, p_pair=float((1 + np.sum(nulc >= obsc)) / (1 + len(nulc))), nul_sd=float(nulc.std()), nul_mean=float(nulc.mean()))
        for fam in FAMS:
            pse = A.pooled_pseudo(np.array([res[s][kind]["pse"][fam][1] for s in sel]), w)
            d[fam] = dict(median=float(np.median(pse)), q95=float(np.quantile(pse, .95)), n=len(pse),
                          p=float((1 + np.sum(pse >= obsc)) / (1 + len(pse))))
        core[kind] = d
    # ---- trio contribution: patient bootstrap within cohort, shared --------------------------------------------------
    rng = np.random.default_rng([336, 1])
    nb = 5 if chk else 200
    boot = {k: [] for k in kinds}
    for kind in kinds:
        sel, w = rep[kind]["sel"], rep[kind]["w"]
        for _ in range(nb):
            idx = {s: rng.integers(0, res[s][kind]["n"], res[s][kind]["n"]) for s in sel}
            v = []
            for which in ("z14", "zc"):
                st = np.array([A.raw_stat(res[s][kind][which][0][idx[s]], res[s][kind][which][1][idx[s]]) for s in sel])
                ok = np.isfinite(st)  # a cohort whose resample is degenerate (zero variance) is left out, weights renormalised
                v.append(float(np.sum(w[ok] * st[ok]) / w[ok].sum()))
            v14, vc = v
            boot[kind].append((v14, vc))
    boot = {k: np.array(v) for k, v in boot.items()}
    # ---- gene drops (pooled statistic) and trio alone -------------------------------------------------------------
    dvals, tvals, gates = {}, {}, {}
    for kind in kinds:
        sel, w = rep[kind]["sel"], rep[kind]["w"]
        dvals[kind] = {lab: float(np.sum(w * [res[s][kind]["drops"][lab] for s in sel]) / w.sum()) for lab, _ in dl}
        assert np.isclose(dvals[kind][tl], core[kind]["obs"])
        tvals[kind] = float(np.sum(w * [res[s][kind]["drops"]["TRIO_ALONE"] for s in sel]) / w.sum())
        try:
            tab = pd.concat([res[s][kind]["gate_rows"] for s in sel], ignore_index=True)
            g = cs.Dataset(tab, anchors=None, waveform=TRIO, name="trio").gate(n_perm=50 if chk else 1000, seed=336)
            gates[kind] = (bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), f"trio gate p_pairing {g.p_pairing:.4f}, p_acrophase {g.p_acrophase:.4f}")
        except Exception as e:  # noqa: BLE001
            gates[kind] = ("not computable", f"trio gate not computable: {type(e).__name__}: {str(e)[:80]}")
    if chk:
        print("check ok")
        return
    # ---- registered criteria ---------------------------------------------------------------------------------------
    dwf = dict(dl)
    holm27 = dict(zip(["myeloid", "T_NK"], A.holm([core[k]["p_pair"] for k in ("myeloid", "T_NK")])))
    holm271 = dict(zip(kinds, A.holm([core[k]["responsive"]["p"] for k in kinds])))
    pr = {k: C.drop_percentiles(dvals[k], tl, dwf, direction=1) for k in kinds}

    def base_row(kind, sid, name, crit, status, extra, stored_val, hyp_note):
        sel = rep[kind]["sel"]
        bt = boot[kind]
        diff = bt[:, 0] - bt[:, 1]
        tv = tvals[kind] if gates[kind][0] is True else np.nan
        return dict(id=sid, figure_panel="Fig 6", statistic=name, direction=1, value_14_stored=stored_val, value_14=rep[kind]["obs"], reproduced=True,
                    value_core=core[kind]["obs"], core_criterion=crit, core_status=status,
                    diff_14_minus_core=rep[kind]["obs"] - core[kind]["obs"], diff_lo90=float(np.quantile(diff, .05)),
                    diff_hi90=float(np.quantile(diff, .95)), diff_draws="200 patient bootstraps within cohort (rng [336, 1]), shared", **pr[kind],
                    trio_gate_passes=gates[kind][0], trio_alone_value=tv, reading=C.reading(status, pr[kind]["trio_drop_below_q10"]), blind="core blind",
                    notes=(f"{len(sel)} cohorts, {int(rep[kind]['w'].sum())} patients; core null mean {core[kind]['nul_mean']:.4f} SD {core[kind]['nul_sd']:.4f}; "
                           f"core pairing p {core[kind]['p_pair']:.4g}; matched core median {core[kind]['matched']['median']:.4f} q95 {core[kind]['matched']['q95']:.4f} "
                           f"p_clock {core[kind]['matched']['p']:.4f}; responsive core median {core[kind]['responsive']['median']:.4f} q95 {core[kind]['responsive']['q95']:.4f} "
                           f"p {core[kind]['responsive']['p']:.4f}; {hyp_note}; {gates[kind][1]}; 14-gene companions refitted for the first {N_CHECK} sets only"),
                    script=SCRIPT)

    rows227 = []
    for kind, sid, name in (("myeloid", "A227_H227.1", "pooled tumour myeloid vs own normal myeloid pairing (24 cohorts)"),
                            ("T_NK", "A227_H227.2", "pooled tumour T/NK vs own normal T/NK pairing (27 cohorts)")):
        t, _ = stored(kind)
        conds = [bool(holm27[kind] < 0.05), bool(core[kind]["obs"] > core[kind]["matched"]["q95"])]
        status = "holds" if all(conds) else "at the line" if any(conds) else "falls"
        rows227.append(base_row(kind, sid, name, "pairing p < 0.05 (Holm over myeloid, T/NK) and statistic above matched core q95", status, None,
                                t["matched"].statistic, f"Holm p {holm27[kind]:.4g} (conditions: pairing {conds[0]}, above q95 {conds[1]})"))
    C.write_rows("a227", rows227)
    rows271, comp_ok = [], []
    for kind, name in (("myeloid", "H271.1 myeloid vs responsive"), ("T_NK", "H271.1 T/NK vs responsive"), ("CD8neg_T", "H271.1 CD8-negative T (H228.4) vs responsive")):
        t, _ = stored(kind)
        conds = [bool(core[kind]["obs"] > core[kind]["responsive"]["q95"]), bool(holm271[kind] < 0.05)]
        status = "holds" if all(conds) else "at the line" if any(conds) else "falls"
        comp_ok.append(conds)
        rows271.append(base_row(kind, f"A271_H271.1:{kind}", name, "statistic above responsive core q95 and responsive p Holm (3) < 0.05", status, None,
                                t["responsive"].statistic, f"responsive Holm p {holm271[kind]:.4g} (above q95 {conds[0]}, Holm {conds[1]})"))
    allc = np.array(comp_ok)
    agg = "holds" if allc.all() else "at the line" if allc.any() else "falls"
    rows271.insert(0, dict(id="A271_H271.1", figure_panel="Fig 6", statistic="Pooled tumour-normal pairing of myeloid, T/NK, CD8-negative T vs 200 responsive sets (aggregate)",
                           direction=1, reproduced=True, core_criterion="every component above the responsive core q95 with responsive Holm p < 0.05",
                           core_status=agg, trio_drop_below_q10=bool(all(pr[k]["trio_drop_below_q10"] for k in kinds)),
                           reading=C.reading(agg, bool(all(pr[k]["trio_drop_below_q10"] for k in kinds))), blind="core blind",
                           notes="aggregate of the three component rows (14 genes: FAIL, Holm 0.597 stored); numeric columns on the component rows", script=SCRIPT))
    C.write_rows("a271_h1", rows271)
    pd.set_option("display.width", 250)
    for r in rows227 + rows271[1:]:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k not in ("notes", "script", "figure_panel", "core_criterion", "diff_draws", "statistic")})
        print(r["notes"])
    print("aggregate", agg, rows271[0]["reading"])


if __name__ == "__main__":
    main()
