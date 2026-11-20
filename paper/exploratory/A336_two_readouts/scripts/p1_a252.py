"""A336 part 1, Figure 1b: A252 H252.1, H252.2, H252.3 (adipose) on the core loop.

Original: `A252_elapsed_time/scripts/compute_a252.py` (imported, not edited). Same data, seeds ([252, 1] label permutations,
[252, tag, sum(ord(acc))] companions, 200 sets, tag 3 matched / 4 responsive), same statistic (pooled rho); only the waveform
changes. Companions: the original 14-gene draws restricted to the core positions of each dataset (A329 rule).

    python p1_a252.py [--check]     # --check: tiny run, assertions only, nothing printed about values nor written
Order: stage A (14 genes only) is compared with `score_tests.csv`; stage B (core, trio, gene drops, bootstrap) only if it
reproduces. Output: ../output/part1/a252.csv.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A252_elapsed_time/scripts"))
sys.path.insert(0, str(HERE.parent))
import compute_a252 as E  # noqa: E402
import a336_common as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a252"
STORED = pd.read_csv(PAPER / "exploratory/A252_elapsed_time/output/score_tests.csv").set_index("family")
ROWS = {"A252_H252.1": ("H252.1 GSE118280 (exact times)", "rho of GSE118280 (exact clock times), pooled over persons x visits"),
        "A252_H252.2": ("H252.2 seven datasets, assumed times", "rho, seven repeat-biopsy datasets pooled"),
        "A252_H252.3_adipose": ("H252.3 adipose", "rho, adipose datasets (GSE142401, GSE303983)")}
CRIT = {"A252_H252.1": "rho > 0, p < 0.05 and rho above matched and responsive q95 (clock-specific)",
        "A252_H252.2": "rho > 0, p < 0.05 and rho above matched and responsive q95 (clock-specific)",
        "A252_H252.3_adipose": "rho > 0, Holm p (over muscle, adipose; as stored) < 0.05 and rho above matched and responsive q95"}
N_SETS, N_PERM, N_BOOT = 200, 10_000, 200


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def fam_rho(parts, ds):
    return E.rho([parts[a] for a in ds])


def stats_parts(D, wf, tested):
    return {a: D[a].sums(E.states(D[a].y, wf)) for a in tested}


def null_pass(D, tested, wfs, n_perm):
    """label permutations (original rng stream) for each waveform in wfs: {label: {a: array(n_perm, 2)}}."""
    Z = {k: {a: E.states(D[a].y, w) for a in tested} for k, w in wfs.items()}
    rng = np.random.default_rng([E.SEED, 1])
    nul = {k: {a: [] for a in tested} for k in wfs}
    for _ in range(n_perm):
        for a in tested:
            el = D[a].perm_el(rng)
            for k in wfs:
                nul[k][a].append(D[a].sums(Z[k][a], el))
    return {k: {a: np.array(v) for a, v in nul[k].items()} for k in wfs}


def fam_null(N, ds):
    num = sum(N[a][:, 0] for a in ds)
    den = sum(N[a][:, 1] for a in ds)
    return num / den


def family_tests(obs_parts, N, comp, fams, n_sets):
    """per family: rho, p, matched / responsive q95; Holm over muscle + adipose as in the original."""
    out = {}
    for k, ds in fams.items():
        obs = fam_rho(obs_parts, ds)
        v = fam_null(N, ds)
        row = dict(rho=obs, p=E.p_up(v, obs))
        for kind in comp:
            s = np.array([E.rho([comp[kind][a][j] for a in ds]) for j in range(n_sets)])
            row[f"{kind}_q95"] = float(np.quantile(s, .95))
        out[k] = row
    ph = E.holm([out["H252.3 muscle"]["p"], out["H252.3 adipose"]["p"]])
    out["H252.3 muscle"]["p_holm"], out["H252.3 adipose"]["p_holm"] = float(ph[0]), float(ph[1])
    return out


def status(i, t):
    cs = [t["rho"] > 0, (t["p_holm"] if "p_holm" in t else t["p"]) < 0.05, t["rho"] > t["matched_q95"] and t["rho"] > t["responsive_q95"]]
    return ("holds" if all(cs) else "falls" if not any(cs) else "at the line"), cs


def unit_contrib(d, z):
    """per-unit (num, den) of rho for one dataset and state z."""
    dl = z - z[d.base]
    num = np.bincount(d.unit, weights=np.where(d.later, (dl * np.conj(d.c)).real, 0.0), minlength=len(d.idx))
    den = np.bincount(d.unit, weights=np.where(d.later, abs(d.c) ** 2, 0.0), minlength=len(d.idx))
    return num, den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t_start = time.time()
    n_sets, n_perm, n_boot = (2, 20, 5) if chk else (N_SETS, N_PERM, N_BOOT)
    tested = [a for a in E.DS if E.DS[a][3] != "descriptive"]
    D = {a: E.Data(a) for a in tested}
    fams = {k: [a for a in ds if a in tested] for k, ds in E.FAMILIES.items()}
    WF = H.WF
    # ---- stage A: 14 genes
    draws = {kind: {a: D[a].pseudo(kind, n_sets, tag, WF) for a in tested} for kind, tag in (("matched", 3), ("responsive", 4))}
    comp14 = {kind: {a: [D[a].sums(E.states(D[a].y, w)) for w in draws[kind][a]] for a in tested} for kind in draws}
    obs14 = stats_parts(D, WF, tested)
    N14 = null_pass(D, tested, {"14": WF}, n_perm)["14"]
    T14 = family_tests(obs14, N14, comp14, fams, n_sets)
    log("stage A done")
    ok, bad = {}, {}
    for i, (fam, _) in ROWS.items():
        s = STORED.loc[fam]
        exp = dict(rho=s.rho, p=s.p, matched_q95=s.matched_q95, responsive_q95=s.responsive_q95)
        if "H252.3" in fam:
            exp["p_holm"] = s.p_holm
        ok[i] = all(np.isclose(T14[fam][k], v, rtol=0, atol=1e-6) for k, v in exp.items())
        bad[i] = {k: (T14[fam][k], v) for k, v in exp.items() if not np.isclose(T14[fam][k], v, rtol=0, atol=1e-6)}
    if chk:
        assert D["GSE118280"].sums(E.states(D["GSE118280"].y, WF), None) == obs14["GSE118280"]
        # vectorised unit contributions agree with the original sums
        for a in tested:
            nu, de = unit_contrib(D[a], E.states(D[a].y, WF))
            assert np.isclose(nu.sum(), obs14[a][0]) and np.isclose(de.sum(), obs14[a][1])
    else:
        for i in ROWS:
            log(i, "reproduced" if ok[i] else f"NOT reproduced {bad[i]}")
    if not chk and not all(ok.values()):
        rows = [dict(id=i, figure_panel="Fig 1b", statistic=ROWS[i][1], direction="+1", value_14_stored=STORED.loc[ROWS[i][0]].rho,
                     value_14=T14[ROWS[i][0]]["rho"], reproduced=bool(ok[i]),
                     notes=("14 genes did not reproduce the stored score_tests.csv row; no core value computed. " + str(bad[i])) if not ok[i]
                     else "another statistic of this wrapper did not reproduce; stopped before core",
                     script="paper/exploratory/A336_two_readouts/scripts/p1_a252.py") for i in ROWS]
        H.write_rows(ENTRY, rows)
        return
    # ---- stage B: core, trio, drops, bootstrap
    core, trio = H.core(WF), H.trio(WF)
    compc = {}
    for kind in draws:
        compc[kind] = {}
        for a in tested:
            base = WF.subset([g for g in WF.genes if g in D[a].y.columns])
            keep = H.core(base)
            compc[kind][a] = [D[a].sums(E.states(D[a].y, H.restrict(w, base, keep))) for w in draws[kind][a]]
    obsc = stats_parts(D, core, tested)
    NB = null_pass(D, tested, {"core": core, "trio": trio}, n_perm)
    Tc = family_tests(obsc, NB["core"], compc, fams, n_sets)
    # trio alone: own clock-plane test (permuted labels), no companions needed
    obst = stats_parts(D, trio, tested)
    trio_p = {}
    for fam, ds in fams.items():
        obs = fam_rho(obst, ds)
        trio_p[fam] = (obs, E.p_up(fam_null(NB["trio"], ds), obs))
    # bootstrap of units (shared by both waveforms)
    U14 = {a: unit_contrib(D[a], E.states(D[a].y, WF)) for a in tested}
    Uc = {a: unit_contrib(D[a], E.states(D[a].y, core)) for a in tested}
    rng = np.random.default_rng([336, 1])
    bd = {fam: [] for fam in fams}
    for _ in range(n_boot):
        pick = {a: rng.integers(0, len(D[a].idx), len(D[a].idx)) for a in tested}
        for fam, ds in fams.items():
            f = lambda U: sum(U[a][0][pick[a]].sum() for a in ds) / sum(U[a][1][pick[a]].sum() for a in ds)  # noqa: E731
            bd[fam].append(f(U14) - f(Uc))
    # gene-drop null
    dl, tl = H.drops(WF)
    wfs = dict(dl)
    vals = {fam: {} for fam in fams}
    for lab, w in (dl[:3] if chk else dl):
        try:
            parts = stats_parts(D, w, tested)
            for fam, ds in fams.items():
                vals[fam][lab] = fam_rho(parts, ds)
        except Exception:
            for fam in fams:
                vals[fam][lab] = np.nan
    if chk:
        assert tl in wfs
        assert H.drops(WF, n=3)[1] == tl
        for lab, w in H.drops(WF, n=3)[0]:
            if lab == tl:
                assert np.isclose(fam_rho(stats_parts(D, w, tested), fams["H252.1 GSE118280 (exact times)"]), Tc["H252.1 GSE118280 (exact times)"]["rho"])
        print("check ok")
        return
    rows = []
    for i, (fam, desc) in ROWS.items():
        assert np.isclose(vals[fam][tl], Tc[fam]["rho"]), (fam, vals[fam][tl], Tc[fam]["rho"])
        pc = H.drop_percentiles(vals[fam], tl, wfs, +1)
        st, conds = status(i, Tc[fam])
        d = np.array(bd[fam])
        tg = trio_p[fam][1] < 0.05
        t = Tc[fam]
        pname = "p_holm" if "p_holm" in t else "p"
        note = (f"core: rho {t['rho']:.4f}, p {t['p']:.4g}" + (f", p_holm {t['p_holm']:.4g}" if "p_holm" in t else "") +
                f", matched q95 {t['matched_q95']:.4f}, responsive q95 {t['responsive_q95']:.4f}; conditions [rho>0, {pname}<0.05, clock-specific] = {conds}; "
                f"trio alone: rho {trio_p[fam][0]:.4f}, permutation p {trio_p[fam][1]:.4g}; core drops use the 364-drop null; "
                "14-gene and core use the same label permutations (seed [252, 1]), companions are the original draws restricted to core positions")
        rows.append(dict(id=i, figure_panel="Fig 1b", statistic=desc, direction="+1", value_14_stored=STORED.loc[fam].rho,
                         value_14=T14[fam]["rho"], reproduced=True, value_core=t["rho"], core_criterion=CRIT[i], core_status=st,
                         diff_14_minus_core=T14[fam]["rho"] - t["rho"], diff_lo90=float(np.quantile(d, .05)), diff_hi90=float(np.quantile(d, .95)),
                         diff_draws="own: 200 unit bootstraps (units = person x visit/leg, within dataset), seed [336, 1]", **pc,
                         trio_gate_passes=bool(tg), trio_alone_value=trio_p[fam][0] if tg else np.nan,
                         reading=H.reading(st, pc["trio_drop_below_q10"]), blind="core blind", notes=note,
                         script="paper/exploratory/A336_two_readouts/scripts/p1_a252.py"))
    H.write_rows(ENTRY, rows)
    log("done", f"{(time.time() - t_start) / 60:.1f} min")


if __name__ == "__main__":
    main()
