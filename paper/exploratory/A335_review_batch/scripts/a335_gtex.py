"""A335 (registered), GTEx parts 1d, 1e, 2, 3 (H315.1) and 5: core-loop reruns, covariate balance, equivalence, D3.

Core loop = `cs.gtex().drop(["PER1", "PER2", "NR1D1"])`; matched and responsive sets are drawn exactly as in A315 (same
seeds, 14-gene waveform) and restricted to the genes at the core positions (A329 companion rule). The 14-gene values are
recomputed first and must reproduce A315 / A332 D2 (`tests.csv`, `descriptive.csv`, `d2_per_organ.csv`); otherwise the
script stops before any core statistic is computed. Aggregates only (GTEx controlled access).

    python a335_gtex.py [--check]     # --check: subsample, no reproduction stop, nothing printed or saved
Outputs (`../output/`): repro_14gene.csv, a315_core.csv, a315_core_tissues.csv, d2_core.csv, balance.csv, matched.csv,
equivalence.csv, d3_tertiles.csv.
"""
import argparse
import copy
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
for d in ("A109_coupling_nulls", "A132_brain_death", "A315_braindeath_gain", "A317_braindeath_confirm", "A231_age_coupling"):
    sys.path.insert(0, str(EXP / d / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_brain_death as B  # noqa: E402
import compute_h109_responsive as R  # noqa: E402
import d2_per_organ as D2  # noqa: E402
import matched_label as ML  # noqa: E402
import score_a315 as A315  # noqa: E402
import score_a317 as S317  # noqa: E402
import score_a325 as S325  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
WF = H.WF
DROP = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(DROP)
CORE_IDX = [k for k, g in enumerate(WF.genes) if g in CORE.genes]
assert len(CORE_IDX) == 11
MARGIN = 0.357
HOURS = H.HOURS


def core_wf(w, j=0):
    return cs.Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, f"core{j}")


def build_core(T, Tc_resp_sets):
    """Core states of the clock, the matched sets (and, passed in, the responsive sets): dict kind -> tissue -> (Z, V)."""
    Z, V = {"clock": {}, "pseudo": {}}, {"clock": {}, "pseudo": {}}
    for t in H.BODY:
        dep = T.rows[t].depth.to_numpy(float)
        Z["clock"][t], V["clock"][t] = H._state_set(T.rows[t], [CORE], dep)
        z, v = [], []
        for j, g in enumerate(T.sets[t]):
            w = core_wf(WF, j)
            gg = [g[k] for k in CORE_IDX]
            xx = T.pseudo[t][gg].astype(float)
            L = H.P.projector(xx, cs.Waveform(tuple(gg), w.amp, w.psi, None, None, w.name), counts=True)
            b, vv = H.P.states(xx, L, cs.Waveform(tuple(gg), w.amp, w.psi, None, None, w.name), depth=dep)
            z.append(b[:, 0] + 1j * b[:, 1])
            v.append(vv)
        Z["pseudo"][t], V["pseudo"][t] = np.array(z), np.array(v)
    return Z, V


def add_resp(T, Zc, Vc, n_sets=200, seed=0):
    """Responsive sets drawn as `compute_h109_responsive.add_sets` (same rng), states for 14 genes (T) and core (Zc, Vc)."""
    T.Z["resp"], T.V["resp"], Zc["resp"], Vc["resp"] = {}, {}, {}, {}
    for k, t, X, lm_all, names in R.reread_rows(T):
        dep = T.rows[t].depth.to_numpy(float)
        lm_r = lm_all[~lm_all.index.duplicated()]
        col = {g: i for i, g in enumerate(names)}
        ws = H.P.pseudo_waveforms(lm_r, n_sets, np.random.default_rng([109, 7, seed, k]), WF, kind="responsive")
        wc = [cs.Waveform(tuple(w.genes[i] for i in CORE_IDX), CORE.amp, CORE.psi, None, None, f"core{j}")
              for j, w in enumerate(ws)]
        for ww, Zd, Vd in ((ws, T.Z, T.V), (wc, Zc, Vc)):
            u = sorted({g for w in ww for g in w.genes})
            x = pd.DataFrame(X[:, [col[g] for g in u]].astype(float), columns=u)
            Zd["resp"][t], Vd["resp"][t] = H._state_set(x, ww, dep)
        del X


def arrows(T, nboot):
    """D2 arrows (as `d2_per_organ.main`) on the states of T: estimate, 90 % bootstrap interval per tissue (h)."""
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(S317.SEED)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(S317.N_DRAWS)]
    Z, M = S317.raw(T, "clock", 0, T.R)
    mask = (B.pair_mask(M, fast, pm, *od_draws[:1]) & B.pair_mask(M, od)).astype(bool)
    obs = np.array([D2.per_tissue(Z, M, d, pm, mask) for d in od_draws])
    est = np.angle(np.nanmean(obs / np.abs(obs), 0)) * HOURS
    rb = np.random.default_rng([332, 2])
    boot = []
    for _ in range(nboot):
        p_, o_ = rb.choice(pm, len(pm)), rb.choice(od, len(pm))
        boot.append(np.angle(D2.per_tissue(Z, M, o_, p_, mask)) * HOURS)
    boot = np.array(boot)
    dev = np.angle(np.exp(1j * (boot - est) / HOURS)) * HOURS
    return pd.DataFrame(dict(tissue=H.BODY, phase_change_h=est, lo90=est + np.nanquantile(dev, 0.05, 0),
                             hi90=est + np.nanquantile(dev, 0.95, 0)))


def contrast(T, a, b, rng=None, nboot=0):
    """Q, matched and responsive Q sets, Q*, percentiles and q95 of the matched Q for donors a (vs b)."""
    q = float(A315.Q(A315.L_table(T, "clock", a, b))[0])
    Qm, Qr = A315.Q(A315.L_table(T, "pseudo", a, b)), A315.Q(A315.L_table(T, "resp", a, b))
    out = dict(Q=q, matched_median=float(np.nanmedian(Qm)), matched_q95=float(np.nanquantile(Qm, 0.95)),
               Q_star=float(q - np.nanmedian(Qm)), pct_matched=float(np.nanmean(Qm < q)),
               resp_median=float(np.nanmedian(Qr)), resp_q95=float(np.nanquantile(Qr, 0.95)),
               pct_resp=float(np.nanmean(Qr < q)), tissues=int(A315.L_table(T, "clock", a, b)[0].notna().sum()),
               n_a=len(a), n_b=len(b))
    if rng is not None:  # H315.1 bootstrap exactly as score_a315
        boot = np.array([A315.Q(A315.L_table(T, "clock", rng.choice(a, len(a)), rng.choice(b, len(b))))[0]
                         for _ in range(nboot)])
        out.update(boot_lo90=float(np.nanquantile(boot, 0.05)), boot_hi90=float(np.nanquantile(boot, 0.95)),
                   matched_q05=float(np.nanquantile(Qm, 0.05)))
    return out


def qstar_boot(T, a, b, rng, n):
    """Q* = Q_clock - median Q_matched, donors of both groups resampled; clock and matched sets share the draw."""
    v = []
    for _ in range(n):
        ra, rb = rng.choice(a, len(a)), rng.choice(b, len(b))
        v.append(A315.Q(A315.L_table(T, "clock", ra, rb))[0] - np.nanmedian(A315.Q(A315.L_table(T, "pseudo", ra, rb))))
    return np.array(v)


def smd(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / s) if s > 0 else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    global OUT
    if chk:
        OUT = Path('/tmp/claude-1000/-home-maxine-Documents-paychere-scCircadianOverview/836c0a38-01a1-4fb5-aefe-3ce8718be187/scratchpad/a335_check')
    limit_threads(4)
    t0 = time.time()
    log = (lambda *a, **k: None) if chk else (lambda *a, **k: print(*a, **k, flush=True))
    nboot = 3 if chk else 1000
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    Tc = copy.copy(T)
    Zc, Vc = build_core(T, None)
    add_resp(T, Zc, Vc)
    Tc.Z, Tc.V = Zc, Vc
    print(f"states built {time.time() - t0:.0f}s", flush=True)
    ph = gtex.phenotypes(B.PH_COLS + ["AGE", "SEX"]).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    print(f"organ donors {len(od)}, post-mortem {len(pm)}, fast {len(fast)}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- reproduction of the 14-gene values (stop before any core statistic) ---------------------------------------
    rep = []
    r14 = contrast(T, od, pm, np.random.default_rng(315), nboot)
    i14 = contrast(T, fast, pm)
    a14 = arrows(T, nboot)
    if not chk:
        t_ref = pd.read_csv(EXP / "A315_braindeath_gain/output/tests.csv").iloc[0]
        d_ref = pd.read_csv(EXP / "A315_braindeath_gain/output/descriptive.csv").set_index("contrast")
        a_ref = pd.read_csv(EXP / "A317_braindeath_confirm/output/d2_per_organ.csv").set_index("tissue")
        rep = [("H315.1 Q", r14["Q"], t_ref.Q), ("H315.1 Q*", r14["Q_star"], t_ref.Q_star),
               ("H315.1 pct matched", r14["pct_matched"], t_ref.percentile_vs_matched),
               ("H315.1 boot_lo90", r14["boot_lo90"], t_ref.boot_lo90), ("H315.1 boot_hi90", r14["boot_hi90"], t_ref.boot_hi90),
               ("ICU Q", i14["Q"], d_ref.loc["fast_vs_postmortem", "Q"]),
               ("ICU Q*", i14["Q_star"], d_ref.loc["fast_vs_postmortem", "Q_star"]),
               ("ICU pct matched", i14["pct_matched"], d_ref.loc["fast_vs_postmortem", "percentile_vs_matched"])]
        for _, r in a14.iterrows():
            for c, rc in (("phase_change_h", "phase_change_h"), ("lo90", "phase_lo90_h"), ("hi90", "phase_hi90_h")):
                rep.append((f"D2 {r.tissue} {c}", r[c], a_ref.loc[r.tissue, rc]))
        rep = pd.DataFrame(rep, columns=["quantity", "recomputed", "stored"])
        rep["match"] = np.isclose(rep.recomputed, rep.stored, rtol=1e-6, atol=1e-6)
        rep.to_csv(OUT / "repro_14gene.csv", index=False)
        print(rep.round(4).to_string(index=False), flush=True)
        if not rep.match.all():
            raise SystemExit("A315 / A332 14-gene values not reproduced: stop before reading the core loop")
    print(f"14-gene reproduction done {time.time() - t0:.0f}s", flush=True)

    # ---- 1d core --------------------------------------------------------------------------------------------------
    rc = contrast(Tc, od, pm, np.random.default_rng(315), nboot)
    ic = contrast(Tc, fast, pm)
    rows = []
    for gs, r, i in (("14 genes", r14, i14), ("core", rc, ic)):
        rows.append(dict(gene_set=gs, test="H315.1 organ donors vs ventilated PM", **r,
                         passes=bool(r["Q"] < 0 and r["Q"] < r["matched_q05"] and r["boot_hi90"] < 0)))
        rows.append(dict(gene_set=gs, test="ICU fast vs ventilated PM", **i,
                         icu_loss_holds=bool(i["Q_star"] > 0 and i["Q"] > i["matched_q95"])))
    pd.DataFrame(rows).to_csv(OUT / "a315_core.csv", index=False)
    log(pd.DataFrame(rows).round(4).T.to_string())
    tab = pd.concat({"core_od_vs_pm": A315.L_table(Tc, "clock", od, pm)[0], "core_fast_vs_pm": A315.L_table(Tc, "clock", fast, pm)[0]},
                    axis=1)
    tab.to_csv(OUT / "a315_core_tissues.csv")
    print(f"1d done {time.time() - t0:.0f}s", flush=True)

    # ---- 1e core arrows -------------------------------------------------------------------------------------------
    ac = arrows(Tc, nboot)
    o = a14.merge(ac, on="tissue", suffixes=("_14", "_core"))
    o["holds"] = (np.sign(o.phase_change_h_core) == np.sign(o.phase_change_h_14)) & ((o.lo90_core > 0) | (o.hi90_core < 0))
    o.to_csv(OUT / "d2_core.csv", index=False)
    log(o.round(3).to_string())
    print(f"1e done {time.time() - t0:.0f}s", flush=True)

    # ---- 2 covariate balance and matched contrasts ----------------------------------------------------------------
    attrs = S317.sample_attrs(T)
    med = lambda col: np.array([np.nanmedian([attrs[t][col].to_numpy(float)[T.R[t][d]] for t in H.BODY if T.R[t][d] >= 0]
                                             or [np.nan]) for d in range(len(T.donors))])
    q = pd.DataFrame(dict(age=pd.to_numeric(ph.AGE, errors="coerce").to_numpy(), male=(pd.to_numeric(ph.SEX, errors="coerce") == 1).astype(int).to_numpy(),
                          rin=med("rin"), ischaemic_min=med("ischemic_time"), hardy=0), index=np.arange(len(T.donors)))
    q["ischaemic_h"] = q.ischaemic_min / 60
    brow, mrow, cases_used = [], [], {}
    for name, ctrl in (("organ donors vs ventilated PM", od), ("fast vs ventilated PM", fast)):
        cas = pd.Series(np.isin(q.index, pm), index=q.index)
        con = pd.Series(np.isin(q.index, ctrl), index=q.index)
        sets = ML.match_ratio(q, cas, con)
        cd = [s[0] for s in sets]
        ctl = [x for s in sets for x in s[1:]]
        for v in ("age", "male", "rin", "ischaemic_min"):
            brow.append(dict(contrast=name, variable=v, smd_before=smd(q.loc[pm, v], q.loc[ctrl, v]),
                             smd_after=smd(q.loc[cd, v], q.loc[ctl, v])))
        cd, ctl = np.array(cd), np.array(sorted(set(ctl)))
        cases_used[name] = (ctl, cd)
        for gs, TT in (("14 genes", T), ("core", Tc)):
            r = contrast(TT, ctl, cd)
            mrow.append(dict(contrast=name, gene_set=gs, n_cases=len(cd), n_controls_distinct=len(ctl),
                             n_cases_pool=len(pm), n_controls_pool=len(ctrl), **r,
                             icu_robust=bool(name.startswith("fast") and r["Q_star"] > 0 and r["Q"] > r["matched_q95"])))
    pd.DataFrame(brow).to_csv(OUT / "balance.csv", index=False)
    pd.DataFrame(mrow).to_csv(OUT / "matched.csv", index=False)
    log(pd.DataFrame(brow).round(3).to_string(index=False))
    log(pd.DataFrame(mrow).round(4).T.to_string())
    print(f"2 done {time.time() - t0:.0f}s", flush=True)

    # ---- 3 H315.1 equivalence --------------------------------------------------------------------------------------
    eq = []
    odm, pmm = cases_used["organ donors vs ventilated PM"]
    for gs, TT in (("14 genes", T), ("core", Tc)):
        for dn, (a, b) in (("pooled", (od, pm)), ("matched donors", (odm, pmm))):
            rng = np.random.default_rng([335, 3, len(a), len(b)])
            v = qstar_boot(TT, a, b, rng, nboot)
            lo, hi = float(np.nanquantile(v, 0.05)), float(np.nanquantile(v, 0.95))
            eq.append(dict(gene_set=gs, donors=dn, n_a=len(a), n_b=len(b),
                           Q_star=contrast(TT, a, b)["Q_star"], lo90=lo, hi90=hi, margin=MARGIN,
                           equal=bool(lo > -MARGIN and hi < MARGIN)))
    pd.DataFrame(eq).to_csv(OUT / "equivalence.csv", index=False)
    log(pd.DataFrame(eq).round(4).to_string(index=False))
    print(f"3 done {time.time() - t0:.0f}s", flush=True)

    # ---- 5 D3 remainder --------------------------------------------------------------------------------------------
    dur = pd.Series(pd.to_numeric(ph.DTHVNTD, errors="coerce").to_numpy(), index=np.arange(len(T.donors)))
    day = (T.tod >= 8) & (T.tod < 20)
    strata = {"08-20 h": pm[day[pm]], "20-08 h": pm[~day[pm]]}
    top, bot = S325.tertiles(dur, strata)
    ZZ, MM = S317.raw(T, "clock", 0, T.R)
    mask = B.pair_mask(MM, fast, top, bot)
    d3 = dict(n_pm=len(pm), n_with_duration=int(dur.reindex(pm).notna().sum()), n_dropped=int(dur.reindex(pm).isna().sum()),
              n_top=len(top), n_bottom=len(bot), pairs=int(mask.sum()),
              n_top_day=int(np.isin(top, strata["08-20 h"]).sum()), n_top_night=int(np.isin(top, strata["20-08 h"]).sum()),
              n_bottom_day=int(np.isin(bot, strata["08-20 h"]).sum()), n_bottom_night=int(np.isin(bot, strata["20-08 h"]).sum()),
              dur_median_top=float(dur.reindex(top).median()), dur_median_bottom=float(dur.reindex(bot).median()))
    rows = []
    for gs, TT in (("14 genes", T), ("core", Tc)):
        Z_, M_ = S317.raw(TT, "clock", 0, TT.R)
        c = contrast(TT, top, bot)
        rows.append(dict(gene_set=gs, rho_top=S317.rho(Z_, M_, fast, top, mask), rho_bottom=S317.rho(Z_, M_, fast, bot, mask),
                         Q=c["Q"], matched_median=c["matched_median"], Q_star=c["Q_star"], pct_matched=c["pct_matched"],
                         tissues=c["tissues"], **d3))
    pd.DataFrame(rows).to_csv(OUT / "d3_tertiles.csv", index=False)
    log(pd.DataFrame(rows).round(3).T.to_string())
    print(f"all done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
