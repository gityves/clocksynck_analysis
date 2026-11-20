"""A336 part 1, Figure 6: A183 H183.1, H183.2 and A308 H308 (psoriatic lesion delay) on the core loop.

Originals (imported, not edited): `A183_skin_gse154200/scripts/score_a183.py`, `A308_psoriasis_calibration/scripts/score_a308.py`,
the A335 core port `A335_review_batch/scripts/a335_1b_psoriasis.py`. Same data, seeds (A183: SEED 20261101, companions
[SEED, 7|8], 200 sets, 10,000 swaps; A308: SEED 308, 2,000 swaps, set bootstrap [308, 1]); only the waveform changes. Companions:
the original 14-gene draws restricted to the core positions (A329 rule). The offset of H183.1 / H308 uses the fixed
clock-fit weights of the 14-gene cohorts (A308 and A335 form); A335's core numbers are reproduced and checked.

    python p1_a183.py [--check]     # --check: tiny run, assertions only, nothing printed about values nor written
Output: ../output/part1/a183.csv (A183_H183.1, A183_H183.2, A308_H308).
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, PAPER / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M8 = _load("p1_a308", "exploratory/A308_psoriasis_calibration/scripts/score_a308.py")
T335 = _load("p1_a335_1b", "exploratory/A335_review_batch/scripts/a335_1b_psoriasis.py")
A = M8.A
import a336_common as H  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a183"
SCRIPT = "paper/exploratory/A336_two_readouts/scripts/p1_a183.py"
A183 = PAPER / "exploratory/A183_skin_gse154200/output"
T183 = pd.read_csv(A183 / "tests.csv").set_index("hypothesis")
T308 = pd.read_csv(PAPER / "exploratory/A308_psoriasis_calibration/output/tests.csv").set_index("family")
T335C = pd.read_csv(PAPER / "exploratory/A335_review_batch/output/a335_1b_tests.csv")
ST = pd.read_csv(A183 / "h183_1_paired_pooled.csv", index_col=0).value


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


class Ctx:
    def __init__(self):
        self.s = A.load()
        self.means = A.gene_means(self.s)
        self.D = A.Design(self.s)
        self.lab = self.D.labels()
        coh = pd.read_csv(A183 / "h183_1_paired_cohorts.csv", dtype={"cohort": str}).set_index("cohort")
        self.used = str(ST["cohorts"]).split(";")
        self.w = 1 / coh.loc[self.used, "se"].to_numpy() ** 2
        ps = self.s[self.s.diagnosis == "psoriasis"]
        self.data = {}
        for k in self.used:
            case = ps[(ps.cohort == k) & (ps.condition == "lesional")].set_index("unit")
            ref = ps[(ps.cohort == k) & (ps.condition == "nonlesional")].set_index("unit")
            u = sorted(set(case.index) & set(ref.index))
            self.data[k] = (case.loc[u], ref.loc[u])
        self.ps = ps

    def states1(self, wk):
        return [P.cohort_states(*self.data[k], True, wk[k]) for k in self.used]

    def off1(self, wk):
        return P.pooled_offset([P.ccorr(*b) for b in self.states1(wk)], self.w)

    def st2(self, wk):
        return {k: P.cohort_states(self.D.b[k][0], self.D.b[k][1], True, wk[k]) for k in self.D.keys}

    def diff2(self, wk):
        return A.pooled_diff(self.st2(wk), self.lab)[0]


def companions(C, n_sets, restrict_to=None):
    """A183's companion draws (14-gene waveform; stream [SEED, 7|8]); per kind a list of {batch: Waveform}."""
    out = {}
    for kind, sub in (("matched", 7), ("responsive", 8)):
        crng = np.random.default_rng([A.SEED, sub])
        sets = {k: P.pseudo_waveforms(C.means[k], n_sets, crng, H.WF, kind) for k in C.D.keys}
        out[kind] = [{k: sets[k][j] for k in C.D.keys} for j in range(n_sets)]
    return out


def restricted(sets):
    keep = H.core(H.WF)
    return {kind: [{k: H.restrict(w, H.WF, keep) for k, w in wj.items()} for wj in v] for kind, v in sets.items()}


def hyp_tests(C, wfc, sets, n_perm, p1):
    """A183's H183.1 / H183.2 for clock waveform wfc (same waveform for every batch) and companions `sets`."""
    wk = {k: wfc for k in C.D.keys}
    o1, o2 = C.off1(wk), None
    st = C.D.states(wfc)
    o2, offs = A.pooled_diff(st, C.lab)
    null2 = A.perm_null(st, C.lab, n_perm, np.random.default_rng(A.SEED))
    p2 = S.perm_p(null2, o2, "less", tie=0.0)
    rows = []
    comp = {kind: np.array([[C.off1(wj), C.diff2(wj)] for wj in v]) for kind, v in sets.items()}
    for h, obs, p, col in (("H183.1", o1, p1, 0), ("H183.2", o2, p2, 1)):
        row = dict(statistic_h=obs, p=p)
        for kind, v in comp.items():
            row[f"{kind}_q05"] = float(np.quantile(v[:, col], 0.05))
            row[f"{kind}_pctile_of_clock"] = float(np.mean(v[:, col] <= obs))
            row[f"beyond_{kind}_q05"] = bool(obs < row[f"{kind}_q05"])
        rows.append(row)
    order = np.argsort([r["p"] for r in rows])
    run, hol = 0.0, [0, 0]
    for r_, i in enumerate(order):
        run = max(run, min(1.0, (2 - r_) * rows[i]["p"]))
        hol[i] = run
    for r, h in zip(rows, hol):
        r["p_holm"] = h
    return {"H183.1": rows[0], "H183.2": rows[1]}, offs


def q_boot(tab, O, kinds, draws):
    """A308's H308 q (matched) and its set-resample bootstrap draws (stream [308, 1], 2000)."""
    o = O[kinds == "matched"][:, 0]
    sd = tab.loc[tab.kind == "matched", "noise_sd_h"].to_numpy()
    sc, clock = tab.noise_sd_h.iloc[0], O[0][0]

    def q_of(idx):
        tau2 = max(0.0, np.var(o[idx], ddof=1) - np.mean(sd[idx] ** 2))
        return norm.cdf((clock - o[idx].mean()) / np.sqrt(tau2 + sc ** 2))

    brng = np.random.default_rng([M8.SEED, 1])
    qb = np.array([q_of(brng.integers(0, len(o), len(o))) for _ in range(2000)])
    return float(q_of(np.arange(len(o)))), qb


def a308_run(C, n_sets, draws, core_mode):
    """O (sets x draws+1 offsets) for the clock and A183's 200+200 sets (14 genes or restricted to the core)."""
    wf = H.WF
    core = H.core(wf)
    sets = [("clock", 0, {k: (core if core_mode else wf) for k in C.used})]
    for kind, sub in (("matched", 7), ("responsive", 8)):
        crng = np.random.default_rng([A.SEED, sub])
        per = {k: P.pseudo_waveforms(C.means[k], n_sets, crng, wf, kind) for k in C.D.keys}
        sets += [(kind, j, {k: (H.restrict(per[k][j], wf, core) if core_mode else per[k][j]) for k in C.used}) for j in range(n_sets)]
    rng = np.random.default_rng(M8.SEED)
    F = {k: np.vstack([np.zeros(len(C.data[k][0]), bool), rng.random((draws, len(C.data[k][0]))) < 0.5]) for k in C.used}
    wd = dict(zip(C.used, C.w))
    kinds = np.array([k for k, _, _ in sets])
    O = np.empty((len(sets), draws + 1))
    for i, (_, _, wfs) in enumerate(sets):
        O[i] = M8.offsets({k: M8.complex_pairs(*C.data[k], wfs[k]) for k in C.used}, wd, F)
    tab = pd.DataFrame(dict(kind=kinds, set=[j for _, j, _ in sets], offset_h=O[:, 0], noise_sd_h=[M8.wsd(o[1:]) for o in O]))
    return tab, O, kinds


def trio_gate(C):
    """The design's reference gate (200 acrophase permutations, the A183 stream per batch) on the trio, used batches."""
    trio = H.trio()
    p = {}
    for k in C.used:
        ref = C.data[k][1]
        yr = P.log_expr(ref[list(trio.genes)].astype(float), ref.depth.to_numpy(float)).to_numpy()
        try:
            p[k] = P.reference_gate(yr, trio, 200, np.random.default_rng([A.SEED, 1]))["p"]
        except Exception:
            p[k] = np.nan
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(1)
    t0 = time.time()
    n_perm, n_sets, draws, n_boot = (30, 3, 20, 20) if chk else (10_000, 200, 2000, 1000)
    C = Ctx()
    WF, core, trio = H.WF, H.core(H.WF), H.trio(H.WF)
    # ================= stage A: 14 genes =================
    res = P.paired(C.ps, unit="unit", condition="condition", levels=("nonlesional", "lesional"), cohort="cohort", counts=True,
                   pseudo_means=C.means, n_pseudo=n_sets, n_perm=n_perm, seed=A.SEED, pseudo_responsive=True)
    p1_14 = float(res.pooled.p_offset_lower)
    sets14 = companions(C, n_sets)
    T14, offs14 = hyp_tests(C, WF, sets14, n_perm, p1_14)
    tab14, O14, kinds = a308_run(C, n_sets, draws, False)
    q14, qb14 = q_boot(tab14, O14, kinds, draws)
    log("stage A done")
    bad = {}
    for h in ("H183.1", "H183.2"):
        for k, ref_ in (("statistic_h", "statistic_h"), ("p", "p"), ("matched_q05", "matched_q05"), ("responsive_q05", "responsive_q05"), ("p_holm", "p_holm")):
            if not np.isclose(T14[h][k], T183.loc[h, ref_], rtol=0, atol=1e-6):
                bad[f"{h}.{k}"] = (T14[h][k], T183.loc[h, ref_])
    r14 = T308.loc["matched"]
    for k, v in (("clock_offset", O14[0][0]), ("q", q14), ("q05_obs", float(np.quantile(O14[kinds == "matched"], .05, axis=0)[0])),
                 ("q_lo90", float(np.quantile(qb14, .05))), ("q_hi90", float(np.quantile(qb14, .95)))):
        if not np.isclose(v, r14[k], rtol=0, atol=1e-6):
            bad[f"A308.{k}"] = (v, r14[k])
    if chk:
        assert len(sets14["matched"]) == n_sets
    else:
        log("reproduced" if not bad else f"NOT reproduced {bad}")
    if not chk and bad:
        rows = [dict(id=i, figure_panel="Fig 6", statistic=s, direction="-1", value_14_stored=v, value_14=v14, reproduced=not any(k.startswith(pre) for k in bad),
                     notes=("14 genes did not reproduce; no core value computed. " + str({k: x for k, x in bad.items() if k.startswith(pre)})) if any(k.startswith(pre) for k in bad) else "stopped before core (another statistic of this wrapper did not reproduce)",
                     script=SCRIPT)
                for i, s, v, v14, pre in (("A183_H183.1", "psoriasis lesional minus non-lesional offset (h)", T183.loc["H183.1", "statistic_h"], T14["H183.1"]["statistic_h"], "H183.1"),
                                          ("A183_H183.2", "psoriasis minus mean eczema offset (h)", T183.loc["H183.2", "statistic_h"], T14["H183.2"]["statistic_h"], "H183.2"),
                                          ("A308_H308", "q of the clock among matched sets", r14.q, q14, "A308"))]
        H.write_rows(ENTRY, rows)
        return
    # ================= stage B: core, trio, drops =================
    sets_c = restricted(sets14)
    # swap null of H183.1 on the core with a fixed stream (the original's P.paired stream cannot be shared across waveforms)
    cohorts_c = C.states1({k: core for k in C.used})
    null_c = P.swap_offset_null(cohorts_c, C.w, n_perm=n_perm, rng=np.random.default_rng([A.SEED, 3]))
    p1_c = S.perm_p(null_c, P.pooled_offset([P.ccorr(*b) for b in cohorts_c], C.w), "less", tie=0.0)
    cohorts_14 = C.states1({k: WF for k in C.used})
    null_14 = P.swap_offset_null(cohorts_14, C.w, n_perm=n_perm, rng=np.random.default_rng([A.SEED, 3]))
    p1_14_own = S.perm_p(null_14, P.pooled_offset([P.ccorr(*b) for b in cohorts_14], C.w), "less", tie=0.0)
    Tc, offsc = hyp_tests(C, core, sets_c, n_perm, p1_c)
    tabc, Oc, _ = a308_run(C, n_sets, draws, True)
    qc, qbc = q_boot(tabc, Oc, kinds, draws)
    tf, tcm = T335C[T335C.waveform == "14gene"].set_index("family"), T335C[T335C.waveform == "core"].set_index("family")
    chk335 = {"clock_offset": (Oc[0][0], tcm.loc["matched", "clock_offset"]), "q": (qc, tcm.loc["matched", "q"]),
              "q05_obs": (float(np.quantile(Oc[kinds == "matched"], .05, axis=0)[0]), tcm.loc["matched", "q05_obs"]),
              "q_lo90": (float(np.quantile(qbc, .05)), tcm.loc["matched", "q_lo90"]), "q_hi90": (float(np.quantile(qbc, .95)), tcm.loc["matched", "q_hi90"])}
    if not chk:
        for k, (a_, b_) in chk335.items():
            assert np.isclose(a_, b_, rtol=0, atol=1e-6), (k, a_, b_)
    # trio contribution: H183.1 on 200 person bootstraps (stream [336, 1]); H183.2 on A183's own bootstrap (stream [SEED, 9], shared)
    rng = np.random.default_rng([336, 1])
    d1 = []
    b14, bc_ = cohorts_14, cohorts_c
    for _ in range(n_boot if chk else 200):
        pick = [rng.integers(0, len(b[0]), len(b[0])) for b in b14]
        o14 = P.pooled_offset([P.ccorr(*[x[i] for x in b]) for b, i in zip(b14, pick)], C.w)
        oc = P.pooled_offset([P.ccorr(*[x[i] for x in b]) for b, i in zip(bc_, pick)], C.w)
        d1.append(float((o14 - oc + 12) % 24 - 12))
    st14, stc = C.D.states(WF), C.D.states(core)
    r_a, r_b = np.random.default_rng([A.SEED, 9]), np.random.default_rng([A.SEED, 9])
    d2 = []
    for _ in range(n_boot):
        a1 = [A.offset_of(st14, C.lab, d, boot=True, rng=r_a)[0] for d in A.DIAGS]
        a2 = [A.offset_of(stc, C.lab, d, boot=True, rng=r_b)[0] for d in A.DIAGS]
        d2.append(float(A.wrap(A.wrap(a1[0] - 0.5 * (a1[1] + a1[2])) - A.wrap(a2[0] - 0.5 * (a2[1] + a2[2])))))
    # A308: q difference on the shared set resamples
    d3 = qb14 - qbc
    # gene-drop null
    dl, tl = H.drops(WF)
    wfs = dict(dl)
    v1, v2 = {}, {}
    for lab_, w in (dl[:3] if chk else dl):
        try:
            v1[lab_] = C.off1({k: w for k in C.used})
            v2[lab_] = A.pooled_diff(C.D.states(w), C.lab)[0]
        except Exception:
            v1[lab_] = v2[lab_] = np.nan
    if chk:
        assert tl in wfs and H.drops(WF, n=3)[1] == tl
        tg_ = trio_gate(C)
        print("check ok")
        return
    assert np.isclose(v1[tl], Tc["H183.1"]["statistic_h"]) and np.isclose(v2[tl], Tc["H183.2"]["statistic_h"]), (v1[tl], v2[tl])
    pc1, pc2 = H.drop_percentiles(v1, tl, wfs, -1), H.drop_percentiles(v2, tl, wfs, -1)
    tg = trio_gate(C)
    gate_ok = bool(all(np.isfinite(list(tg.values()))) and all(v < 0.05 for v in tg.values())) if np.isfinite(list(tg.values())).any() else "not computable"
    gnote = f"trio gate: reference gate p per used batch {({k: round(v, 3) for k, v in tg.items()})} (3 genes allow 6 acrophase arrangements: p cannot fall below about 1/6)"

    def conds(t, holm):
        return [t["p_holm"] < 0.05, t["beyond_matched_q05"], t["beyond_responsive_q05"]]

    def status(c):
        return "holds" if all(c) else "falls" if not any(c) else "at the line"

    rows = []
    for h, i, desc_, dd, pc in (("H183.1", "A183_H183.1", "psoriasis lesional minus non-lesional clock offset, GSE154200 (h)", d1, pc1),
                                ("H183.2", "A183_H183.2", "psoriasis minus mean eczema offset (h)", d2, pc2)):
        t = Tc[h]
        c = conds(t, None)
        st = status(c)
        dd = np.array(dd)
        note = (f"core: statistic {t['statistic_h']:.4f} h, p {t['p']:.4g}, Holm {t['p_holm']:.4g}, matched q05 {t['matched_q05']:.4f} (clock at percentile {t['matched_pctile_of_clock']:.3f}), "
                f"responsive q05 {t['responsive_q05']:.4f} (percentile {t['responsive_pctile_of_clock']:.3f}); conditions [Holm p<0.05, beyond matched q05, beyond responsive q05] = {c}; "
                + (f"H183.1 core swap null drawn with stream [SEED, 3] (14-gene p with the same stream {p1_14_own:.4g}; stored p from P.paired's stream); " if h == "H183.1" else "")
                + ("A335 core values reproduced to 1e-6 (offset, q05s); " if h == "H183.1" else "") + gnote)
        rows.append(dict(id=i, figure_panel="Fig 6 (psoriasis delay)", statistic=desc_, direction="-1", value_14_stored=T183.loc[h, "statistic_h"],
                         value_14=T14[h]["statistic_h"], reproduced=True, value_core=t["statistic_h"],
                         core_criterion="Holm p < 0.05 (per-person swaps / labels within batch, 10,000) and beyond both companion q05 (clock-specific); holds if all, at the line if some",
                         core_status=st, diff_14_minus_core=T14[h]["statistic_h"] - t["statistic_h"], diff_lo90=float(np.quantile(dd, .05)), diff_hi90=float(np.quantile(dd, .95)),
                         diff_draws=("own: A183 descriptives bootstrap (persons within batch x diagnosis, 1000, stream [SEED, 9], shared)" if h == "H183.2"
                                     else "200 person bootstraps within used batches, seed [336, 1] (A183 has no bootstrap of this statistic)"),
                         **pc, trio_gate_passes=gate_ok, trio_alone_value=np.nan, reading=H.reading(st, pc["trio_drop_below_q10"]),
                         blind="core not blind (A335 / A329)" if h == "H183.1" else "core blind", notes=note, script=SCRIPT))
    rc = tcm.loc["matched"]
    strict = bool(qc < 0.05 and rc.v308_rejection <= 0.0685)
    st3 = "at the line"
    d3 = np.array(d3)
    rows.append(dict(id="A308_H308", figure_panel="Fig 6 (psoriasis delay, calibrated)", statistic="q: clock quantile among noise-free matched sets (lesion-wide rotation null)", direction="-1 (small q = delay beyond the matched family); drop null on the clock offset",
                     value_14_stored=r14.q, value_14=q14, reproduced=True, value_core=qc,
                     core_criterion="q < 0.05 (V308 calibrated: clock below matched q05 in <= 6.85% of complete-null draws); at the line when the 90% interval of q straddles 0.05 (A335's reading)",
                     core_status=st3, diff_14_minus_core=q14 - qc, diff_lo90=float(np.quantile(d3, .05)), diff_hi90=float(np.quantile(d3, .95)),
                     diff_draws="own: 2000 resamples of the matched sets (stream [308, 1]), shared by both waveforms; statistic = q", **pc1, trio_gate_passes=gate_ok,
                     trio_alone_value=np.nan, reading=H.reading(st3, pc1["trio_drop_below_q10"]), blind="core not blind (A335 / A329)",
                     notes=(f"core q {qc:.4f} (90% {np.quantile(qbc, .05):.5f} to {np.quantile(qbc, .95):.4f}), clock offset {Oc[0][0]:.4f} h, matched q05 {chk335['q05_obs'][0]:.4f}, V308 rejection {rc.v308_rejection}; "
                            f"strict reading of q<0.05 and V308: {'holds' if strict else 'falls'}, status 'at the line' because the interval of q reaches above 0.05 (A335); A335 core values reproduced to 1e-6; "
                            "gene-drop null on the clock offset (q needs sets); " + gnote),
                     script=SCRIPT))
    H.write_rows(ENTRY, rows)
    log("done", f"{(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
