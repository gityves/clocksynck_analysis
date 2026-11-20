"""A336 part 1 wrapper: A315 ICU loudness loss (fast deaths vs ventilated post-mortem) and A315 H315.1 (brain-death
loudness, Q*, equivalence margin 0.357 as A335) on the 14 genes, the core loop, the trio contribution, the 364 gene drops
and the trio alone.

Reuses a335_gtex (`contrast`, `qstar_boot`, core companions). Trio contribution: 1,000 donor bootstraps of Q* (clock minus
median matched), same resamples for both waveforms; seed as A335's equivalence for H315.1 ([335, 3, n_a, n_b], which
reproduces equivalence.csv) and the same form for the ICU contrast. Gene-drop null: clock Q only (no matched median),
direction +1. Trio alone: Q* with matched sets at the trio positions, reported if the gate on the contrast's donors passes.
Aggregates only (GTEx controlled access).

    python p1_a315.py [--check]
"""
import argparse
import time

from p1_util import *  # noqa: F401,F403
import p1_util as U
import a336_common as C
import score_a315 as A315
import compute_brain_death as B
from clocksync.config import limit_threads
from common import gtex

OUTD = U.EXP / "A335_review_batch" / "output"
MARGIN = U.G.MARGIN


def qstar(T, a, b):
    return float(A315.Q(A315.L_table(T, "clock", a, b))[0] - np.nanmedian(A315.Q(A315.L_table(T, "pseudo", a, b))))


def paired_boot(T, Tc, a, b, n):
    rng = np.random.default_rng([335, 3, len(a), len(b)])
    out = []
    for _ in range(n):
        ra, rb = rng.choice(a, len(a)), rng.choice(b, len(b))
        out.append((qstar(T, ra, rb), qstar(Tc, ra, rb)))
    return np.array(out)


def trio_companion_T(T):
    """Trio states and matched-set states at the trio positions (G.build_core with CORE swapped), as a T copy."""
    import copy
    keep = (U.G.CORE, U.G.CORE_IDX)
    tw = C.trio()
    U.G.CORE, U.G.CORE_IDX = tw, [k for k, g in enumerate(U.WF.genes) if g in tw.genes]
    try:
        Z, V = U.G.build_core(T, None)
    finally:
        U.G.CORE, U.G.CORE_IDX = keep
    Tt = copy.copy(T)
    Tt.Z, Tt.V = Z, V
    return Tt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t0 = time.time()
    n_sets, sub, n_boot = (20, 300, 3) if chk else (200, None, 1000)
    T, Tc = U.build(n_sets, sub)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    print(f"states built {time.time() - t0:.0f}s: od {len(od)}, pm {len(pm)}, fast {len(fast)}", flush=True)
    stat = {"A315_H315.1": ("organ donors vs ventilated PM", od), "A315_ICU": ("fast vs ventilated PM", fast)}

    # ---- 14 genes, reproduced first (A315 tests.csv, descriptive.csv, A335 equivalence.csv) ------------------------
    t_ref = pd.read_csv(U.EXP / "A315_braindeath_gain/output/tests.csv").iloc[0]
    d_ref = pd.read_csv(U.EXP / "A315_braindeath_gain/output/descriptive.csv").set_index("contrast").loc["fast_vs_postmortem"]
    eq = pd.read_csv(OUTD / "equivalence.csv")
    eq = eq[eq.donors == "pooled"].set_index("gene_set")
    boots = {k: paired_boot(T, Tc, a, pm, n_boot) for k, (_, a) in stat.items()}
    c14 = {k: G_ for k, G_ in ((k, U.G.contrast(T, a, pm)) for k, (_, a) in stat.items())}
    chk14 = {"H315.1 Q": (c14["A315_H315.1"]["Q"], t_ref.Q), "H315.1 Q*": (c14["A315_H315.1"]["Q_star"], t_ref.Q_star),
             "H315.1 pct": (c14["A315_H315.1"]["pct_matched"], t_ref.percentile_vs_matched),
             "H315.1 eq lo90": (U.ci90(boots["A315_H315.1"][:, 0])[0], eq.loc["14 genes", "lo90"]),
             "H315.1 eq hi90": (U.ci90(boots["A315_H315.1"][:, 0])[1], eq.loc["14 genes", "hi90"]),
             "ICU Q": (c14["A315_ICU"]["Q"], d_ref.Q), "ICU Q*": (c14["A315_ICU"]["Q_star"], d_ref.Q_star),
             "ICU pct": (c14["A315_ICU"]["pct_matched"], d_ref.percentile_vs_matched)}
    bad = {k: v for k, v in chk14.items() if not np.isclose(*v, atol=1e-6, rtol=0)}
    if not chk:
        print("14-gene reproduction:", "OK" if not bad else bad, flush=True)
    if bad and not chk:
        rows = [dict(id=i, statistic="Q*", value_14_stored=t_ref.Q_star if "H315" in i else d_ref.Q_star,
                     value_14=c14[i]["Q_star"], reproduced=False, notes=f"14-gene values not reproduced: {bad}",
                     script=str(U.HERE.parent / "p1_a315.py")) for i in stat]
        C.write_rows("a315", rows)
        return

    # ---- core, reproduced against A335 a315_core.csv / equivalence.csv --------------------------------------------
    cc = {k: U.G.contrast(Tc, a, pm) for k, (_, a) in stat.items()}
    if not chk:
        ref = pd.read_csv(OUTD / "a315_core.csv")
        r1 = ref[(ref.gene_set == "core") & ref.test.str.startswith("H315.1")].iloc[0]
        r2 = ref[(ref.gene_set == "core") & ref.test.str.startswith("ICU")].iloc[0]
        for k, r in (("A315_H315.1", r1), ("A315_ICU", r2)):
            assert np.isclose(cc[k]["Q"], r.Q, atol=1e-9) and np.isclose(cc[k]["Q_star"], r.Q_star, atol=1e-9), "A335 core not reproduced"
        assert np.isclose(U.ci90(boots["A315_H315.1"][:, 1])[0], eq.loc["core", "lo90"], atol=1e-9)

    # ---- gene drops (clock Q only) and trio alone ----------------------------------------------------------------
    dl, tl = C.drops()
    if chk:
        dl = dl[:4] + [x for x in dl if x[0] == tl]
    wfs = dict(dl)
    U.add_states(T, "drop", [w for _, w in dl])
    Tt = trio_companion_T(T)
    out = []
    for k, (name, a) in stat.items():
        vals = {lab: float(A315.Q(A315.L_table(T, "drop", a, pm))[j]) for j, (lab, _) in enumerate(dl)}
        assert np.isclose(vals[tl], cc[k]["Q"], atol=1e-9)
        dp = C.drop_percentiles(vals, tl, wfs, +1)
        gp = U.gate(T, np.concatenate([a, pm]), C.trio(), 1000 if not chk else 50)
        v_trio = qstar(Tt, a, pm)
        lo, hi = U.ci90(boots[k][:, 0] - boots[k][:, 1])
        base = dict(id=k, figure_panel="Fig 5 (brain-death loudness)" if "H315" in k else "Fig 5 (ICU loudness loss)",
                    direction=1, value_14=c14[k]["Q_star"], reproduced=True, value_core=cc[k]["Q_star"],
                    diff_14_minus_core=c14[k]["Q_star"] - cc[k]["Q_star"], diff_lo90=lo, diff_hi90=hi,
                    diff_draws="1000 donor bootstraps of Q* (clock minus median matched), shared resamples, seed [335,3,n_a,n_b]",
                    trio_gate_passes=gp[0], trio_alone_value=v_trio if gp[0] else np.nan, blind="core not blind (A335 / A329)",
                    script=str(U.HERE.parent / "p1_a315.py"), **{x: dp[x] for x in (
                        "n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")})
        if k == "A315_H315.1":
            lo_c, hi_c = U.ci90(boots[k][:, 1])
            st = "holds" if (lo_c > -MARGIN and hi_c < MARGIN) else "falls"
            base.update(statistic="Q* (clock minus median matched Q), organ donors vs ventilated PM", value_14_stored=t_ref.Q_star,
                        core_criterion=f"core 90% bootstrap interval of Q* inside +-{MARGIN}", core_status=st,
                        notes=f"core Q* 90% {lo_c:.4f} to {hi_c:.4f}; 14 genes {U.ci90(boots[k][:, 0])[0]:.4f} to {U.ci90(boots[k][:, 0])[1]:.4f}; "
                              f"drop statistic = clock Q, direction +1 (original Q* positive); gate {gp[1]:.4f}/{gp[2]:.4f} on od + pm donors")
        else:
            c1, c2 = cc[k]["Q_star"] > 0, cc[k]["Q"] > cc[k]["matched_q95"]
            st = U.status([c1, c2])
            base.update(statistic="Q* fast deaths vs ventilated PM (ICU loudness loss)", value_14_stored=d_ref.Q_star,
                        core_criterion="Q* > 0 AND Q above matched q95", core_status=st,
                        notes=f"core matched q95 {cc[k]['matched_q95']:.4f}, pct matched {cc[k]['pct_matched']:.3f}; "
                              f"bootstrap = same form as H315.1 (fast vs pm); gate {gp[1]:.4f}/{gp[2]:.4f} on fast + pm donors")
        base["core_status"] = st
        base["reading"] = C.reading(st, dp["trio_drop_below_q10"])
        out.append(base)
    if chk:
        print("check ok")
        return
    df = C.write_rows("a315", out)
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
