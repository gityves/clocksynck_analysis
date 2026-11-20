"""A336 part 1 wrapper: A293 H293.1 / H293.2 (organ order of slow deaths and of ventilated post-mortem donors, Fisher-Lee
rho with fast deaths) on the 14 genes, the core loop, the trio contribution, the 364 gene drops and the trio alone.

Original: A293_organisation_levels/scripts/part1_death_types.py (imported, not edited). Core companions as A335
(`build_core` / `add_resp`). Trio contribution: no bootstrap in the original -> 200 donor bootstraps (fast deaths and the
group resampled, pair mask fixed), seed [336, 1], same resamples for both waveforms. Gate: model-free gate of the trio
on the donors of fast deaths plus the group (both groups' donors pooled). Aggregates only (GTEx controlled access).

    python p1_a293.py [--check]
"""
import argparse
import importlib.util
import time

from p1_util import *  # noqa: F401,F403
import p1_util as U
import a336_common as C
from clocksync import stats as ST
from clocksync.config import limit_threads
import compute_brain_death as B
from common import gtex

spec = importlib.util.spec_from_file_location("part1", U.EXP / "A293_organisation_levels/scripts/part1_death_types.py")
P1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P1)
import compute_coupling_nulls as N  # noqa: E402

REF = pd.read_csv(U.EXP / "A293_organisation_levels/output/part1/groups.csv").set_index("group")
ENTRIES = {"slow_hardy4": "A293_H293.1", "ventilated_postmortem": "A293_H293.2"}


def evaluate(T, groups, fast, masks, n_perm):
    """rho, tissue-label p and q95 of matched and responsive sets, as part1_death_types (same rng order)."""
    rng = np.random.default_rng([293, 1])
    Zc, Mc = P1.raw(T, "clock", 0)
    out = {}
    for g, d in groups.items():
        m = masks[g]
        r0 = P1.rho(Zc, Mc, fast, d, m)
        Hg = N._H_of(Zc, Mc, d)
        th_f = B.phase(Zc, Mc, fast, m)
        null = np.array([ST.fisher_lee(th_f, P1.H.permute_labels(Hg, rng.permutation(len(P1.H.BODY)))[m]) for _ in range(n_perm)])
        row = dict(rho_clock=r0, p_labels=N.pval(r0, null))
        for kind, lab in (("pseudo", "matched"), ("resp", "responsive")):
            rs = np.array([P1.rho(*P1.raw(T, kind, k), fast, d, m) for k in range(T.Z[kind][P1.H.BODY[0]].shape[0])])
            row[f"{lab}_q95"] = float(np.nanquantile(rs, 0.95))
        out[g] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t0 = time.time()
    n_sets, sub, n_perm, n_boot = (20, 300, 50, 3) if chk else (200, None, P1.N_FL, U.N_BOOT_UNIT)
    T, Tc = U.build(n_sets, sub)
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    ns = min(len(od), len(pm))
    rd = np.random.default_rng(132 + P1.SEED)
    od_ = rd.choice(od, ns, replace=False) if len(od) > ns else od
    pm_ = rd.choice(pm, ns, replace=False) if len(pm) > ns else pm
    groups = {"slow_hardy4": T.group_donors("H4"), "ventilated_postmortem": pm_, "hardy3": T.group_donors("H3"),
              "organ_donor": od_}
    Zc0, Mc0 = P1.raw(T, "clock", 0)
    masks = {g: B.pair_mask(Mc0, fast, d) for g, d in groups.items()}
    print(f"states built {time.time() - t0:.0f}s", flush=True)

    # ---- 14 genes, reproduced first -------------------------------------------------------------------------------
    r14 = evaluate(T, groups, fast, masks, n_perm)
    bad = []
    if not chk:
        for g in ENTRIES:
            for c in ("rho_clock", "p_labels", "matched_q95", "responsive_q95"):
                if not np.isclose(r14[g][c], REF.loc[g, c], atol=1e-6, rtol=0):
                    bad.append((g, c, r14[g][c], REF.loc[g, c]))
        print("14-gene reproduction:", "OK" if not bad else bad, flush=True)
    if bad:
        rows = [dict(id=i, figure_panel="Fig 5", statistic="Fisher-Lee rho of tissue-pair phases vs fast deaths",
                     value_14_stored=REF.loc[g, "rho_clock"], value_14=r14[g]["rho_clock"], reproduced=False,
                     notes=f"14-gene values not reproduced (tolerance 1e-6): {bad}", script=str(U.HERE.parent / "p1_a293.py"))
                for g, i in ENTRIES.items()]
        C.write_rows("a293", rows)
        return
    rc = evaluate(Tc, groups, fast, masks, n_perm)

    # ---- trio contribution (200 donor bootstraps, shared resamples) -----------------------------------------------
    Zc14, Mc14 = P1.raw(T, "clock", 0)
    Zcc, Mcc = P1.raw(Tc, "clock", 0)
    rng = np.random.default_rng(U.SEED_BOOT)
    diffs = {g: [] for g in ENTRIES}
    for _ in range(n_boot):
        for g in ENTRIES:
            f, d = rng.choice(fast, len(fast)), rng.choice(groups[g], len(groups[g]))
            diffs[g].append(P1.rho(Zc14, Mc14, f, d, masks[g]) - P1.rho(Zcc, Mcc, f, d, masks[g]))

    # ---- gene-drop null and trio alone ---------------------------------------------------------------------------
    dl, tl = C.drops(n=None if not chk else None)
    if chk:
        keep = dl[:4] + [x for x in dl if x[0] == tl]
        dl = keep
    wfs = dict(dl)
    U.add_states(T, "drop", [w for _, w in dl])
    U.add_states(T, "trio", [C.trio()])
    vals = {g: {} for g in ENTRIES}
    for k, (lab, _) in enumerate(dl):
        Z, M = P1.raw(T, "drop", k)
        for g in ENTRIES:
            vals[g][lab] = P1.rho(Z, M, fast, groups[g], masks[g])
    out = []
    for g, i in ENTRIES.items():
        assert np.isclose(vals[g][tl], rc[g]["rho_clock"], atol=1e-9)
        dp = C.drop_percentiles(vals[g], tl, wfs, +1)
        gp = U.gate(T, np.concatenate([fast, groups[g]]), C.trio(), 1000 if not chk else 50)
        Zt, Mt = P1.raw(T, "trio", 0)
        v_trio = P1.rho(Zt, Mt, fast, groups[g], masks[g])
        c1, c2 = rc[g]["rho_clock"] > rc[g]["matched_q95"], rc[g]["p_labels"] < 0.05
        st = U.status([c1, c2])
        lo, hi = U.ci90(diffs[g])
        out.append(dict(
            id=i, figure_panel="Fig 5 (organ order, slow deaths)" if g == "slow_hardy4" else "Fig 5 (ICU order kept)",
            statistic="Fisher-Lee rho of tissue-pair phases with fast deaths", direction=1,
            value_14_stored=REF.loc[g, "rho_clock"], value_14=r14[g]["rho_clock"], reproduced=True,
            value_core=rc[g]["rho_clock"], core_criterion="rho_core > core matched q95 (decisive) AND tissue-label p < 0.05 (10,000)",
            core_status=st, diff_14_minus_core=r14[g]["rho_clock"] - rc[g]["rho_clock"], diff_lo90=lo, diff_hi90=hi,
            diff_draws="200 donor bootstraps (fast + group), seed [336,1], pair mask fixed", **{k: dp[k] for k in (
                "n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")},
            trio_gate_passes=gp[0], trio_alone_value=v_trio if gp[0] else np.nan,
            reading=C.reading(st, dp["trio_drop_below_q10"]), blind="core blind",
            notes=(f"core matched_q95 {rc[g]['matched_q95']:.4f}, p_labels {rc[g]['p_labels']:.4f}, responsive_q95 {rc[g]['responsive_q95']:.4f} "
                   f"(reported, not decisive); gate p_pairing {gp[1]:.4f}, p_acrophase {gp[2]:.4f} on fast + group donors"),
            script=str(U.HERE.parent / "p1_a293.py")))
    if chk:
        assert len(out) == 2
        print("check ok")
        return
    df = C.write_rows("a293", out)
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
