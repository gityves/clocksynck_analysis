"""A336 part 1 wrapper: A331 H331.1 (person loudness, MuTHER fat x skin leverage) on the 14 genes, the core loop (A335 1c),
the trio contribution, the 364 gene drops and the trio alone.

Original: P1_person_loudness/scripts/score_a331.py (imported) and A335_review_batch/scripts/a335_1c_muther.py. The
permutation null and the twin-pair bootstrap replay `score_a331` (one rng(331): 2,000 persons permutations of the skin
states, then 2,000 twin-pair-unit resamples); both waveforms share the permutations and the resamples. Trio contribution:
the 2,000 twin-pair bootstrap draws. Companions: matched / responsive sets (tags 3 / 4) of the cohort waveform,
restricted to the core positions. Gene drop: clock leverage only, direction +1. Trio alone: A255's reference gate
(`paired.reference_gate`, clock-plane share of the covariance within the cohort's samples, 1,000 permutations) of the
trio in fat and skin, both p < 0.05 (MuTHER is an array cohort, not a donor x compartment pseudobulk).

    python p1_a331.py [--check]
"""
import argparse
import importlib.util
import time
import warnings

from p1_util import *  # noqa: F401,F403
import p1_util as U
import a336_common as C
from clocksync import elasticity as E
from clocksync import paired as P
from clocksync.config import limit_threads

spec = importlib.util.spec_from_file_location("a331", U.EXP / "P1_person_loudness/scripts/score_a331.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
OUTD = U.EXP / "P1_person_loudness" / "output"
SCRIPT = str(U.HERE.parent / "p1_a331.py")


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t0 = time.time()
    m = M.a255()
    c, persons = m.load_muther(False)
    both = persons[(persons.has_fat == "True") & (persons.has_skin == "True")].person
    idx = np.array([i in set(both) for i in c.ids])
    wf = c.wf
    assert wf.genes == U.WF.genes
    cw = C.core(wf)
    n_sets = 3 if chk else m.N_SETS
    n_perm, n_boot = (20, 5) if chk else (2000, M.N_BOOT)
    lev = lambda w: E.radius_covariance(*[c.Z(w)[t][idx] for t in ("fat", "skin")])["leverage"]  # noqa: E731

    # companions (matched tag 3, responsive tag 4), 14 genes and restricted
    comp = {"14": {}, "core": {}}
    for kind, tag in (("matched", 3), ("responsive", 4)):
        ws = c.pseudo(kind, n_sets, tag)
        comp["14"][kind] = np.array([lev(w) for w in ws])
        comp["core"][kind] = np.array([lev(C.restrict(w, base=wf)) for w in ws])
    Z14, Zc = c.Z(wf), c.Z(cw)
    zf14, zs14, zfc, zsc = Z14["fat"][idx], Z14["skin"][idx], Zc["fat"][idx], Zc["skin"][idx]
    L = lambda a, b: E.radius_covariance(a, b)["leverage"]  # noqa: E731
    l14, lc = L(zf14, zs14), L(zfc, zsc)

    # permutation null and twin-pair bootstrap, one rng as score_a331, shared by both waveforms
    rng = np.random.default_rng(M.SEED)
    nul14, nulc = [], []
    for _ in range(n_perm):
        p = rng.permutation(len(zs14))
        nul14.append(L(zf14, zs14[p]))
        nulc.append(L(zfc, zsc[p]))
    ids = np.array(c.ids)[idx]
    pp = persons.set_index("person")
    unit = {}
    for i in ids:
        tw = pp.cotwin.get(i)
        unit[i] = "|".join(sorted([i, tw])) if isinstance(tw, str) and tw in set(ids) else i
    u = pd.Series(unit)
    groups = {k: np.flatnonzero(np.isin(ids, v.index)) for k, v in u.groupby(u)}
    keys = list(groups)
    b14, bc = [], []
    for _ in range(n_boot):
        pick = np.concatenate([groups[k] for k in rng.choice(keys, len(keys))])
        b14.append(L(zf14[pick], zs14[pick]))
        bc.append(L(zfc[pick], zsc[pick]))
    b14, bc = np.array(b14), np.array(bc)
    stat = lambda l, nul, mt: dict(p=float((1 + np.sum(np.array(nul) >= l)) / (len(nul) + 1)), q95=float(np.quantile(mt, 0.95)))  # noqa: E731
    s14 = stat(l14, nul14, comp["14"]["matched"])
    sc_ = stat(lc, nulc, comp["core"]["matched"])

    # ---- 14 genes reproduced first ----------------------------------------------------------------------------------
    if not chk:
        o = pd.read_csv(OUTD / "tests.csv").iloc[0]
        rep = {"leverage": (l14, o.leverage), "p_perm": (s14["p"], o.p_perm), "matched_q95": (s14["q95"], o.matched_q95),
               "responsive_q95": (float(np.quantile(comp["14"]["responsive"], 0.95)), o.responsive_q95),
               "twin lo90": (float(np.quantile(b14, 0.05)), o.lo90_twin_pairs), "twin hi90": (float(np.quantile(b14, 0.95)), o.hi90_twin_pairs)}
        bad = {k: v for k, v in rep.items() if not np.isclose(*v, atol=1e-6, rtol=0)}
        print("14-gene reproduction:", "OK" if not bad else bad, flush=True)
        if bad:
            C.write_rows("a331", [dict(id="A331_H331.1", figure_panel="Fig 4f", value_14_stored=o.leverage, value_14=l14, reproduced=False,
                                        notes=f"not reproduced (1e-6): {bad}", script=SCRIPT)])
            return
        ref = pd.read_csv(U.EXP / "A335_review_batch/output/a335_1c_tests.csv").set_index("waveform").loc["core"]
        assert np.isclose(lc, ref.leverage, atol=1e-9) and np.isclose(sc_["p"], ref.p_perm, atol=1e-9), "A335 1c core not reproduced"

    # ---- gene drops, trio alone -----------------------------------------------------------------------------------------
    dl, tl = C.drops(wf)
    if chk:
        dl = dl[:4] + [x for x in dl if x[0] == tl]
    wfs = dict(dl)
    vals = {lab: lev(w) for lab, w in dl}
    assert np.isclose(vals[tl], lc, atol=1e-9)
    dp = C.drop_percentiles(vals, tl, wfs, +1)
    tw_ = C.trio(wf)
    gates = [P.reference_gate(c.t[k][list(tw_.genes)].to_numpy(float), tw_, 1000 if not chk else 50, np.random.default_rng([M.SEED, 1, i]))["p"]
             for i, k in enumerate(("fat", "skin"))]
    gp = bool(max(gates) < 0.05)
    v_trio = lev(tw_)
    st = U.status([lc > sc_["q95"], sc_["p"] < 0.05])
    lo, hi = U.ci90(b14 - bc)
    if chk:
        print("check ok")
        return
    mc = comp["core"]
    row = dict(id="A331_H331.1", figure_panel="Fig 4f", statistic="MuTHER leverage of fat and skin clock-plane radii (680 women)",
               direction=1, value_14_stored=o.leverage, value_14=l14, reproduced=True, value_core=lc,
               core_criterion="leverage > core matched q95 AND p vs persons-permuted null < 0.05", core_status=st,
               diff_14_minus_core=l14 - lc, diff_lo90=lo, diff_hi90=hi,
               diff_draws="2,000 twin-pair-unit bootstrap draws (score_a331), shared by both waveforms",
               trio_gate_passes=gp, trio_alone_value=v_trio if gp else np.nan, reading=C.reading(st, dp["trio_drop_below_q10"]),
               blind="core not blind (A335 / A329)", script=SCRIPT,
               notes=(f"core p_perm {sc_['p']:.5f}, matched q95 {sc_['q95']:.4f} (pct matched {np.mean(mc['matched'] < lc) * 100:.1f}), "
                      f"responsive q95 {np.quantile(mc['responsive'], 0.95):.4f} (pct {np.mean(mc['responsive'] < lc) * 100:.1f}); "
                      f"trio reference-gate p fat {gates[0]:.4f}, skin {gates[1]:.4f} (A255 gate, 1,000 permutations)"),
               **{x: dp[x] for x in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")})
    df = C.write_rows("a331", [row])
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
