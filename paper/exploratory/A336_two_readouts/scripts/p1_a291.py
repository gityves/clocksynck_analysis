"""A336 part 1 wrapper: A291 H291.1 (is "stronger clocks are more coupled" what measurement alone produces: Spearman rho of
own amplitude with coupling / simulated ceiling over the 27 peripheral GTEx tissues) on the 14 genes, the core loop, the
trio contribution, the 364 gene drops and the trio alone.

Original: A291_attenuation/scripts/score_a291.py (imported: `simulate_tissue`, `perm_rho`, seeds). No companions. One
waveform = one task: Dataset (fast deaths, periphery, every tissue an anchor), `A104.per_tissue`, and per tissue the
ceiling simulation `simulate_tissue` (100 replicates at truth 1 and at the estimate, original seeds). Trio contribution:
200 tissue bootstraps (the statistic's unit; rng [336, 1]) of the paired per-tissue (A, c / c_ceil), Spearman rho of
both waveforms on the same resamples. Gene drops: rho(A, c / c_ceil) per drop (the permutation p is not recomputed),
direction +1. Trio alone: the trio's rho, if the model-free gate of the trio on the periphery rows passes.
Aggregates only. Long run: `RUN_CAP_GB=30 paper/common/run.sh <this> --jobs 4`.

    python p1_a291.py [--jobs 4] [--check] [--n-drops N]
"""
import argparse
import importlib.util
import sys
import time
import warnings
from multiprocessing import Pool

from p1_util import *  # noqa: F401,F403
import p1_util as U
import a336_common as C
from scipy.stats import spearmanr
sys.path.insert(0, str(U.PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads
from common import gtex

spec = importlib.util.spec_from_file_location("a291", U.EXP / "A291_attenuation/scripts/score_a291.py")
A291 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A291)
META = ["donor", "compartment", "depth"]
SCRIPT = str(U.HERE.parent / "p1_a291.py")
Gl = {}


def table(i, n_rep=None):
    """Per tissue: own amplitude, coupling, ceiling coupling and c / c_ceil of waveform i (NaN table on failure)."""
    warnings.filterwarnings("ignore")
    limit_threads(1)
    lab, w = Gl["sets"][i]
    pb = Gl["pb"]
    try:
        comps = sorted(c for c in pb.compartment.unique() if not c.startswith("Brain"))
        ds = cs.Dataset(pb[pb.compartment.isin(comps)][META + list(w.genes)], anchors=comps, waveform=w, name=lab)
        obs = A104.per_tissue(ds)
        rows = [A291.simulate_tissue(ds, "periphery", k, j, n_rep or Gl["n_rep"]) for j, k in enumerate(obs.index)]
        t = obs.join(pd.DataFrame(rows).set_index("tissue")[["c_ceil"]])
        t["c_over_ceil"] = t.coupling / t.c_ceil
        return lab, t[["own_amplitude", "coupling", "c_ceil", "c_over_ceil"]], float(A104.rho(obs))
    except Exception as e:  # noqa: BLE001
        return lab, f"{type(e).__name__}: {e}"[:200], np.nan


def rho_of(t):
    ok = np.isfinite(t.own_amplitude) & np.isfinite(t.c_over_ceil)
    return float(spearmanr(t.own_amplitude[ok], t.c_over_ceil[ok]).statistic)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--n-drops", type=int, default=None)
    a = ap.parse_args()
    chk = a.check
    t0 = time.time()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    limit_threads(1)
    tw = C.trio()
    dl, tl = C.drops(n=a.n_drops)
    if chk:
        dl = dl[:2] + [x for x in dl if x[0] == tl]
    sets = [("clock14", U.WF), ("trio", tw)] + dl
    Gl.update(pb=pb, sets=sets, n_rep=5 if chk else 100)

    # ---- 14 genes first -------------------------------------------------------------------------------------------
    lab, t14, rho_obs = table(0)
    assert isinstance(t14, pd.DataFrame), t14
    r14, p14, n14 = A291.perm_rho(t14.own_amplitude.to_numpy(float), t14.c_over_ceil.to_numpy(float), A291.SEED)
    ref = pd.read_csv(U.EXP / "A291_attenuation/output/tests.csv")
    ref = ref[ref.block == "periphery"].iloc[0]
    rep = {"rho_obs": (rho_obs, ref.rho_obs), "rho_corrected": (r14, ref.rho_corrected), "p_one_sided": (p14, ref.p_one_sided),
           "median_c_ceil": (float(t14.c_ceil.median()), ref.median_c_ceil)}
    bad = {k: v for k, v in rep.items() if not np.isclose(*v, atol=1e-6, rtol=0)}
    if bad and not chk:
        print("14-gene NOT reproduced", bad, flush=True)
        C.write_rows("a291", [dict(id="A291_H291.1", figure_panel="Fig 4e", value_14_stored=ref.rho_corrected, value_14=r14, reproduced=False,
                                    notes=f"14-gene values not reproduced (1e-6): {bad}", script=SCRIPT)])
        return
    print("14-gene reproduction OK", flush=True)

    # ---- core (the trio drop), trio, drops ----------------------------------------------------------------------------
    res = {}
    with Pool(a.jobs) as pool:
        for k, (lab, t, _) in enumerate(pool.imap(table, range(1, len(sets)), chunksize=1)):
            res[lab] = t
            if k % 10 == 0:
                print(f"{k + 1}/{len(sets) - 1} waveforms {time.time() - t0:.0f}s", flush=True)
    errs = {k: v for k, v in res.items() if not isinstance(v, pd.DataFrame)}
    tc = res[tl]
    rc, pc, _ = A291.perm_rho(tc.own_amplitude.to_numpy(float), tc.c_over_ceil.to_numpy(float), A291.SEED)
    # trio contribution: 200 tissue bootstraps, paired
    j = t14.index.intersection(tc.index)
    A14, c14_, Ac, cc_ = (t14.loc[j].own_amplitude.to_numpy(), t14.loc[j].c_over_ceil.to_numpy(),
                          tc.loc[j].own_amplitude.to_numpy(), tc.loc[j].c_over_ceil.to_numpy())
    rng = np.random.default_rng(U.SEED_BOOT)
    diffs = []
    for _ in range(3 if chk else U.N_BOOT_UNIT):
        i = rng.integers(0, len(j), len(j))
        diffs.append(spearmanr(A14[i], c14_[i]).statistic - spearmanr(Ac[i], cc_[i]).statistic)
    vals = {lab: (rho_of(t) if isinstance(t, pd.DataFrame) else np.nan) for lab, t in res.items() if lab in dict(dl)}
    assert np.isclose(vals[tl], rc, atol=1e-9)
    dp = C.drop_percentiles(vals, tl, dict(dl), +1)
    try:
        comps = sorted(c for c in pb.compartment.unique() if not c.startswith("Brain"))
        ds = cs.Dataset(pb[pb.compartment.isin(comps)][META + list(tw.genes)], anchors=comps, waveform=tw, name="trio")
        g = ds.gate(n_perm=1000 if not chk else 50, seed=336)
        gp = (bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), float(g.p_pairing), float(g.p_acrophase))
    except Exception as e:  # noqa: BLE001
        gp = ("not computable", np.nan, np.nan)
        print("gate error", e, flush=True)
    tt = res["trio"]
    v_trio = rho_of(tt) if isinstance(tt, pd.DataFrame) else np.nan
    if chk:
        print("check ok")
        return
    st = U.status([rc > 0, pc < 0.05])
    lo, hi = U.ci90(diffs)
    row = dict(id="A291_H291.1", figure_panel="Fig 4e", statistic="Spearman rho(own amplitude, coupling / simulated ceiling), 27 periphery tissues",
               direction=1, value_14_stored=ref.rho_corrected, value_14=r14, reproduced=True, value_core=rc,
               core_criterion="rho_corrected > 0 AND one-sided permutation p < 0.05 (10,000)", core_status=st,
               diff_14_minus_core=r14 - rc, diff_lo90=lo, diff_hi90=hi,
               diff_draws="200 tissue bootstraps of paired (A, c/c_ceil), seed [336,1]", trio_gate_passes=gp[0],
               trio_alone_value=v_trio if gp[0] is True else np.nan, reading=C.reading(st, dp["trio_drop_below_q10"]),
               blind="core blind", script=SCRIPT,
               notes=(f"core p {pc:.5f}, median c_ceil {float(tc.c_ceil.median()):.4f}; 14-gene rho_obs reproduced; gate p_pairing {gp[1]:.4f}, "
                      f"p_acrophase {gp[2]:.4f}; waveform fit errors {len(errs)}; drops: rho_corrected without permutation p"),
               **{x: dp[x] for x in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")})
    df = C.write_rows("a291", [row])
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
