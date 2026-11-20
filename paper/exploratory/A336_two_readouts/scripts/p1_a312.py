"""A336 part 1 wrapper: A312 H312.1 (amplitude elasticity of the 27 peripheral GTEx tissues to the persons' shared clock
amplitude falls with gain: Spearman rho(elasticity, log S), one-sided tissue-label permutation) on the 14 genes, the core
loop, the trio contribution, the gene drops and the trio alone.

Original: A312_elasticity/scripts/score_a312.py (same calls: `E.split_hosts(ds, others, 20, seed=[312, k])`,
`E.from_hosts`; permutation null rng [312, 2], 2,000 permutations). No companions. S = z2 V per tissue is refitted for each
waveform (`A104.per_tissue`; the 14-gene S must equal A294's table). The donor bootstrap of the original only feeds the
elasticity SE and is not run. Trio contribution: 200 tissue bootstraps (the statistic's unit; rng [336, 1]) of the paired
(elasticity, log S), Spearman rho of both waveforms on the same resamples. Gene drops: rho per drop, direction -1 (the
registered effect is negative). Trio alone: rho of the trio if the model-free gate on the periphery rows passes.
Aggregates only. Long run: `RUN_CAP_GB=30 paper/common/run.sh <this> --jobs 4`.

    python p1_a312.py [--jobs 4] [--check] [--n-drops N]
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
from clocksync import elasticity as E
from clocksync.config import limit_threads
from common import gtex

META = ["donor", "compartment", "depth"]
SCRIPT = str(U.HERE.parent / "p1_a312.py")
Gl = {}


def table(i):
    """Elasticity and S per tissue (sorted tissue names) of waveform i."""
    warnings.filterwarnings("ignore")
    limit_threads(1)
    lab, w = Gl["sets"][i]
    pb = Gl["pb"]
    try:
        comps = sorted(pb.compartment.unique())
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=lab)
        pt = A104.per_tissue(ds)
        el = {}
        for k, t in enumerate(comps):
            hosts = E.split_hosts(ds, [c for c in comps if c != t], 20, seed=[312, k])
            el[t] = E.from_hosts(ds, t, hosts)["elasticity"]
        out = pd.DataFrame(dict(elasticity=pd.Series(el), S=(pt.z2 * pt.V))).sort_index()
        return lab, out
    except Exception as e:  # noqa: BLE001
        return lab, f"{type(e).__name__}: {e}"[:200]


def rho_p(t, n=2000):
    """rho(elasticity, log S) and the original one-sided tissue-label permutation p (rng [312, 2])."""
    x, ls = t.elasticity.to_numpy(), np.log(t.S.to_numpy())
    rho = spearmanr(x, ls).statistic
    rng = np.random.default_rng([312, 2])
    null = np.array([spearmanr(rng.permutation(x), ls).statistic for _ in range(n)])
    return float(rho), float((1 + np.sum(null <= rho)) / (n + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--n-drops", type=int, default=None)
    a = ap.parse_args()
    chk = a.check
    t0 = time.time()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    limit_threads(1)
    tw = C.trio()
    dl, tl = C.drops(n=a.n_drops)
    if chk:
        dl = dl[:2] + [x for x in dl if x[0] == tl]
    sets = [("clock14", U.WF), ("trio", tw)] + dl
    Gl.update(pb=pb, sets=sets)

    # ---- 14 genes first -------------------------------------------------------------------------------------------
    lab, t14 = table(0)
    assert isinstance(t14, pd.DataFrame), t14
    ref = pd.read_csv(U.EXP / "A312_elasticity/output/tests.csv").iloc[0]
    S_ref = pd.read_csv(U.EXP / "A294_floor/output/tissues.csv", index_col=0).S.reindex(t14.index)
    r14, p14 = rho_p(t14, 2000)
    rep = {"rho": (r14, ref.rho), "p_one_sided": (p14, ref.p_one_sided)}
    bad = {k: v for k, v in rep.items() if not np.isclose(*v, atol=1e-6, rtol=0)}
    if not np.allclose(t14.S, S_ref, atol=1e-6, rtol=0):
        bad["S"] = "A104 S differs from A294 table"
    if bad and not chk:
        print("14-gene NOT reproduced", bad, flush=True)
        C.write_rows("a312", [dict(id="A312_H312.1", figure_panel="Fig 4f", value_14_stored=ref.rho, value_14=r14, reproduced=False,
                                    notes=f"14-gene values not reproduced (1e-6): {bad}", script=SCRIPT)])
        return
    print("14-gene reproduction OK", flush=True)

    res = {}
    with Pool(a.jobs) as pool:
        for k, (lab, t) in enumerate(pool.imap(table, range(1, len(sets)), chunksize=1)):
            res[lab] = t
            if k % 5 == 0:
                print(f"{k + 1}/{len(sets) - 1} waveforms {time.time() - t0:.0f}s", flush=True)
    errs = {k: v for k, v in res.items() if not isinstance(v, pd.DataFrame)}
    tc = res[tl]
    rc, pc = rho_p(tc, 2000)
    rng = np.random.default_rng(U.SEED_BOOT)
    j = t14.index
    diffs = []
    for _ in range(3 if chk else U.N_BOOT_UNIT):
        i = rng.integers(0, len(j), len(j))
        diffs.append(spearmanr(t14.elasticity.to_numpy()[i], np.log(t14.S.to_numpy())[i]).statistic
                     - spearmanr(tc.elasticity.to_numpy()[i], np.log(tc.S.to_numpy())[i]).statistic)
    vals = {lab: (float(spearmanr(t.elasticity, np.log(t.S)).statistic) if isinstance(t, pd.DataFrame) else np.nan)
            for lab, t in res.items() if lab in dict(dl)}
    assert np.isclose(vals[tl], rc, atol=1e-9)
    dp = C.drop_percentiles(vals, tl, dict(dl), -1)
    try:
        g = cs.Dataset(pb[META + list(tw.genes)], anchors=sorted(pb.compartment.unique()), waveform=tw, name="trio").gate(
            n_perm=1000 if not chk else 50, seed=336)
        gp = (bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), float(g.p_pairing), float(g.p_acrophase))
    except Exception as e:  # noqa: BLE001
        gp = ("not computable", np.nan, np.nan)
        print("gate error", e, flush=True)
    tt = res["trio"]
    v_trio = float(spearmanr(tt.elasticity, np.log(tt.S)).statistic) if isinstance(tt, pd.DataFrame) else np.nan
    if chk:
        print("check ok")
        return
    st = U.status([rc < 0, pc < 0.05])
    lo, hi = U.ci90(diffs)
    row = dict(id="A312_H312.1", figure_panel="Fig 4f", statistic="Spearman rho(elasticity, S) over 27 tissues", direction=-1,
               value_14_stored=ref.rho, value_14=r14, reproduced=True, value_core=rc,
               core_criterion="rho < 0 AND one-sided tissue-label permutation p < 0.05 (2,000)", core_status=st,
               diff_14_minus_core=r14 - rc, diff_lo90=lo, diff_hi90=hi,
               diff_draws="200 tissue bootstraps of paired (elasticity, log S), seed [336,1]", trio_gate_passes=gp[0],
               trio_alone_value=v_trio if gp[0] is True else np.nan, reading=C.reading(st, dp["trio_drop_below_q10"]),
               blind="core blind", script=SCRIPT, notes=(f"core p {pc:.5f}; gate p_pairing {gp[1]:.4f}, p_acrophase {gp[2]:.4f}; "
                                                          f"waveform fit errors {len(errs)}; drops n={len(dl)} direction -1"),
               **{x: dp[x] for x in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")})
    df = C.write_rows("a312", [row])
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
