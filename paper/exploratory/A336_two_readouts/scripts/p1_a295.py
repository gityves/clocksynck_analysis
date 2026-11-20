"""A336 part 1 wrapper: A295 H295.1 / H295.2 (composition does not change the gain, nor the phase, of a tissue's clock
state; GTEx fast deaths, 27 peripheral tissues, random-effects pooling) on the 14 genes, the core loop, the trio
contribution, the 364 gene drops and the trio alone.

Original: A295_composition/scripts/score_a295.py (imported: `fit`, SCORES, BOUND_H, `A292.prepared`'s steps). No
companions (marker-score interaction). Each waveform is one task: Dataset fit (`A104.per_tissue`), per tissue and score
the interaction fit with its 200 donor-bootstrap draws (rng [295, tissue, score], as the original). Trio contribution:
the same draws for both waveforms, pooled estimate per draw with each waveform's own tissue SEs (random-effects
meta). Gene drops: pooled gain change of the lymphocyte score (direction -1, original negative) and, for H295.2, the
larger of the two scores' 90 % half-extents max(|lo|, |hi|) (direction -1: tighter = stronger). Trio alone: the gate on
the same rows and the trio's pooled values. Aggregates only. Long run: `RUN_CAP_GB=30 paper/common/run.sh <this> --jobs 4`.

    python p1_a295.py [--jobs 4] [--check]
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
from scipy.stats import norm
sys.path.insert(0, str(U.PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync import stats as ST
from clocksync.config import limit_threads
from common import gtex

spec = importlib.util.spec_from_file_location("a295", U.EXP / "A295_composition/scripts/score_a295.py")
A295 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A295)
A292 = A295.A292
SCORES, BOUND_H = A295.SCORES, A295.BOUND_H
META = ["donor", "compartment", "depth"]
OUTD = U.EXP / "A295_composition" / "output"
SCRIPT = str(U.HERE.parent / "p1_a295.py")
Gl = {}


def evaluate(i):
    """Per tissue and score: gain and phase estimates, and the n_boot bootstrap draws (NaN where the fit failed)."""
    warnings.filterwarnings("ignore")
    limit_threads(1)
    lab, w = Gl["sets"][i]
    pb = Gl["pb"]
    try:
        comps0 = sorted(pb.compartment.unique())
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=comps0, waveform=w, name=lab)
        pt = A104.per_tissue(ds)
        comps = [k for k in ds.anchors if k in pt.index]
        D = A292.index_donors(A292.designs(ds, Gl["tod"]), comps)
        rows = {}
        for k in comps:
            e = D[k]
            for name in SCORES:
                s = Gl["sc"][name].reindex(list(zip(e["donor"], [k] * len(e["donor"])))).to_numpy(float)
                ok = np.isfinite(s)
                x, h, G, s_ = e["x"][ok], e["h"][ok], e["G"][ok], s[ok]
                g0, p0 = A295.fit(x, h, s_, G)
                rng = np.random.default_rng([295, comps.index(k), list(SCORES).index(name)])
                bs = np.full((Gl["n_boot"], 2), np.nan)
                for b in range(Gl["n_boot"]):
                    ii = rng.integers(0, len(x), len(x))
                    try:
                        bs[b] = A295.fit(x[ii], h[ii], s_[ii], G[ii])
                    except np.linalg.LinAlgError:
                        continue
                rows[(k, name)] = (g0, p0, bs)
        return lab, rows
    except Exception as e:  # noqa: BLE001
        return lab, f"{type(e).__name__}: {e}"[:200]


def pooled(rows):
    """Per score: pooled gain change, its SE and one-sided p, pooled phase change and its 90% interval (as score_a295)."""
    out = {}
    for name in SCORES:
        ks = [k for (k, n) in rows if n == name]
        g = np.array([rows[(k, name)][0] for k in ks])
        p = np.array([rows[(k, name)][1] for k in ks])
        gse = np.array([np.nanstd(rows[(k, name)][2][:, 0], ddof=1) for k in ks])
        pse = np.array([np.nanstd(rows[(k, name)][2][:, 1], ddof=1) for k in ks])
        mg, mp = ST.meta(g, gse), ST.meta(p, pse)
        lo, hi = mp["random"] - 1.645 * mp["se_random"], mp["random"] + 1.645 * mp["se_random"]
        out[name] = dict(gain=mg["random"], gain_se=mg["se_random"], p=float(norm.cdf(mg["random"] / mg["se_random"])),
                         phase=mp["random"], lo=lo, hi=hi, ks=ks, gse=gse, pse=pse)
    return out


def pooled_draws(rows, P_):
    """Pooled estimate per bootstrap draw with the waveform's own tissue SEs: {score: (gain (n_boot,), phase (n_boot,))}."""
    out = {}
    for name in SCORES:
        ks = P_[name]["ks"]
        nb = rows[(ks[0], name)][2].shape[0]
        G_ = np.array([[rows[(k, name)][2][b, 0] for k in ks] for b in range(nb)])
        P2 = np.array([[rows[(k, name)][2][b, 1] for k in ks] for b in range(nb)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out[name] = (np.array([ST.meta(g, P_[name]["gse"])["random"] for g in G_]),
                         np.array([ST.meta(p, P_[name]["pse"])["random"] for p in P2]))
    return out


def extent(P_):
    return float(max(max(abs(P_[n]["lo"]), abs(P_[n]["hi"])) for n in SCORES))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    chk = a.check
    t0 = time.time()
    warnings.filterwarnings("ignore")
    markers = sorted({g for v in SCORES.values() for g in v})
    clock = list(cs.gtex().genes)
    pbm, donors = gtex.load(min_donors=A104.MIN_DONORS, genes=clock + markers)
    pbm = pbm[~pbm.compartment.str.startswith("Brain")].reset_index(drop=True)
    sc = {}
    for name, g in SCORES.items():
        v = np.log((pbm[g].sum(1) + 0.5) / pbm.depth)
        sc[name] = v.groupby(pbm.compartment).transform(lambda y: (y - y.mean()) / y.std())
    sc = pd.DataFrame(sc).set_index([pbm.donor.astype(str), pbm.compartment])
    pb = pbm[META + clock]
    limit_threads(1)
    tw, wfc = C.trio(), C.core()
    dl, tl = C.drops()
    if chk:
        dl = dl[:2] + [x for x in dl if x[0] == tl]
    sets = [("clock14", U.WF), ("core", wfc), ("trio", tw)] + dl
    Gl.update(pb=pb, tod=donors.tod, sc=sc, sets=sets, n_boot=5 if chk else 200)

    # ---- 14 genes first -------------------------------------------------------------------------------------------
    lab, r14 = evaluate(0)
    assert isinstance(r14, dict), r14
    P14 = pooled(r14)
    ref = pd.read_csv(OUTD / "tests.csv").set_index("score")
    bad = {}
    for n in SCORES:
        for k, c in (("gain", "pooled_gain_change"), ("gain_se", "gain_se"), ("p", "gain_p_one_sided_negative"),
                     ("phase", "pooled_phase_change_h"), ("lo", "phase_lo90"), ("hi", "phase_hi90")):
            if not np.isclose(P14[n][k], ref.loc[n, c], atol=1e-6, rtol=0):
                bad[(n, k)] = (P14[n][k], ref.loc[n, c])
    if bad and not chk:
        print("14-gene NOT reproduced", bad, flush=True)
        C.write_rows("a295", [dict(id=i, figure_panel="Fig 4b", value_14_stored=v, value_14=w_, reproduced=False,
                                    notes=f"14-gene values not reproduced (1e-6): {bad}", script=SCRIPT)
                              for i, v, w_ in (("A295_H295.1", ref.loc["lymphocyte", "pooled_gain_change"], P14["lymphocyte"]["gain"]),
                                               ("A295_H295.2", ref.loc["lymphocyte", "pooled_phase_change_h"], P14["lymphocyte"]["phase"]))])
        return
    print("14-gene reproduction OK", flush=True)

    # ---- core, trio, drops (all other waveforms) ---------------------------------------------------------------------
    res = {}
    with Pool(a.jobs) as pool:
        for k, (lab, r) in enumerate(pool.imap(evaluate, range(1, len(sets)), chunksize=1)):
            res[lab] = r
            if k % 10 == 0:
                print(f"{k + 1}/{len(sets) - 1} waveforms {time.time() - t0:.0f}s", flush=True)
    res["clock14"] = r14
    errs = {k: v for k, v in res.items() if not isinstance(v, dict)}
    Pc = pooled(res["core"])
    assert isinstance(res[tl], dict) and np.isclose(pooled(res[tl])["lymphocyte"]["gain"], Pc["lymphocyte"]["gain"], atol=1e-9)
    D14, Dc = pooled_draws(r14, P14), pooled_draws(res["core"], Pc)
    P_trio = pooled(res["trio"]) if isinstance(res["trio"], dict) else None
    # trio alone gate
    try:
        comps = sorted(pb.compartment.unique())
        ds = cs.Dataset(pb[META + list(tw.genes)], anchors=comps, waveform=tw, name="trio")
        g = ds.gate(n_perm=1000 if not chk else 50, seed=336)
        gp = (bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), float(g.p_pairing), float(g.p_acrophase))
    except Exception as e:  # noqa: BLE001
        gp = ("not computable", np.nan, np.nan)
        print("gate error", e, flush=True)

    wfs = dict(dl)
    v1, v2 = {}, {}
    for lab, _ in dl:
        r = res[lab]
        if isinstance(r, dict):
            P_ = pooled(r)
            v1[lab], v2[lab] = P_["lymphocyte"]["gain"], extent(P_)
        else:
            v1[lab] = v2[lab] = np.nan
    if chk:
        C.drop_percentiles(v1, tl, wfs, -1)
        assert len(D14["lymphocyte"][0]) == Gl["n_boot"]
        print("check ok")
        return

    def passes1(P_):
        return P_["lymphocyte"]["gain"] < 0 and P_["lymphocyte"]["p"] < 0.05

    def inside(P_):
        return all(-BOUND_H < P_[n]["lo"] and P_[n]["hi"] < BOUND_H for n in SCORES)

    out = []
    dp1, dp2 = C.drop_percentiles(v1, tl, wfs, -1), C.drop_percentiles(v2, tl, wfs, -1)
    lo1, hi1 = U.ci90(D14["lymphocyte"][0] - Dc["lymphocyte"][0])
    d2 = np.maximum.reduce([np.abs(D14[n][1]) for n in SCORES]) - np.maximum.reduce([np.abs(Dc[n][1]) for n in SCORES])
    lo2, hi2 = U.ci90(d2)
    common = dict(figure_panel="Fig 4b", direction=-1, reproduced=True, blind="core blind", script=SCRIPT,
                  diff_draws="200 donor bootstraps per tissue and score (original rng [295, tissue, score]), pooled per draw with each waveform's own tissue SEs",
                  trio_gate_passes=gp[0])
    # H295.1: original verdict FAIL; core verdict compared
    v14_pass = passes1(P14)
    core_pass = passes1(Pc)
    st1 = "holds" if core_pass == v14_pass else "falls" if v14_pass else "at the line"
    out.append(dict(id="A295_H295.1", statistic="pooled lymphocyte-score interaction on gain (Re(b/a) per SD), 27 tissues, random effects",
                    value_14_stored=ref.loc["lymphocyte", "pooled_gain_change"], value_14=P14["lymphocyte"]["gain"],
                    value_core=Pc["lymphocyte"]["gain"],
                    core_criterion="pooled gain change < 0 AND one-sided p < 0.05; status = verdict unchanged from the 14-gene one (FAIL: no gain change)",
                    core_status=st1, diff_14_minus_core=P14["lymphocyte"]["gain"] - Pc["lymphocyte"]["gain"], diff_lo90=lo1, diff_hi90=hi1,
                    trio_alone_value=P_trio["lymphocyte"]["gain"] if gp[0] is True and P_trio else np.nan,
                    reading=C.reading(st1, dp1["trio_drop_below_q10"]),
                    notes=(f"14-gene verdict {'PASS' if v14_pass else 'FAIL'} (pooled {P14['lymphocyte']['gain']:.5f}, p {P14['lymphocyte']['p']:.4f}); "
                           f"core pooled {Pc['lymphocyte']['gain']:.5f} (se {Pc['lymphocyte']['gain_se']:.5f}), p {Pc['lymphocyte']['p']:.4f}: verdict {'PASS' if core_pass else 'FAIL'}; "
                           f"gate p_pairing {gp[1]:.4f}, p_acrophase {gp[2]:.4f}; waveform fit errors {len(errs)}; drops: pooled lymphocyte gain, direction -1"),
                    **common, **{x: dp1[x] for x in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")}))
    ok14, okc = inside(P14), inside(Pc)
    st2 = "holds" if okc else "falls"
    out.append(dict(id="A295_H295.2", statistic="larger 90% half-extent max(|lo|, |hi|) (h) of the pooled phase change per SD, both scores",
                    value_14_stored=max(max(abs(ref.loc[n, "phase_lo90"]), abs(ref.loc[n, "phase_hi90"])) for n in SCORES), value_14=extent(P14),
                    value_core=extent(Pc),
                    core_criterion=f"90% interval of pooled phase change inside +-{BOUND_H} h for both scores",
                    core_status=st2, diff_14_minus_core=extent(P14) - extent(Pc), diff_lo90=lo2, diff_hi90=hi2,
                    trio_alone_value=extent(P_trio) if gp[0] is True and P_trio else np.nan,
                    reading=C.reading(st2, dp2["trio_drop_below_q10"]),
                    notes=("core intervals: " + "; ".join(f"{n} {Pc[n]['lo']:.4f} to {Pc[n]['hi']:.4f} (pooled {Pc[n]['phase']:.4f})" for n in SCORES) +
                           f"; 14 genes within bound {ok14}; direction -1 (tighter = stronger)"),
                    **common, **{x: dp2[x] for x in ("n_drops", "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median", "trio_drop_below_q10")}))
    df = C.write_rows("a295", out)
    pd.set_option("display.width", 250)
    print(df.T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
