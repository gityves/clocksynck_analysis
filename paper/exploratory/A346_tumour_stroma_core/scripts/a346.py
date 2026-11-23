"""A346: tumour stroma pairing with the patient's normal stroma on 14 genes and on the core loop, calibrated companions.

Registered in paper/docs/registered_experiments.md (A346). Units: A227 cohorts and tables, tumour vs adjacent-normal
fibroblast and endothelial (primary); A149 cohorts, malignant (or tumour epithelium) vs adjacent-normal epithelium (D346,
descriptive). Statistic: per cohort the uncorrected |complex correlation| of the stratum-centred states (A227 `raw_stat`),
pooled with patient-count weights. Every unit holds one stacked state matrix (rows = gene sets):

    [0, 601)       14 genes: clock, 200 matched, 200 narrow responsive, 200 broad responsive (A227 / A285 / A320 seeds)
    [601, 1202)    core loop: core clock, the same 600 sets restricted to the core positions (a336_common.restrict)

  power  14-gene reproduction of the stored A227 anchors; on both readouts the within-cohort permutation nulls (10,000,
         shuffled pairings only), the companions on the real pairing and 2,000 placebo draws (normal rows permuted within
         chemistry stratum, seeds [346, k, d]); the identity-pairing statistic of the core clock row is NOT computed.
         p_cal that hypothetical core values (14-gene value - delta) would obtain.   -> output/power/*.csv, output/cache/
  score  identity-pairing statistics on both readouts, pairing p, p_cal per family, Holm over the two compartments, V346,
         H346.1-3, D346, trio contribution (200 patient bootstraps within cohort).   -> output/score/*.csv

    python a346.py --stage power|score [--jobs 2] [--smoke]     (long: paper/common/run.sh, RUN_CAP_GB=40)
"""
from __future__ import annotations

import argparse
import pickle
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(PAPER))
for p in ("A227_tme_compartment_pairing/scripts", "A149_paired_tumour_normal/scripts", "A151_fibroblast_pairing/scripts",
          "A226_stroma_pairing_confirm/scripts", "A336_two_readouts/scripts"):
    sys.path.insert(0, str(EXP / p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a149_lib as L  # noqa: E402
import a336_common as C  # noqa: E402
import compute_a227 as M  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.pairwise import calibrated_p, percentile  # noqa: E402

OUT = HERE.parents[1] / "output"
BASE = C.WF
CORE = C.core(BASE)
N_SET, N_PERM, N_DRAW, N_BOOT = 200, 10_000, 2000, 200
FAMS = {"matched": [44], "responsive": [44], "responsive_broad": [45]}  # seed tags (A227 / A285 / A320)
KCODE = {"fibroblast": 1, "endothelial": 2, "epithelium": 3}
PRIMARY = ("fibroblast", "endothelial")
NB = 1 + 3 * N_SET                       # rows per readout block
BLOCK = {"14": 0, "core": NB}
FAM_SL = {f: slice(1 + i * N_SET, 1 + (i + 1) * N_SET) for i, f in enumerate(FAMS)}
DELTAS = (0.05, 0.10, 0.15)
BAND = 0.0646


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


# ---- units -----------------------------------------------------------------------------------------------------

def _states(G, a, b, pts, n, wfs):
    za, zb = [], []
    for w in wfs:
        p = L.pair_states(G, a, b, wf=w, pts=pts)
        za.append(p["za"] if p else np.full(n, np.nan + 0j))
        zb.append(p["zb"] if p else np.full(n, np.nan + 0j))
    return za, zb


def unit(study, kind):
    """Stacked states (14 + core blocks), strata and permutation nulls of one cohort x compartment, or None."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    from clocksync.paired import pseudo_waveforms
    if kind == "epithelium":
        try:
            pb, mal = L.table(study)
        except Exception as e:  # noqa: BLE001  (as compute_a149: a cohort that cannot be built is reported, not fatal)
            print(f"epithelium unit {study} not built: {type(e).__name__}: {e}", flush=True)
            return None
        a, b = mal, L.NORM
        if a not in set(pb.compartment) or b not in set(pb.compartment):
            return None
        u = L.pair_states(pb, a, b)
        if u is None:
            return None
        G = L.genome_epithelium(study, mal, u["pts"])
        seed = [346, KCODE[kind], sum(map(ord, study))]
    else:
        pb = M.table(study)
        a, b = M.PAIRS[kind]
        if a not in set(pb.compartment) or b not in set(pb.compartment):
            return None
        u = L.pair_states(pb, a, b)
        if u is None:
            return None
        G = M.genome_pair(study, (a, b), u["pts"])
        seed = M.seed_of(kind, study)
    pts, n = u["pts"], u["n"]
    uc = L.pair_states(pb, a, b, wf=CORE, pts=pts)
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    sets = [w for f, tag in FAMS.items() for w in pseudo_waveforms(lm, N_SET, np.random.default_rng(seed + tag), BASE, f)]
    za14, zb14 = _states(G, a, b, pts, n, sets)
    zac, zbc = _states(G, a, b, pts, n, [C.restrict(w, BASE, CORE) for w in sets])
    za = np.array([u["za"]] + za14 + [uc["za"]] + zac)
    zb = np.array([u["zb"]] + zb14 + [uc["zb"]] + zbc)
    null = {r: M.perm_null_raw(dict(za=za[BLOCK[r]], zb=zb[BLOCK[r]]), N_PERM, np.random.default_rng(seed + [1]))
            for r in BLOCK}
    return dict(study=study, kind=kind, n=n, za=za, zb=zb, stratum=np.asarray(u["stratum"]).astype(str), null=null)


def _raw(za, zb):
    za, zb = za - za.mean(1, keepdims=True), zb - zb.mean(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(np.sum(za * np.conj(zb), 1)) / np.sqrt(np.sum(np.abs(za) ** 2, 1) * np.sum(np.abs(zb) ** 2, 1))


def pooled(units, rows, rng=None):
    """Pooled statistic (patient weights, NaN sets drop their cohort) of the given rows; rng: one placebo draw."""
    num, den = 0.0, 0.0
    for u in units:
        zb = u["zb"][rows]
        if rng is not None:
            idx = np.arange(u["n"])
            for st in np.unique(u["stratum"]):
                j = np.flatnonzero(u["stratum"] == st)
                idx[j] = j[rng.permutation(len(j))]
            zb = zb[:, idx]
        r = _raw(u["za"][rows], zb)
        ok = np.isfinite(r)
        num = num + np.where(ok, r, 0.0) * u["n"]
        den = den + ok * u["n"]
    return num / den


def draws(units, kind, rows):
    return np.array([pooled(units, rows, np.random.default_rng([346, KCODE[kind], d])) for d in range(1, N_DRAW + 1)])


# ---- stages ----------------------------------------------------------------------------------------------------

def build(args):
    from joblib import Parallel, delayed
    kinds = list(KCODE)
    a227 = M.cohort_list(args.smoke)
    a149 = L.cohorts()[:2] if args.smoke else L.cohorts()
    keys = [(s, k) for k in PRIMARY for s in a227] + [(s, "epithelium") for s in a149]
    log("units to build", len(keys))
    res = Parallel(n_jobs=args.jobs, return_as="generator_unordered")(delayed(unit)(s, k) for s, k in keys)
    U = {k: [] for k in kinds}
    for i, r in enumerate(res):
        if r is not None:
            U[r["kind"]].append(r)
        if (i + 1) % 10 == 0:
            log("built", i + 1, "of", len(keys))
    for k in kinds:
        U[k].sort(key=lambda u: u["study"])
    return U


def stored_anchor(kind):
    """Stored 14-gene A227 anchor rows (matched and narrow responsive)."""
    d = EXP / "A227_tme_compartment_pairing/output"
    m = pd.read_csv(d / "descriptive.csv").query("role == 'anchor' and kind == @kind").iloc[0]
    r = pd.read_csv(d / "responsive/descriptive.csv").query("role == 'anchor' and kind == @kind").iloc[0]
    return m, r


def stage_power(args):
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out = OUT / ("power_smoke" if args.smoke else "power")
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    U = build(args)
    with open(cache / "units.pkl", "wb") as f:
        pickle.dump(U, f)
    rep, pw, dr = [], [], {}
    for kind, units in U.items():
        if not units:
            continue
        w = np.array([u["n"] for u in units], float)
        i14 = np.arange(0, NB)
        icomp = np.arange(NB + 1, 2 * NB)                      # core companions only: core clock row excluded
        obs14 = pooled(units, i14)
        obsc_comp = pooled(units, icomp)
        D14 = draws(units, kind, i14)
        Dc = draws(units, kind, np.arange(NB, 2 * NB))           # placebo draws of the core block (no identity draw)
        np.savez_compressed(cache / f"draws_{kind}.npz", d14=D14, dcore=Dc, obs14=obs14, obsc_comp=obsc_comp)
        # 14 genes: reproduction of the stored A227 anchor (not for the epithelium: new seeds, descriptive)
        nul14 = M.pool_null(w, [u["null"]["14"] for u in units])
        r = dict(kind=kind, cohorts=len(units), patients=int(w.sum()), clock_14=obs14[0],
                 p_pairing_14=float((1 + np.sum(nul14 >= obs14[0])) / (1 + len(nul14))),
                 **{f"{f}_q95_14": float(np.nanquantile(obs14[sl], .95)) for f, sl in FAM_SL.items()},
                 **{f"{f}_median_14": float(np.nanmedian(obs14[sl])) for f, sl in FAM_SL.items()})
        if kind in PRIMARY:
            m, s = stored_anchor(kind)
            r.update(stored_clock=m.statistic, stored_patients=m.patients, stored_matched_q95=m.pseudo_q95,
                     stored_responsive_q95=s.pseudo_q95,
                     reproduced=bool(np.isclose(obs14[0], m.statistic, atol=1e-6) and int(w.sum()) == m.patients
                                     and np.isclose(r["matched_q95_14"], m.pseudo_q95, atol=1e-6)
                                     and np.isclose(r["responsive_q95_14"], s.pseudo_q95, atol=1e-6)))
        rep.append(r)
        # core: null SD and MDE, companions on the real pairing, p_cal of hypothetical values
        nulc = M.pool_null(w, [u["null"]["core"] for u in units])
        comp = np.r_[np.nan, obsc_comp]                           # column 0 placeholder for the clock
        for delta in DELTAS:
            v = obs14[0] - delta
            vec = comp.copy()
            vec[0] = v
            row = dict(kind=kind, delta=delta, hypothetical_core=v, core_null_sd=float(nulc.std()),
                       core_mde=float(2.49 * nulc.std()), core_p_pairing_hyp=float((1 + np.sum(nulc >= v)) / (1 + len(nulc))))
            for f, sl in FAM_SL.items():
                row[f"{f}_core_q95"] = float(np.nanquantile(comp[sl], .95))
                row[f"{f}_core_median"] = float(np.nanmedian(comp[sl]))
                row[f"{f}_pcal"] = float(calibrated_p(percentile(vec, sl)[0], percentile(Dc, sl))[0])
                row[f"{f}_placebo_rate_core"] = float(np.mean(Dc[:, 0] > np.nanquantile(Dc[:, sl], .95, axis=1)))
            pw.append(row)
    rep = pd.DataFrame(rep)
    pw = pd.DataFrame(pw)
    # Holm over the two primary compartments, per delta and family
    for f in FAMS:
        pw[f"{f}_pcal_holm"] = np.nan
        for delta in DELTAS:
            sel = pw.index[(pw.delta == delta) & pw.kind.isin(PRIMARY)]
            pw.loc[sel, f"{f}_pcal_holm"] = M.holm(pw.loc[sel, f"{f}_pcal"].to_numpy())
    pw["decidable"] = (pw.matched_pcal_holm < 0.05) & (pw.responsive_pcal_holm < 0.05)
    rep.to_csv(out / "repro14.csv", index=False)
    pw.to_csv(out / "power.csv", index=False)
    pd.set_option("display.width", 250)
    log("\n" + rep.round(4).to_string())
    log("\n" + pw.round(4).to_string())


def stage_score(args):
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out = OUT / ("score_smoke" if args.smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    with open(cache / "units.pkl", "rb") as f:
        U = pickle.load(f)
    if args.smoke:  # code check only: the pairing is scrambled so that no real core value is ever read
        rng = np.random.default_rng(0)
        for units in U.values():
            for u in units:
                u["zb"] = u["zb"][:, rng.permutation(u["n"])]
    rows, trio = [], []
    for kind, units in U.items():
        if not units:
            continue
        z = np.load(cache / f"draws_{kind}.npz")
        w = np.array([u["n"] for u in units], float)
        for rd, D in (("14", z["d14"]), ("core", z["dcore"])):
            rws = np.arange(BLOCK[rd], BLOCK[rd] + NB)
            obs = pooled(units, rws)
            if args.smoke:
                pass
            elif rd == "14":
                assert np.allclose(obs, z["obs14"], equal_nan=True)
            else:
                assert np.allclose(obs[1:], z["obsc_comp"], equal_nan=True)
            nul = M.pool_null(w, [u["null"][rd] for u in units])
            r = dict(kind=kind, readout=rd, cohorts=len(units), patients=int(w.sum()), clock=float(obs[0]),
                     p_pairing=float((1 + np.sum(nul >= obs[0])) / (1 + len(nul))), null_mean=float(nul.mean()), null_sd=float(nul.std()))
            for f, sl in FAM_SL.items():
                q95 = np.nanquantile(D[:, sl], .95, axis=1)
                r.update({f"{f}_median": float(np.nanmedian(obs[sl])), f"{f}_q95": float(np.nanquantile(obs[sl], .95)),
                          f"{f}_pi": float(percentile(obs, sl)[0]),
                          f"{f}_pcal": float(calibrated_p(percentile(obs, sl)[0], percentile(D, sl))[0]),
                          f"{f}_placebo_rate": float(np.mean(D[:, 0] > q95)),
                          f"{f}_V346": bool(np.mean(D[:, 0] > q95) <= BAND)})
            rows.append(r)
        # trio contribution: patient bootstrap within cohort, shared resamples for both readouts
        rng = np.random.default_rng([346, 9, KCODE[kind]])
        bt = []
        for _ in range(N_BOOT):
            v = []
            for rd in ("14", "core"):
                st, ww = [], []
                for u in units:
                    idx = rng.integers(0, u["n"], u["n"]) if rd == "14" else u["_idx"]
                    u["_idx"] = idx
                    s = M.raw_stat(u["za"][BLOCK[rd]][idx], u["zb"][BLOCK[rd]][idx])
                    if np.isfinite(s):
                        st.append(s)
                        ww.append(u["n"])
                v.append(float(np.average(st, weights=ww)))
            bt.append(v)
        bt = np.array(bt)
        d = bt[:, 0] - bt[:, 1]
        trio.append(dict(kind=kind, diff_14_minus_core=float(rows[-2]["clock"] - rows[-1]["clock"]),
                         lo90=float(np.quantile(d, .05)), hi90=float(np.quantile(d, .95)), n_boot=N_BOOT))
    T = pd.DataFrame(rows)
    for col in ("p_pairing", "matched_pcal", "responsive_pcal", "responsive_broad_pcal"):
        T[f"{col}_holm"] = np.nan
        for rd in ("14", "core"):
            sel = T.index[(T.readout == rd) & T.kind.isin(PRIMARY)]
            T.loc[sel, f"{col}_holm"] = M.holm(T.loc[sel, col].to_numpy())
    prim = T.kind.isin(PRIMARY)
    T["H1_rule"] = prim & (T.p_pairing_holm < 0.05) & (T.matched_pcal_holm < 0.05) & (T.responsive_pcal_holm < 0.05)
    T["H2_rule"] = prim & (T.responsive_broad_pcal_holm < 0.05)
    T["hypothesis"] = np.where(~prim, "D346", np.where(T.readout == "core", "H346.1 / H346.2", "H346.3"))
    T.to_csv(out / "tests.csv", index=False)
    pd.DataFrame(trio).to_csv(out / "trio.csv", index=False)
    pd.set_option("display.width", 250)
    log("\n" + T.round(4).to_string())
    log("\n" + pd.DataFrame(trio).round(4).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    global N_PERM, N_DRAW, N_BOOT, N_SET
    if args.smoke:
        N_PERM, N_DRAW, N_BOOT = 200, 50, 5
    t0 = time.time()
    (stage_power if args.stage == "power" else stage_score)(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
