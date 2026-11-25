"""A389 score stage: phase order of stromal-vascular classes within organs (normal and tumour samples).

Registered in paper/docs/registered_experiments.md (A389, note 1 and addendum, power notes). Units: the tumour /
adjacent-normal cohorts of the power stage and the healthy-atlas sites of note 1 (`a389.ATLASES`). Companion moments are
read from the power-stage caches (`output/cache/units.pkl`, `units_atlas.pkl`), so clock and companions use the same
units and sets.

Per unit x sample type x pair (A minus B; peak time, + = A later; clocksync's angle is negated):
  clock states on 14 genes and on the core loop (`a149_lib.pair_states`), cross-moment and energies; 2,000 half-samples
  without replacement (shared by both readouts; seed [389, 3, unit]); leave-one-gene-out over the 14 genes and the 364
  three-gene drops (primary pair only); C3 gain |C| / V_B and amplitude ratio sqrt(V_A / V_B).
Per organ: moments summed over units, offset and R, 90 % interval and SE from the half-samples.
Pooled over organs (organs with >= 8 donors): inverse-variance random effects (DerSimonian-Laird) on the offsets in hours,
90 % interval, Q and its p, tau with a Q-profile 90 % interval; the matched sets pooled with the clock's weights (central
90 % band). H389.1 (normal), H389.2 (consistency), H389.3 (tumour), per the registration; decidability from the power note
(descriptive if not decidable). C1: per-gene class offsets from `Dataset.gene_loadings` of both classes on the same host
(the unit's other compartments, normal samples, primary pair) pooled over units; leave-one-gene-out pooled offsets.
C2: the lung atlas subtypes (A272 pairs; listed). D389: 14 minus core on the shared half-samples and among the drops.

    python a389_score.py [--jobs 8] [--smoke]     (paper/common/run.sh, RUN_CAP_GB 30)
"""
from __future__ import annotations

import argparse
import itertools
import pickle
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

import a389 as Q
from a389 import BASE, CORE, HOURS, L, M, OUT, PAIRS, SAMPLES, log
from clocksync.config import limit_threads

N_HALF, MIN_N, BOUND_DECIDE = 2000, 8, 0.5
DROPS = list(itertools.combinations(BASE.genes, 3))
TRIO = tuple(g for g in BASE.genes if g in ("PER1", "PER2", "NR1D1"))
PRIMARY = "endo-fib"


def sums(za, zb):
    za, zb = za - za.mean(), zb - zb.mean()
    return complex(np.sum(za * np.conj(zb))), float(np.sum(np.abs(za) ** 2)), float(np.sum(np.abs(zb) ** 2))


def peak_h(num):
    return -np.angle(num) * HOURS


def load_atlas_tables(atlas, site, classes):
    import importlib.util
    import sys
    spec = importlib.util.spec_from_file_location("compute_a150", Q.EXP / "A150_atlas_specificity/scripts/compute_a150.py")
    A150 = importlib.util.module_from_spec(spec)
    sys.modules["compute_a150"] = A150
    other = sys.modules.pop("_shared", None)
    sys.path.insert(0, str(Q.PAPER / "results/figure_2/scripts/atlas_fits"))
    try:
        spec.loader.exec_module(A150)
    finally:
        sys.path.pop(0)
        if other is not None:
            sys.modules["_shared"] = other
    A150.ATLASES.update(Q.GENOME)
    pb = A150.load_atlas(atlas, None, 0)["pb"]
    genes = list(BASE.genes)
    parts = []
    for cls, comps in classes.items():
        d = pb[pb.compartment.isin(comps)]
        g = d.groupby("donor")
        t = g[genes + ["depth", "n_cells"]].sum(min_count=1).reset_index()
        t["stratum"] = g.stratum.first().to_numpy()
        t["compartment"] = "adjacent_normal:" + cls
        parts.append(t)
    used = {c for v in classes.values() for c in v}
    site_pre = next(iter(classes.values()))[0].split(":")[0] + ":"
    rest = pb[~pb.compartment.isin(used) & pb.compartment.str.startswith(site_pre)][["donor", "compartment", "stratum", "depth", "n_cells"] + genes]
    full = pd.concat(parts + [rest], ignore_index=True)
    return pd.concat(parts, ignore_index=True), full


def score_unit(kind, spec, smoke):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    if kind == "cohort":
        study, cancer = spec
        key, organ = study, Q.ORGAN.get(cancer, cancer)
        try:
            pb = M.table(study)
        except Exception as e:  # noqa: BLE001
            return dict(key=key, error=str(e)[:200])
        samples = {st: (pb, pre) for st, pre in SAMPLES.items()}
        full = {st: pb for st in SAMPLES}
    else:
        atlas, site, organ, classes = spec
        key = f"{atlas}|{site}"
        try:
            tab, fl = load_atlas_tables(atlas, site, classes)
        except Exception as e:  # noqa: BLE001
            return dict(key=key, error=str(e)[:200])
        samples = {"normal": (tab, "adjacent_normal:")}
        full = {"normal": fl}
    out = dict(key=key, organ=organ, kind=kind, res={})
    rng = np.random.default_rng([389, 3, sum(map(ord, key))])
    for st, (pb, pre) in samples.items():
        for pn, (ka, kb) in PAIRS.items():
            a, b = pre + ka, pre + kb
            if a not in set(pb.compartment) or b not in set(pb.compartment):
                continue
            u = L.pair_states(pb, a, b, wf=BASE)
            if u is None:
                continue
            pts, n = u["pts"], u["n"]
            perm = np.random.default_rng([0, sum(map(ord, key))]).permutation(n) if smoke else np.arange(n)
            st_ = {}
            for rd, w in (("14", BASE), ("core", CORE)):
                p = u if rd == "14" else L.pair_states(pb, a, b, wf=w, pts=pts)
                st_[rd] = (p["za"], p["zb"][perm], p["vA"], p["vB"][perm])
            r = dict(n=n, **{f"sums_{rd}": sums(st_[rd][0], st_[rd][1]) for rd in st_})
            halves = [rng.choice(n, max(n // 2, 2), replace=False) for _ in range(N_HALF)]
            r["half"] = {rd: np.array([sums(st_[rd][0][i], st_[rd][1][i]) for i in halves], dtype=object) for rd in st_}
            za, zb = st_["14"][0] - st_["14"][0].mean(), st_["14"][1] - st_["14"][1].mean()
            VA, VB = np.mean(np.abs(za) ** 2) - np.mean(st_["14"][2]), np.mean(np.abs(zb) ** 2) - np.mean(st_["14"][3])
            r["c3"] = dict(gain=float(np.abs(np.mean(za * np.conj(zb))) / VB) if VB > 0 else np.nan,
                           amp_ratio=float(np.sqrt(VA / VB)) if VA > 0 and VB > 0 else np.nan)
            r["logo"] = {}
            for g in BASE.genes:
                p = L.pair_states(pb, a, b, wf=BASE.drop([g]), pts=pts)
                r["logo"][g] = sums(p["za"], p["zb"][perm]) if p is not None else None
            if pn == PRIMARY:
                r["drops"] = {}
                for d in DROPS:
                    p = L.pair_states(pb, a, b, wf=BASE.drop(list(d)), pts=pts)
                    r["drops"][d] = sums(p["za"], p["zb"][perm]) if p is not None else None
                if st == "normal" and not smoke:
                    r["c1"] = gene_offsets(full[st], a, b)
            out["res"][(st, pn)] = r
    return out


def gene_offsets(pb, a, b):
    """C1: per-gene class offset (peak A minus B, h) from both classes' gene loadings on the same host (the unit's
    other compartments, >= 2 needed); None when not available."""
    import clocksync as cs
    hosts = [c for c in pb.compartment.unique() if c not in (a, b)]
    if len(hosts) < 2:
        return None
    try:
        ds = cs.Dataset(pb, anchors=hosts, mixed_strata=True)
        ga, sa = ds.gene_loadings(a, n_perm=0, hosts=hosts)
        gb, sb = ds.gene_loadings(b, n_perm=0, hosts=hosts)
    except Exception:  # noqa: BLE001
        return None
    # gene peak in host time = -z_offset_h (clocksync sign: > 0 advanced) + dev_h (> 0 later) + the gene's reference phase
    d = (-sa.z_offset_h + ga.dev_h) - (-sb.z_offset_h + gb.dev_h)
    return {g: float((v + 12) % 24 - 12) for g, v in d.items()}


# ---- pooling ---------------------------------------------------------------------------------------------------------

def organ_table(units, st, pn, rd):
    rows = {}
    for u in units:
        r = u.get("res", {}).get((st, pn))
        if r is None:
            continue
        o = rows.setdefault(u["organ"], dict(n=0, num=0j, ea=0.0, eb=0.0, half=None, units=[]))
        c, x, y = r[f"sums_{rd}"]
        o["n"] += r["n"]
        o["num"] += c
        o["ea"] += x
        o["eb"] += y
        h = np.array([complex(v[0]) for v in r["half"][rd]])
        o["half"] = h if o["half"] is None else o["half"] + h
        o["units"].append(u["key"])
    out = []
    for organ, o in rows.items():
        if o["n"] < MIN_N:
            continue
        off = peak_h(o["num"])
        dev = (peak_h(o["half"]) - off + 12) % 24 - 12
        out.append(dict(sample=st, pair=pn, readout=rd, organ=organ, units=len(o["units"]), n=o["n"],
                        R=abs(o["num"]) / np.sqrt(o["ea"] * o["eb"]), peak_h=off,
                        lo90=off + np.quantile(dev, .05), hi90=off + np.quantile(dev, .95), se_h=float(np.std(dev))))
    return pd.DataFrame(out)


def random_effects(y, se):
    w = 1 / se ** 2
    mu_f = np.sum(w * y) / np.sum(w)
    q = float(np.sum(w * (y - mu_f) ** 2))
    k = len(y)
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))) if k > 1 else 0.0
    ws = 1 / (se ** 2 + tau2)
    mu = np.sum(ws * y) / np.sum(ws)
    se_mu = 1 / np.sqrt(np.sum(ws))

    def qgen(t2):
        wt = 1 / (se ** 2 + t2)
        m = np.sum(wt * y) / np.sum(wt)
        return np.sum(wt * (y - m) ** 2)

    def bound(target):
        if k < 2 or qgen(0.0) < target:
            return 0.0
        lo, hi = 0.0, 100.0
        for _ in range(80):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if qgen(mid) > target else (lo, mid)
        return (lo + hi) / 2
    tau_hi = np.sqrt(bound(stats.chi2.ppf(0.05, k - 1))) if k > 1 else np.nan
    tau_lo = np.sqrt(bound(stats.chi2.ppf(0.95, k - 1))) if k > 1 else np.nan
    return dict(mu=float(mu), lo90=float(mu - 1.645 * se_mu), hi90=float(mu + 1.645 * se_mu), se=float(se_mu),
                Q=q, Q_p=float(1 - stats.chi2.cdf(q, k - 1)) if k > 1 else np.nan, tau=float(np.sqrt(tau2)),
                tau_lo90=float(tau_lo), tau_hi90=float(tau_hi), k=k, weights=ws / ws.sum())


def companion_band(cache_units, organs, weights, st, pn, rd):
    """Matched sets pooled with the clock's weights: per set, organ offsets from the summed moments."""
    offs = []
    for organ in organs:
        num = None
        for r in cache_units:
            if "comp" not in r or (st, pn) not in r["comp"] or Q.ORGAN.get(r["cancer"], r["cancer"]) != organ:
                continue
            c, x, _ = r["comp"][(st, pn)][rd]
            c = np.where(np.isfinite(x), c, 0)
            num = c if num is None else num + c
        offs.append(peak_h(num) if num is not None else np.full(len(weights) and 200, np.nan))
    offs = np.array(offs)
    pooled = np.sum(weights[:, None] * ((offs - offs[[0]] + 12) % 24 - 12), 0) + offs[0]
    pooled = (pooled + 12) % 24 - 12
    return pooled


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    global N_HALF, DROPS
    if args.smoke:
        N_HALF, DROPS = 50, DROPS[:5] + [TRIO]
    t0 = time.time()
    out = OUT / ("score_smoke" if args.smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    cache = OUT / "cache"
    with open(cache / "units.pkl", "rb") as f:
        cu = pickle.load(f)
    with open(cache / "units_atlas.pkl", "rb") as f:
        au = pickle.load(f)
    cohorts = [("cohort", (r["study"], r["cancer"])) for r in cu if "rows" in r and r["rows"]]
    atl = [("atlas", x) for x in Q.ATLASES if any(a.get("study") == f"{x[0]}|{x[1]}" and a.get("rows") for a in au)]
    jobs = (cohorts[:3] + atl[:1]) if args.smoke else cohorts + atl
    log(f"scoring {len(jobs)} units")
    units = Parallel(n_jobs=args.jobs)(delayed(score_unit)(k, s, args.smoke) for k, s in jobs)
    errs = [u for u in units if "error" in u]
    units = [u for u in units if "error" not in u]
    with open(out / "units_scored.pkl", "wb") as f:
        pickle.dump(units, f)
    log(f"units scored ({len(errs)} errors)")
    cache_units = cu + au
    org_rows, pooled_rows, comp_rows = [], [], []
    for st in SAMPLES:
        for pn in PAIRS:
            for rd in ("14", "core"):
                T = organ_table(units, st, pn, rd)
                if T.empty:
                    continue
                org_rows.append(T)
                re = random_effects(T.peak_h.to_numpy(), T.se_h.to_numpy())
                band = companion_band(cache_units, list(T.organ), re["weights"], st, pn, rd)
                q05, q95 = np.nanquantile(band, [.05, .95])
                same_sign = int(np.sum(np.sign(T.peak_h) == np.sign(re["mu"])))
                includes = int(np.sum((T.lo90 <= re["mu"]) & (T.hi90 >= re["mu"]) | (np.sign(T.peak_h) == np.sign(re["mu"]))))
                pooled_rows.append(dict(sample=st, pair=pn, readout=rd, organs=re["k"], donors=int(T.n.sum()),
                                        pooled_peak_h=re["mu"], lo90=re["lo90"], hi90=re["hi90"], Q=re["Q"], Q_p=re["Q_p"],
                                        tau_h=re["tau"], tau_lo90=re["tau_lo90"], tau_hi90=re["tau_hi90"],
                                        comp_q05=q05, comp_q95=q95, outside_comp=bool(re["mu"] < q05 or re["mu"] > q95),
                                        organs_same_sign=same_sign, organs_compatible=includes))
                comp_rows.append(pd.DataFrame(dict(sample=st, pair=pn, readout=rd, pooled_comp_peak_h=band)))
    O = pd.concat(org_rows)
    P = pd.DataFrame(pooled_rows)
    # verdicts (registered); the decidability flags come from the power notes
    P["excludes0_mouse_sign"] = (P.hi90 < 0)
    P["H389_1_rule"] = (P["sample"] == "normal") & P.excludes0_mouse_sign & P.outside_comp
    P["H389_2_rule"] = P.H389_1_rule & (P.tau_hi90 < np.abs(P.pooled_peak_h) / 2) & (P.organs_compatible == P.organs)
    P["H389_3_rule"] = (P["sample"] == "tumour") & P.excludes0_mouse_sign & P.outside_comp
    O.to_csv(out / "organs.csv", index=False)
    P.to_csv(out / "pooled.csv", index=False)
    pd.concat(comp_rows).to_csv(out / "companion_pooled.csv", index=False)
    pd.DataFrame(errs).to_csv(out / "errors.csv", index=False)
    # C1: leave-one-gene-out and per-gene offsets (primary pair)
    c1 = []
    for st in SAMPLES:
        base = P[(P["sample"] == st) & (P.pair == PRIMARY) & (P.readout == "14")]
        if base.empty:
            continue
        mu0 = float(base.pooled_peak_h.iloc[0])
        for g in BASE.genes:
            rows = {}
            for u in units:
                r = u["res"].get((st, PRIMARY))
                if r is None or r["logo"].get(g) is None:
                    continue
                o = rows.setdefault(u["organ"], [0j, 0])
                o[0] += r["logo"][g][0]
                o[1] += r["n"]
            T = pd.DataFrame([dict(organ=k, peak=peak_h(v[0]), n=v[1]) for k, v in rows.items() if v[1] >= MIN_N])
            mu = float(np.angle(np.sum(T.n * np.exp(1j * T.peak / HOURS))) * HOURS) if len(T) else np.nan
            c1.append(dict(sample=st, gene=g, logo_pooled_peak_h=mu, shift_vs_all=mu - mu0,
                           moves_more_than_half=bool(abs(mu - mu0) > abs(mu0) / 2)))
        if st == "normal":
            per_gene = {}
            for u in units:
                r = u["res"].get((st, PRIMARY))
                if r is None or not r.get("c1"):
                    continue
                for g, v in r["c1"].items():
                    per_gene.setdefault(g, []).append((v, r["n"]))
            for g, v in per_gene.items():
                x = np.array([a for a, _ in v]); w = np.array([b for _, b in v], float)
                ok = np.isfinite(x)
                c1.append(dict(sample=st, gene=g, gene_offset_peak_h=float(np.angle(np.sum(w[ok] * np.exp(1j * x[ok] / HOURS))) * HOURS)
                               if ok.any() else np.nan, gene_units=int(ok.sum())))
    pd.DataFrame(c1).to_csv(out / "c1_genes.csv", index=False)
    # C3 and D389
    c3, d389 = [], []
    for u in units:
        for (st, pn), r in u["res"].items():
            c3.append(dict(unit=u["key"], organ=u["organ"], sample=st, pair=pn, n=r["n"], **r["c3"]))
    for st in SAMPLES:
        T14, Tc = organ_table(units, st, PRIMARY, "14"), organ_table(units, st, PRIMARY, "core")
        if T14.empty:
            continue
        w = random_effects(T14.peak_h.to_numpy(), T14.se_h.to_numpy())["weights"]
        organs = list(T14.organ)
        hd = {rd: [] for rd in ("14", "core")}
        for rd in hd:
            for organ in organs:
                h = None
                for u in units:
                    r = u["res"].get((st, PRIMARY))
                    if r is None or u["organ"] != organ:
                        continue
                    x = np.array([complex(v[0]) for v in r["half"][rd]])
                    h = x if h is None else h + x
                hd[rd].append(peak_h(h))
        p14 = np.sum(w[:, None] * np.array(hd["14"]), 0)
        pc = np.sum(w[:, None] * np.array(hd["core"]), 0)
        diff = p14 - pc
        obs14 = float(np.sum(w * T14.set_index("organ").loc[organs].peak_h))
        obsc = float(np.sum(w * Tc.set_index("organ").loc[organs].peak_h))
        dv = []
        for d in DROPS:
            num = {}
            for u in units:
                r = u["res"].get((st, PRIMARY))
                if r is None or u["organ"] not in organs or r["drops"].get(d) is None:
                    continue
                num[u["organ"]] = num.get(u["organ"], 0j) + r["drops"][d][0]
            dv.append(float(np.sum(w * np.array([peak_h(num.get(o, np.nan)) for o in organs]))))
        dv = np.array(dv)
        tr = dv[DROPS.index(TRIO)]
        d389.append(dict(sample=st, pooled_14=obs14, pooled_core=obsc, diff=obs14 - obsc,
                         diff_lo90=obs14 - obsc + np.quantile(diff - diff.mean(), .05),
                         diff_hi90=obs14 - obsc + np.quantile(diff - diff.mean(), .95),
                         trio_drop_percentile=float(np.mean(np.abs(dv - obs14) <= abs(tr - obs14))),
                         drops_median=float(np.nanmedian(dv))))
    pd.DataFrame(c3).to_csv(out / "c3_amplitude.csv", index=False)
    pd.DataFrame(d389).to_csv(out / "d389_gene_drop.csv", index=False)
    pd.read_csv(Q.A272_PAIRS).pipe(lambda a: a[(a.a.str.contains("EC_") | a.b.str.contains("EC_"))
                                               & (a.family_a.str.contains("fibro") | a.family_b.str.contains("fibro"))]) \
        .to_csv(out / "c2_lung_subtypes.csv", index=False)
    pd.set_option("display.width", 250)
    log("organs\n" + O[O.pair == PRIMARY].round(3).to_string(index=False))
    log("pooled\n" + P.round(3).to_string(index=False))
    log("D389\n" + pd.DataFrame(d389).round(3).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
