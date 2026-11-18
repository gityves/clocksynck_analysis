"""A313 (registered): do anatomically adjacent GTEx tissues share a person's departure from the body's time more than
distant tissues, at equal composition? Aggregates only (GTEx controlled access).

A104 periphery rows without Whole Blood and Testis (25 tissues). Per pair: A302 estimator (host without both members,
strata intercepts, host-error-corrected |correlation| of residuals; >= 30 shared donors). Cells from fixed regions and
families. T_adj = median r(AD) - median r(DD); T_fam = median r(DS) - median r(DD). Nulls: region (or family)
reassignment of the tissues (10,000), and 100 matched pseudo-clock sets (seed [313, 1]).

    python score_a313.py --stage power [--jobs 4]   # 20 matched sets only: SDs of T_adj and the power of H313.1
    python score_a313.py --stage score [--jobs 4]   # clock and 100 matched sets; H313.1, H313.2
Outputs (`../output/`): power.csv, power_sets.csv | pairs.csv, sets.csv, tests.csv.
"""
import argparse
import itertools
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
sys.path.insert(0, str(HERE.parents[2] / "A302_departure_groups" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
import score_a302 as A302  # noqa: E402
from clocksync._core import host as HC  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
EXCLUDE = ["Whole Blood", "Testis"]
REGIONS = {"leg": ["Skin - Sun Exposed (Lower leg)", "Adipose - Subcutaneous", "Muscle - Skeletal", "Artery - Tibial",
                   "Nerve - Tibial"],
           "heart": ["Heart - Atrial Appendage", "Heart - Left Ventricle", "Artery - Coronary", "Artery - Aorta"],
           "upper_gi": ["Esophagus - Mucosa", "Esophagus - Muscularis", "Esophagus - Gastroesophageal Junction", "Stomach"],
           "abdomen": ["Colon - Transverse", "Colon - Sigmoid", "Adipose - Visceral (Omentum)"]}
FAMILIES = {"squamous": ["Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)", "Esophagus - Mucosa"],
            "adipose": ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Breast - Mammary Tissue"],
            "artery": ["Artery - Tibial", "Artery - Aorta", "Artery - Coronary"],
            "heart": ["Heart - Atrial Appendage", "Heart - Left Ventricle"],
            "gi_smooth_muscle": ["Esophagus - Muscularis", "Esophagus - Gastroesophageal Junction", "Colon - Sigmoid"],
            "gi_mucosa": ["Stomach", "Colon - Transverse"]}
META = A302.META
N_DRAWS, DELTA = 10_000, 0.10
S = {}


def labels(groups: dict, tissues) -> dict:
    lab = {t: f"_{t}" for t in tissues}
    for g, ts in groups.items():
        for t in ts:
            lab[t] = g
    return lab


def cells(pairs: pd.DataFrame, region: dict, family: dict) -> np.ndarray:
    adj = np.array([region[a] == region[b] for a, b in zip(pairs.a, pairs.b)])
    sim = np.array([family[a] == family[b] for a, b in zip(pairs.a, pairs.b)])
    return np.where(adj & sim, "AS", np.where(adj, "AD", np.where(sim, "DS", "DD")))


def T_of(r, cell, x):
    return float(np.nanmedian(r[cell == x]) - np.nanmedian(r[cell == "DD"]))


def pair_r(d, comps, a, b):
    h = HC.host(d, [c for c in comps if c not in (a, b)])
    ea, aa, ta = A302.residuals(d, a, h)
    eb, ab, _ = A302.residuals(d, b, h)
    sh = ea.index.intersection(eb.index)
    if len(sh) < A302.MIN_SHARED:
        return dict(a=a, b=b, n=len(sh), r=np.nan)
    x, y, t = ea[sh].to_numpy(), eb[sh].to_numpy(), ta[sh].to_numpy()
    m = np.mean(x * y.conj()) - aa * np.conj(ab) * t.mean()
    va = np.mean(np.abs(x) ** 2) - abs(aa) ** 2 * t.mean()
    vb = np.mean(np.abs(y) ** 2) - abs(ab) ** 2 * t.mean()
    return dict(a=a, b=b, n=len(sh), r=float(abs(m) / np.sqrt(va * vb)) if va > 0 and vb > 0 else np.nan)


def all_pairs(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    w = S["sets"][i]
    pb, comps = S["pb"], S["comps"]
    try:
        d = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=f"s{i}").data
        return pd.DataFrame([pair_r(d, comps, a, b) for a, b in itertools.combinations(comps, 2)]).assign(set=i)
    except Exception as e:  # noqa: BLE001
        return pd.DataFrame([dict(set=i, error=f"{type(e).__name__}: {e}"[:200])])


def reassign_null(p, groups, fixed, x, rng, n):
    """T for random reassignments of tissues to the slots of `groups` (sizes kept); `fixed` labels the other axis."""
    tissues = S["comps"]
    out = np.empty(n)
    for k in range(n):
        perm = dict(zip(tissues, rng.permutation(tissues)))
        g = {name: [perm[t] for t in ts] for name, ts in groups.items()}
        lab = labels(g, tissues)
        c = cells(p, lab, fixed) if groups is REGIONS else cells(p, fixed, lab)
        out[k] = T_of(p.r.to_numpy(), c, x)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t0 = time.time()
    import anndata as ad
    h5 = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(h5.var_names)
    h5.file.close()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    limit_threads(1)
    pb = pb[~pb.compartment.str.startswith("Brain") & ~pb.compartment.isin(EXCLUDE)].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    gn = [c for c in pb.columns if c not in META]
    xg = pb[gn].to_numpy(np.float32)
    lm = pd.Series(np.log((xg / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=gn)[xg.mean(0) >= 1]
    del xg
    wf = cs.gtex()
    matched = pseudo_waveforms(lm, 100, np.random.default_rng([313, 1]), wf, kind="matched")
    sets = ([] if a.stage == "power" else [wf]) + (matched[:20] if a.stage == "power" else matched)
    keep = sorted(set(g for w in sets for g in w.genes))
    S.update(pb=pb[META + keep].copy(), sets=sets, comps=comps)
    del pb
    region, family = labels(REGIONS, comps), labels(FAMILIES, comps)
    print(f"{len(comps)} tissues, {len(sets)} gene sets ({time.time() - t0:.0f}s)", flush=True)
    res = []
    with Pool(a.jobs) as pool:
        for r in pool.imap_unordered(all_pairs, range(len(sets))):
            res.append(r)
            print(f"{len(res)}/{len(sets)} sets ({time.time() - t0:.0f}s)", flush=True)
    allp = pd.concat(res, ignore_index=True).sort_values(["set", "a", "b"])
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng([313, 2])
    if a.stage == "power":
        rows = []
        for i, p in allp.groupby("set"):
            c = cells(p, region, family)
            null = reassign_null(p, REGIONS, family, "AD", rng, 1000)
            rows.append(dict(set=i, T_adj=T_of(p.r.to_numpy(), c, "AD"), null_sd=float(np.nanstd(null)),
                             n_AD=int(((c == "AD") & p.r.notna()).sum()), n_DD=int(((c == "DD") & p.r.notna()).sum()),
                             n_DS=int(((c == "DS") & p.r.notna()).sum()), n_AS=int(((c == "AS") & p.r.notna()).sum()),
                             median_r=float(p.r.median())))
        ps = pd.DataFrame(rows)
        ps.to_csv(OUT / "power_sets.csv", index=False)
        sd = max(float(ps.T_adj.std(ddof=1)), float(ps.null_sd.median()))
        res = dict(sets=len(ps), sd_T_over_sets=float(ps.T_adj.std(ddof=1)), median_null_sd=float(ps.null_sd.median()),
                   sd_used=sd, delta=DELTA, power=float(norm.cdf(DELTA / sd - 1.645)), median_r_matched=float(ps.median_r.median()),
                   n_AD=int(ps.n_AD.median()), n_DD=int(ps.n_DD.median()), n_DS=int(ps.n_DS.median()), n_AS=int(ps.n_AS.median()))
        res["gate"] = bool(res["power"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).to_string(), flush=True)
        return
    clock = allp[allp.set == 0].copy()
    clock["cell"] = cells(clock, region, family)
    clock.to_csv(OUT / "pairs.csv", index=False)
    r = clock.r.to_numpy()
    T_adj, T_fam = T_of(r, clock.cell.to_numpy(), "AD"), T_of(r, clock.cell.to_numpy(), "DS")
    n_adj = reassign_null(clock, REGIONS, family, "AD", rng, N_DRAWS)
    n_fam = reassign_null(clock, FAMILIES, region, "DS", rng, N_DRAWS)
    srows = []
    for i, p in allp[allp.set > 0].groupby("set"):
        c = cells(p, region, family)
        srows.append(dict(set=i, T_adj=T_of(p.r.to_numpy(), c, "AD"), T_fam=T_of(p.r.to_numpy(), c, "DS"),
                          **{f"median_{x}": float(np.nanmedian(p.r.to_numpy()[c == x])) for x in ("AS", "AD", "DS", "DD")}))
    st = pd.DataFrame(srows)
    st.to_csv(OUT / "sets.csv", index=False)
    tests = []
    for name, T, null, col in (("H313.1", T_adj, n_adj, "T_adj"), ("H313.2", T_fam, n_fam, "T_fam")):
        p_re = float((1 + np.sum(null >= T)) / (1 + len(null)))
        q95 = float(st[col].quantile(0.95))
        tests.append(dict(test=name, T=T, p_reassign=p_re, matched_q95=q95, matched_median=float(st[col].median()),
                          percentile_vs_matched=float(np.mean(st[col] < T)), passes=bool(T > 0 and p_re < 0.05 and T > q95)))
    tests.append(dict(test="descriptive", **{f"median_r_{x}": float(clock[clock.cell == x].r.median()) for x in ("AS", "AD", "DS", "DD")},
                      **{f"n_{x}": int(clock[(clock.cell == x)].r.notna().sum()) for x in ("AS", "AD", "DS", "DD")}))
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    print(clock[clock.cell != "DD"].round(3).sort_values(["cell", "r"]).to_string(index=False))
    print(pd.DataFrame(tests).round(4).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
