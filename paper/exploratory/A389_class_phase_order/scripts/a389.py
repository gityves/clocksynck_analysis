"""A389: phase order of stromal-vascular cell classes within organs (normal and tumour samples).

Registered in paper/docs/registered_experiments.md (A389). Cohorts: the A388 feasibility inventory (A149 / A226 / A227
lists; the Wu and Pal breast atlases are not included here: their tables need the author-label builder), tables
`compute_a227.table` (A149 rules). Sample types: `adjacent_normal:` and `tumour:`. Pairs: endothelial - fibroblast
(primary), endothelial - mural, mural - fibroblast. Organ: the cohort's tissue of origin (normal) or cancer (tumour).

  power  (no clock offset computed: the clock states of the pairs are never cross-multiplied)
         donors with both classes per cohort, sample type, pair and organ; organs with >= 8 donors; analytic 90 %
         half-width of each organ's offset at R 0.5 and 0.3 and the pooled MDE (2.49 x pooled SE at R 0.5, inverse-
         variance pooling over organs); companions: 200 expression-matched sets per cohort x sample type x pair (seeds
         [389, sample type, pair, cohort]) on genome-wide pseudobulks (`compute_a227.genome_pair`), cross-moments
         pooled per organ as the clock will be, 14 genes and core positions; the lung-atlas subtype pairs available for C2
         (A272 pair table: counts only).                                  -> output/power/*.csv, output/cache/
  score  (written after the power note)

    python a389.py --stage power [--jobs 4] [--genome-jobs 3] [--smoke]     (paper/common/run.sh, RUN_CAP_GB 30)
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
sys.path.insert(0, str(EXP / "A387_tumour_full_model/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a387 as R  # noqa: E402  (sets the A149 / A227 / A336 paths)
from a387 import C, HOURS, M, log  # noqa: E402
import a149_lib as L  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
INVENTORY = EXP / "A388_malignant_tumour_host/output/feasibility/inventory.csv"
A272_PAIRS = EXP / "A272_backbone_pairwise/output/pairs.csv"
BASE = C.WF
CORE = C.core(BASE)
SAMPLES = {"normal": "adjacent_normal:", "tumour": "tumour:"}
PAIRS = {"endo-fib": ("endothelial", "fibroblast"), "endo-mural": ("endothelial", "mural"),
         "mural-fib": ("mural", "fibroblast")}
ORGAN = {"colorectal": "colon", "gastric": "stomach", "esophageal": "oesophagus", "pancreatic": "pancreas",
         "head and neck": "head and neck", "liver": "liver", "lung": "lung", "prostate": "prostate", "thyroid": "thyroid",
         "breast": "breast", "kidney": "kidney", "bladder": "bladder", "ovarian": "ovary", "cervical": "cervix",
         "skin": "skin", "ampullary": "ampulla"}
N_SET, MIN_N, R_GRID = 200, 8, (0.5, 0.3)


def studies(smoke):
    inv = pd.read_csv(INVENTORY)
    inv = inv[~inv.source.str.startswith("breast atlas") & inv.error.isna()]
    s = inv[["study", "cancer"]].drop_duplicates()
    return s.head(4) if smoke else s


def _pairs(out, pb, genome, key, si_list):
    """Counts and companion cross-moments of every pair x sample type of one table (no clock cross-moment).
    genome(a, b, pts) -> genome-wide table of the two compartments."""
    for si, (st, pre) in si_list:
        for pi, (pn, (ka, kb)) in enumerate(PAIRS.items()):
            a, b = pre + ka, pre + kb
            if a not in set(pb.compartment) or b not in set(pb.compartment):
                continue
            u = L.pair_states(pb, a, b, wf=BASE)
            if u is None:
                continue
            out["rows"].append(dict(study=out["study"], cancer=out["cancer"], sample=st, pair=pn, n=u["n"]))
            try:
                G = genome(a, b, u["pts"])
            except Exception as e:  # noqa: BLE001
                out["rows"][-1]["genome_error"] = f"{type(e).__name__}: {e}"[:200]
                continue
            ref = G[G.compartment == b]
            genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
            frac = ref[genes].div(ref.depth, axis=0)
            lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
            rng = np.random.default_rng([389, si, pi, sum(map(ord, key))])
            sets = pseudo_waveforms(lm, N_SET, rng, BASE, "matched")
            res = {}
            for rd, wfs in (("14", sets), ("core", [C.restrict(w, BASE, CORE) for w in sets])):
                num, ea, eb = np.full(N_SET, np.nan + 0j), np.full(N_SET, np.nan), np.full(N_SET, np.nan)
                for i, w in enumerate(wfs):
                    p = L.pair_states(G, a, b, wf=w, pts=u["pts"])
                    if p is None:
                        continue
                    c, x, y = R._sums(p["za"][None], p["zb"][None])
                    num[i], ea[i], eb[i] = c[0], x[0], y[0]
                res[rd] = (num, ea, eb)
            out["comp"][(st, pn)] = res
    return out


def unit(study, cancer):
    """A tumour / adjacent-normal cohort: both sample types."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        pb = M.table(study)
    except Exception as e:  # noqa: BLE001
        return dict(study=study, error=f"{type(e).__name__}: {e}"[:200])
    out = dict(study=study, cancer=cancer, rows=[], comp={}, kind="cohort")
    return _pairs(out, pb, lambda a, b, pts: M.genome_pair(study, (a, b), pts), study, list(enumerate(SAMPLES.items())))


# healthy atlases (A389 note 1): atlas, site, organ, class -> compartments (summed per donor when several)
ATLASES = [
    ("heart", "heart", "heart", {"endothelial": ["heart:endothelial"], "fibroblast": ["heart:fibroblast"], "mural": ["heart:mural"]}),
    ("lung_hlca", "lung", "lung", {"endothelial": ["lung:EC_arterial", "lung:EC_capillary", "lung:EC_venous"],
                                   "fibroblast": ["lung:fibroblast"], "mural": ["lung:smooth_muscle"]}),
    ("organ_breast", "breast", "breast", {"endothelial": ["breast:endothelial"], "fibroblast": ["breast:fibroblast_stromal"],
                                          "mural": ["breast:mural_smooth_muscle"]}),
    ("organ_colon", "colon", "colon", {"endothelial": ["colon:endothelial"], "fibroblast": ["colon:fibroblast_stromal"]}),
    ("kong_colon", "colon", "colon", {"endothelial": ["colon:endothelial"], "fibroblast": ["colon:fibroblast_stromal"]}),
    ("kong_ileum", "ileum", "ileum", {"endothelial": ["ileum:endothelial"], "fibroblast": ["ileum:fibroblast_stromal"]}),
    ("kidney_sc", "kidney", "kidney", {"endothelial": ["kidney:endothelial"], "fibroblast": ["kidney:interstitial"]}),
    ("organ_eye_anterior", "ciliary", "anterior eye", {"endothelial": ["eye_anterior:ciliary_endothelial"],
                                                       "fibroblast": ["eye_anterior:ciliary_fibroblast_stromal"],
                                                       "mural": ["eye_anterior:ciliary_mural_smooth_muscle"]}),
    ("organ_eye_anterior", "trabecular", "anterior eye", {"endothelial": ["eye_anterior:trabecular_endothelial"],
                                                          "fibroblast": ["eye_anterior:trabecular_fibroblast_stromal"]}),
    ("organ_optic_nerve", "head", "optic nerve", {"endothelial": ["optic_nerve:head_endothelial"],
                                                  "fibroblast": ["optic_nerve:head_fibroblast_stromal"],
                                                  "mural": ["optic_nerve:head_mural_smooth_muscle"]}),
    ("organ_optic_nerve", "nerve", "optic nerve", {"endothelial": ["optic_nerve:nerve_endothelial"],
                                                   "fibroblast": ["optic_nerve:nerve_fibroblast_stromal"],
                                                   "mural": ["optic_nerve:nerve_mural_smooth_muscle"]}),
    ("endometrium", "endometrium", "endometrium", {"endothelial": ["endometrium:endothelial"], "fibroblast": ["endometrium:stromal"]}),
    ("skin_ad_control", "skin", "skin", {"endothelial": ["skin:endothelial"], "fibroblast": ["skin:fibroblast"], "mural": ["skin:pericyte"]}),
    ("lung_natri_control", "lung", "lung", {"endothelial": ["lung:EC_venous"], "fibroblast": ["lung:fibroblast"]}),
    ("islets", "islet", "pancreas", {"endothelial": ["islet:endothelial"], "fibroblast": ["islet:stellate"]}),
]
BLIND_ATLASES = {"endometrium", "skin_ad_control", "lung_natri_control", "islets"}


def atlas_unit(atlas, site, organ, classes):
    """A healthy atlas site: genome-wide pseudobulks (`compute_a150.load_atlas`), classes relabelled
    `adjacent_normal:<class>` (summed per donor over several compartments) so the cohort code applies."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    import importlib.util
    spec = importlib.util.spec_from_file_location("compute_a150", EXP / "A150_atlas_specificity/scripts/compute_a150.py")
    A150 = importlib.util.module_from_spec(spec)
    sys.modules["compute_a150"] = A150
    other = sys.modules.pop("_shared", None)          # the cancer `_shared` (figure 5) shadows the atlas_fits one
    sys.path.insert(0, str(PAPER / "results/figure_2/scripts/atlas_fits"))
    try:
        spec.loader.exec_module(A150)
    finally:
        sys.path.pop(0)
        if other is not None:
            sys.modules["_shared"] = other
    key = f"{atlas}|{site}"
    try:
        pb = A150.load_atlas(atlas, None, 0)["pb"]
    except Exception as e:  # noqa: BLE001
        return dict(study=key, error=f"{type(e).__name__}: {e}"[:200])
    parts = []
    for cls, comps in classes.items():
        d = pb[pb.compartment.isin(comps)]
        if not len(d):
            continue
        num = [c for c in d.columns if c not in ("donor", "compartment", "stratum")]
        g = d.groupby("donor")
        t = g[num].sum(min_count=1).reset_index()
        t["stratum"] = g.stratum.first().to_numpy()
        t["compartment"] = "adjacent_normal:" + cls
        parts.append(t)
    tab = pd.concat(parts, ignore_index=True)
    out = dict(study=key, cancer=organ, rows=[], comp={}, kind="atlas", blind=atlas in BLIND_ATLASES)
    return _pairs(out, tab, lambda a, b, pts: tab[tab.compartment.isin([a, b]) & tab.donor.isin(pts)], key,
                  [(2, ("normal", "adjacent_normal:"))])


def halfwidth_h(rho, n):
    return 1.645 * np.sqrt((1 - rho ** 2) / (2 * n * rho ** 2)) * HOURS


def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    S = studies(args.smoke)
    log(f"{len(S)} cohorts")
    res = Parallel(n_jobs=args.genome_jobs)(delayed(unit)(s, c) for s, c in zip(S.study, S.cancer))
    with open(cache / "units.pkl", "wb") as f:
        pickle.dump(res, f)
    summarize(res, out)


def stage_power_atlas(args):
    """A389 note 1: the cohort units of the power stage plus the healthy-atlas units."""
    out = OUT / ("power_extended_smoke" if args.smoke else "power_extended")
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out.mkdir(parents=True, exist_ok=True)
    with open(cache / "units.pkl", "rb") as f:
        res = pickle.load(f)
    specs = ATLASES[:2] if args.smoke else ATLASES
    log(f"{len(specs)} atlas sites")
    at = Parallel(n_jobs=args.genome_jobs)(delayed(atlas_unit)(*x) for x in specs)
    with open(cache / "units_atlas.pkl", "wb") as f:
        pickle.dump(at, f)
    summarize(res + at, out)


def summarize(res, out):
    errs = pd.DataFrame([r for r in res if "error" in r])
    rows = pd.DataFrame([x for r in res if "rows" in r for x in r["rows"]])
    rows["organ"] = rows.cancer.map(ORGAN).fillna(rows.cancer)
    rows.to_csv(out / "counts_cohort.csv", index=False)
    errs.to_csv(out / "errors.csv", index=False)
    org = rows.groupby(["sample", "pair", "organ"]).agg(cohorts=("study", "nunique"), n=("n", "sum")).reset_index()
    org["included"] = org.n >= MIN_N
    for rho in R_GRID:
        org[f"halfwidth90_h_R{rho}"] = halfwidth_h(rho, org.n)
    org.to_csv(out / "counts_organ.csv", index=False)
    # pooled MDE and companion offsets per organ and pooled
    pw, comp = [], []
    for (st, pn), g in org[org.included].groupby(["sample", "pair"]):
        se = np.sqrt((1 - 0.25) / (2 * g.n * 0.25))                       # rad, R 0.5
        pse = 1 / np.sqrt(np.sum(1 / se ** 2))
        pw.append(dict(sample=st, pair=pn, organs=len(g), donors=int(g.n.sum()), pooled_se_h=float(pse * HOURS),
                       mde_h=float(2.49 * pse * HOURS), decidable=bool(2.49 * pse * HOURS <= 0.5)))
        for rd in ("14", "core"):
            org_off, w = [], []
            for organ in g.organ:
                num, ea, eb = np.zeros(N_SET, complex), np.zeros(N_SET), np.zeros(N_SET)
                for r in res:
                    if "comp" not in r or (st, pn) not in r["comp"] or ORGAN.get(r["cancer"], r["cancer"]) != organ:
                        continue
                    c, x, y = r["comp"][(st, pn)][rd]
                    ok = np.isfinite(x)
                    num, ea, eb = num + np.where(ok, c, 0), ea + np.where(ok, x, 0), eb + np.where(ok, y, 0)
                off = np.angle(num) * HOURS * -1                           # peak time, + = first class later
                Rc = np.abs(num) / np.sqrt(ea * eb)
                comp.append(dict(sample=st, pair=pn, readout=rd, organ=organ, comp_R_median=float(np.nanmedian(Rc)),
                                 comp_R_q95=float(np.nanquantile(Rc, .95)), comp_offset_q05=float(np.nanquantile(off, .05)),
                                 comp_offset_median=float(np.nanmedian(off)), comp_offset_q95=float(np.nanquantile(off, .95)),
                                 comp_offset_circ_sd_h=float(np.sqrt(-2 * np.log(np.abs(np.nanmean(np.exp(1j * off / HOURS)))))
                                                             * HOURS)))
                org_off.append(np.exp(1j * off / HOURS))
                w.append(float(g[g.organ == organ].n.iloc[0]))
            pooled = np.angle(np.average(np.array(org_off), axis=0, weights=w)) * HOURS
            comp.append(dict(sample=st, pair=pn, readout=rd, organ="pooled (n-weighted)",
                             comp_offset_q05=float(np.quantile(pooled, .05)), comp_offset_median=float(np.median(pooled)),
                             comp_offset_q95=float(np.quantile(pooled, .95))))
    pd.DataFrame(pw).to_csv(out / "power.csv", index=False)
    pd.DataFrame(comp).to_csv(out / "companions.csv", index=False)
    a = pd.read_csv(A272_PAIRS)
    sub = a[(a.a.str.contains("EC_") | a.b.str.contains("EC_")) & (a.family_a.str.contains("fibro") | a.family_b.str.contains("fibro"))]
    sub[["atlas", "a", "b", "n"]].to_csv(out / "c2_subtypes.csv", index=False)
    pd.set_option("display.width", 250)
    log("organs\n" + org[org.included].round(3).to_string(index=False))
    log("power\n" + pd.DataFrame(pw).round(3).to_string(index=False))
    log("companions\n" + pd.DataFrame(comp).round(3).to_string(index=False))
    log("C2 subtypes\n" + sub[["atlas", "a", "b", "n"]].to_string(index=False))



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "power_atlas"), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    global N_SET
    if args.smoke:
        N_SET = 10
    t0 = time.time()
    (stage_power if args.stage == "power" else stage_power_atlas)(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
