"""D1b feasibility (before any registration): can the Cao 2020 fetal design detect a clock with three readouts?

No fetal clock statistic is computed here. Adult material is cut to the fetal design (number of donors and pseudobulk
depth of the matching Cao organ or cell class), and the share of subsamples where the readout passes is the
sensitivity of that readout at the fetal design (feasible if >= 0.8, the A181 rule).

1. L1 per organ: clock-gene coherence R (`clocksync.decompose.components`, GTEx waveform) of one compartment across
   donors against a per-gene donor permutation null (breaks the gene covariance, keeps each gene's distribution);
   pass = one-sided p < 0.05. Calibrators: GTEx bulk organs, deaths at 08-18 h (fetal material is collected in hospital
   hours), thinned to the Cao median whole-organ depth; and adult nucleus atlases (heart, kidney, lung) per cell class,
   thinned to the Cao median class depth.
2. L2 within organ across cell classes: `Dataset.gate` (1,000 permutations; pass = p_pairing and p_acrophase < 0.05) on
   adult nucleus atlases with the Cao classes of that organ, at the number of fetuses with all those classes.
   The adult sampling times are unknown (no hospital-hours window).
3. Clock maturation markers: counts of the core (ARNTL, CLOCK, NPAS2, CRY1) and output (DBP, NR1D1, NR1D2, PER3, TEF)
   genes in the fetal organ and class pseudobulks (detection only; no ratio, no age trend), and (`--markers`) the
   counts expected at the fetal depth if the class expressed each gene at the adult level.

Controlled access: GTEx time of death is used in memory to select donors; only aggregates are written.

   `--adrenal`: proxy L2 calibration for the fetal adrenal classes (heart and kidney atlases at the adrenal design).

    python feasibility.py [--n-sub 200] [--jobs 4] [--check] [--all-day | --adrenal]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, logfrac  # noqa: E402
from common import data, gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_PERM_L1, N_PERM_GATE = 1810, 200, 1000
MIN_ORGAN_CELLS, MIN_CLASS_CELLS = 1000, 200
WF = cs.gtex()
GENES = list(WF.genes)
GTEX = {"Lung": "Lung", "Kidney": "Kidney - Cortex", "Heart": "Heart - Left Ventricle",
        "Cerebrum": "Brain - Frontal Cortex (BA9)", "Cerebellum": "Brain - Cerebellum", "Liver": "Liver",
        "Muscle": "Muscle - Skeletal", "Adrenal": "Adrenal Gland", "Spleen": "Spleen", "Stomach": "Stomach",
        "Pancreas": "Pancreas", "Intestine": "Small Intestine - Terminal Ileum"}
# adult nucleus atlas, its largest single stratum, and Cao class -> adult compartment
ATLAS = {"Heart": ("heart_dcm_healthy", {"parenchyma": "heart:cardiomyocyte", "stromal": "heart:fibroblast",
                                         "endothelial": "heart:endothelial"}),
         "Lung": ("lung_hlca", {"epithelial": "lung:AT2", "stromal": "lung:fibroblast",
                                "endothelial": "lung:EC_capillary"}),
         "Kidney": ("kidney_sc", {"parenchyma": "kidney:proximal_tubule", "stromal": "kidney:interstitial",
                                  "endothelial": "kidney:endothelial"})}
MIN_FETUSES_L2, POOL_FACTOR = 8, 1.5  # L2 uses the classes present in >= 8 fetuses; adult pool >= 1.5 x fetal n
CORE, OUTPUT = ["ARNTL", "CLOCK", "NPAS2", "CRY1"], ["DBP", "NR1D1", "NR1D2", "PER3", "TEF"]
JOBS = []  # (kind, label, table pool by donor, n donors, {compartment: target depth}, compartments)


def thin(tab, target, rng):
    """Binomial thinning of every gene to the target depth of its compartment (depth scaled by the same p)."""
    t = tab.copy()
    p = np.minimum(1.0, t.compartment.map(target).to_numpy(float) / t.depth.to_numpy(float))
    t[GENES] = rng.binomial(t[GENES].to_numpy(np.int64), p[:, None])
    t["depth"] = t.depth * p
    return t


def l1(tab, rng):
    """Coherence R of one compartment and its one-sided p against per-gene donor permutations."""
    y, pv, k = logfrac(tab, GENES)
    obs = components(y, pv, k, WF)["coherence"]
    null = np.empty(N_PERM_L1)
    for b in range(N_PERM_L1):
        yb = np.column_stack([rng.permutation(c) for c in y.T])
        null[b] = components(yb, pv, k, WF)["coherence"]
    return obs, (1 + np.sum(null >= obs)) / (1 + N_PERM_L1)


def one(args):
    j, i = args
    kind, label, tabs, n, target, comps = JOBS[j]
    rng = np.random.default_rng([SEED, j, i])
    donors = rng.choice(sorted(tabs), size=n, replace=False)
    raw = pd.concat([tabs[d] for d in donors], ignore_index=True)
    rows = []
    for tag, tab in (("thinned", thin(raw, target, rng)), ("unthinned", raw)):
        if kind == "L1":
            r, p = l1(tab.assign(stratum="one"), rng)
            rows.append(dict(kind=kind, label=label, sub=i, depth=tag, n=n, stat=r, p=p, passed=p < .05))
        else:
            g = cs.Dataset(tab, anchors=None, name="d1b").gate(n_perm=N_PERM_GATE, seed=SEED + i)
            rows.append(dict(kind=kind, label=label, sub=i, depth=tag, n=n, stat=g.clock_share, p=g.p_pairing,
                             p_acrophase=g.p_acrophase, passed=bool(g.p_pairing < .05 and g.p_acrophase < .05)))
    return rows


def by_donor(pb):
    return {d: g.reset_index(drop=True) for d, g in pb.groupby("donor")}


def setup(all_day=False, gtex_only=False):
    """Fetal design (n, depth) and the adult pools of every calibration job (`all_day`: GTEx deaths at any hour,
    the bound for fetal clocks not entrained to clock time; `gtex_only`: skip the nucleus atlases)."""
    org, _ = data.load("fetal_cao", level="organ")
    org = org[org.n_cells >= MIN_ORGAN_CELLS]
    cls, _ = data.load("fetal_cao", level="class")
    cls = cls[cls.n_cells >= MIN_CLASS_CELLS]
    design = []
    # 1a. GTEx bulk organs, 08-18 h deaths
    pb, _ = gtex.load(list(GTEX.values()), hardy=(0, 1, 2, 3, 4), depth="library_size")
    _, d = gtex.load(list(GTEX.values()), hardy=(0, 1, 2, 3, 4), depth="library_size")
    if not all_day:
        pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))]
    del d
    for o, t in GTEX.items():
        f = org[org.compartment == o]
        a = pb[pb.compartment == t].assign(compartment=o, stratum="one")
        if len(f) < 5 or a.donor.nunique() < 2 * len(f):
            continue
        JOBS.append(("L1", f"gtex:{o}", by_donor(a), len(f), {o: f.depth.median()}, [o]))
        design.append(dict(kind="L1", label=f"gtex:{o}", fetal_n=len(f), fetal_depth=f.depth.median(),
                           adult_pool=a.donor.nunique(), adult_depth=a.depth.median(), calibrated=True))
    # 1b and 2. adult nucleus atlases (one stratum: the one with most usable donors) per Cao class, and the classes
    # of an organ together
    for o, (name, cmap) in ([] if gtex_only else ATLAS.items()):
        a, _ = data.load(name)
        fc = cls[cls.compartment.isin([f"{o}:{c}" for c in cmap])].copy()
        fc["cls"] = fc.compartment.str.split(":").str[1]
        for c, ac in cmap.items():
            f = fc[fc.cls == c]
            x = a[a.compartment == ac]
            x = x[x.stratum == x.groupby("stratum").donor.nunique().idxmax()].assign(compartment=c, stratum="one")
            ok = len(f) >= 5 and x.donor.nunique() >= POOL_FACTOR * len(f)
            if ok:
                JOBS.append(("L1", f"{name}:{c}", by_donor(x), len(f), {c: f.depth.median()}, [c]))
            design.append(dict(kind="L1", label=f"{name}:{c}", fetal_n=len(f), fetal_depth=f.depth.median(),
                               adult_pool=x.donor.nunique(), adult_depth=x.depth.median(), calibrated=ok))
        use = [c for c in cmap if (fc.cls == c).sum() >= MIN_FETUSES_L2]
        n_f = int((fc[fc.cls.isin(use)].groupby("donor").cls.nunique() == len(use)).sum())
        x = a[a.compartment.isin([cmap[c] for c in use])].copy()
        x["compartment"] = x.compartment.map({cmap[c]: c for c in use})
        full = x.groupby(["stratum", "donor"]).compartment.nunique()
        full = full[full == len(use)].reset_index()
        st = full.groupby("stratum").donor.nunique().idxmax()
        x = x[(x.stratum == st) & x.donor.isin(full.donor[full.stratum == st])].assign(stratum="one")
        target = fc[fc.cls.isin(use)].groupby("cls").depth.median().to_dict()
        ok = len(use) >= 2 and n_f >= 5 and x.donor.nunique() >= POOL_FACTOR * n_f
        if ok:
            JOBS.append(("L2", f"{name}:{'+'.join(use)}", by_donor(x), n_f, target, use))
        design.append(dict(kind="L2", label=f"{name}:{'+'.join(use)}", fetal_n=n_f,
                           fetal_depth=";".join(f"{k}={v:.3g}" for k, v in target.items()),
                           adult_pool=x.donor.nunique(), adult_depth=x.depth.median(), calibrated=ok))
    return pd.DataFrame(design), org, cls


def setup_adrenal():
    """Proxy L2 calibration for the fetal adrenal (parenchyma, stromal, endothelial; no adult adrenal nucleus atlas):
    the heart and kidney atlas classes cut to the adrenal design (fetuses with all three classes, adrenal class
    depths)."""
    cls, _ = data.load("fetal_cao", level="class")
    cls = cls[cls.n_cells >= MIN_CLASS_CELLS]
    use = ["parenchyma", "stromal", "endothelial"]
    fc = cls[cls.compartment.isin([f"Adrenal:{c}" for c in use])].copy()
    fc["cls"] = fc.compartment.str.split(":").str[1]
    n_f = int((fc.groupby("donor").cls.nunique() == len(use)).sum())
    target = fc.groupby("cls").depth.median().to_dict()
    design = []
    for o in ("Heart", "Kidney"):
        name, cmap = ATLAS[o]
        a, _ = data.load(name)
        x = a[a.compartment.isin([cmap[c] for c in use])].copy()
        x["compartment"] = x.compartment.map({cmap[c]: c for c in use})
        full = x.groupby(["stratum", "donor"]).compartment.nunique()
        full = full[full == len(use)].reset_index()
        st = full.groupby("stratum").donor.nunique().idxmax()
        x = x[(x.stratum == st) & x.donor.isin(full.donor[full.stratum == st])].assign(stratum="one")
        ok = x.donor.nunique() >= POOL_FACTOR * n_f
        if ok:
            JOBS.append(("L2", f"adrenal_proxy:{name}", by_donor(x), n_f, target, use))
        design.append(dict(kind="L2", label=f"adrenal_proxy:{name}", fetal_n=n_f,
                           fetal_depth=";".join(f"{k}={v:.3g}" for k, v in target.items()),
                           adult_pool=x.donor.nunique(), adult_depth=x.depth.median(), calibrated=ok))
    return pd.DataFrame(design)


def markers(org, cls):
    """Detection of core and output clock genes in the fetal pseudobulks (counts only)."""
    rows = []
    for level, t in (("organ", org), ("class", cls)):
        for c, g in t.groupby("compartment"):
            if g.donor.nunique() < 5:
                continue
            r = dict(level=level, compartment=c, n_fetuses=g.donor.nunique(), median_depth=g.depth.median(),
                     ga_days_min=g.ga_days.min() if "ga_days" in g else np.nan,
                     ga_days_max=g.ga_days.max() if "ga_days" in g else np.nan)
            for gene in CORE + OUTPUT:
                r[f"{gene}_median_count"] = g[gene].median()
                r[f"{gene}_frac_ge20"] = (g[gene] >= 20).mean()
            rows.append(r)
    return pd.DataFrame(rows)


def expected_markers(cls):
    """Counts of each marker gene expected in a fetal class pseudobulk if it expressed the gene at the adult level:
    adult median fraction (gene count / depth, matching nucleus-atlas class) x fetal median class depth. An absence
    (< 20 counts) is informative only where this expectation is well above 20."""
    rows = []
    for o, (name, cmap) in ATLAS.items():
        a, _ = data.load(name)
        for c, ac in cmap.items():
            f = cls[cls.compartment == f"{o}:{c}"]
            x = a[a.compartment == ac]
            if len(f) < 5 or x.empty:
                continue
            r = dict(compartment=f"{o}:{c}", adult=f"{name}:{ac}", n_fetuses=len(f), fetal_depth=f.depth.median(),
                     adult_depth=x.depth.median())
            for g in CORE + OUTPUT:
                r[f"{g}_expected"] = (x[g] / x.depth).median() * f.depth.median()
                r[f"{g}_fetal_median"] = f[g].median()
            rows.append(r)
    return pd.DataFrame(rows)


def ci(k, n):
    r = binomtest(int(k), int(n)).proportion_ci(0.95, method="exact")
    return r.low, r.high


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--markers", action="store_true", help="only the expected-marker table")
    ap.add_argument("--all-day", action="store_true", help="GTEx L1 only, deaths at any hour (suffix _allday)")
    ap.add_argument("--adrenal", action="store_true", help="proxy L2 calibration at the adrenal design (_adrenal)")
    a = ap.parse_args()
    if a.markers:
        OUT.mkdir(exist_ok=True)
        cls, _ = data.load("fetal_cao", level="class")
        expected_markers(cls[cls.n_cells >= MIN_CLASS_CELLS]).to_csv(OUT / "markers_expected.csv", index=False)
        return
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    limit_threads(1)
    sfx = ("_check" if a.check else "") + ("_allday" if a.all_day else "") + ("_adrenal" if a.adrenal else "")
    if a.adrenal:
        design = setup_adrenal()
    else:
        design, org, cls = setup(a.all_day, gtex_only=a.all_day)
        meta = data.load("fetal_cao", level="organ")[1]
        if meta is not None and "ga_days" in meta:
            org, cls = (t.merge(meta[["ga_days"]], left_on="donor", right_index=True, how="left") for t in (org, cls))
        markers(org, cls).to_csv(OUT / f"markers{sfx}.csv", index=False)
    design.to_csv(OUT / f"design{sfx}.csv", index=False)
    print(design.to_string(index=False), flush=True)
    n = 3 if a.check else a.n_sub
    tasks = [(j, i) for j in range(len(JOBS)) for i in range(n)]
    print(f"{len(JOBS)} calibration jobs x {n} subsamples", flush=True)
    res = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, tasks, chunksize=4)):
            res.extend(r)
            if (k + 1) % 200 == 0:
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res)
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    rows = []
    for (kind, label, depth), r in res.groupby(["kind", "label", "depth"]):
        k, m = int(r.passed.sum()), len(r)
        lo, hi = ci(k, m)
        rows.append(dict(kind=kind, label=label, depth=depth, n=int(r.n.iloc[0]), n_subsamples=m, pass_rate=k / m,
                         ci_lo=lo, ci_hi=hi, median_stat=r.stat.median(), median_p=r.p.median(),
                         feasible=k / m >= 0.8))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False))
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
