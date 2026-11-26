"""A401: does a cancer type's bulk tumour-normal clock pairing follow from cell-level timekeeping and composition?
Power stage. Registered in paper/docs/registered_experiments.md (A401).

NO CANCER-CELL STATISTIC is computed here. Malignant pseudobulks are only counted and carried (genome tables) for the
score stage; every dataset that is fitted or gated contains the tumour STROMAL anchors only.

  1 pooled tumour table of three sources (unit: a malignant / tumour-epithelium pseudobulk plus >= 1 stromal anchor of
    >= 50 cells): (a) the A388 cohorts (`a388.cohort_table`, `min_anchor=1`, its deduplication), (b) the GEO builds
    `$SCCC_DATA/geo_tumour_a401/` (count units), (c) the CELLxGENE builds `$SCCC_DATA/cellxgene_tumour_a401/`
    (donor x compartment pseudobulks of the primary organ for matched cancers). A patient's samples are pooled.
  2 deduplication across sources by normalised patient id within a cancer family (priority A388 > GEO > Census) and by
    source batch (Global Atlas pal_2021 = Pal on disk)                                          -> dedup.csv, excluded.csv
  3 counts per cancer and source                                                               -> counts.csv, counts_source.csv
    matched cancers (>= 5 patients, TCGA code) and the exact one-sided critical Spearman rho   -> matched.csv
  4 P401 per cancer (anchors only): host-only gate (1,000 permutations, seed 401); leave-one-out coupling of each
    stromal class on the pooled host of the cancer's other stromal classes, standardised against 200 matched sets
    (`clocksync.paired.pseudo_waveforms`, seeds [401, 0], drawn on the cancer's pooled anchors, genes measured in every
    cohort of the cancer); Stouffer z over the classes, one-sided p                             -> p401.csv, p401_classes.csv
  5 TCGA mean ABSOLUTE purity of the matched cancers                                           -> purity.csv
  6 decidability: matched cancers passing P401 (>= 7 needed)                                   -> decidability.csv

    python a401.py --stage power [--jobs 4] [--genome-jobs 2] [--smoke]       (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import itertools
import re
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(EXP / "A388_malignant_tumour_host/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a388 as B  # noqa: E402
import a336_common as C  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

OUT = HERE.parents[1] / "output"
CACHE = OUT / "cache"
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
GEO = DATA / "geo_tumour_a401"
CEN = DATA / "cellxgene_tumour_a401"
PURITY = DATA / "tcga_purity_absolute" / "mastercalls.txt"
TCGA_CACHE = PAPER / "results/figure_5/output/cancer/cache"
BASE = C.WF
CORE = C.core(BASE)
READOUTS = {"14": BASE, "core": CORE}
STROMA = ("fibroblast", "endothelial", "mural")
ANCH = [f"tumour:{k}" for k in STROMA]
MALT = B.MALT
SEED = 401
META = B.META
MIN_CELLS = 50
MATCHED_TCGA = {"BRCA": "BRCA", "COAD": "COAD", "LUAD": "LUAD", "LIHC": "LIHC", "PRAD": "PRAD", "KIRC": "KIRC", "KIRP": "KIRP",
                "THCA": "THCA", "BLCA": "BLCA", "UCEC": "UCEC"}
FAMILY = {"LUAD": "lung", "lung_mixed": "lung", "LUSC": "lung", "lung_other": "lung", "SCLC": "lung"}
GEO_CANCER = {"GSE207493": "KIRC", "GSE163203": "THCA", "GSE191288": "THCA", "GSE135337": "BLCA", "GSE202051": "PAAD",
              "GSE155698": "PAAD", "GSE148071": "lung_mixed", "GSE125449": "LIHC", "GSE141445": "PRAD", "GSE197461": "CESC",
              "GSE212527": "SARC"}
DROP_CENSUS = {"allcells_16023185": "copy of GSE178341", "asinglecellandspat_9fddb063": "copy of Wu 2021 (GSE176078)",
               "amultitissuesingle_b617ee1b": "A226 donors (Qian 2020 CRC / LUAD)", "anintegratedcellat_9f222629": "HLCA",
               "extendedgbmap_56c4912d": "GBmap; glioma has no fibroblasts", "singlecellrnaseque_b6b5ea88": "OPSCC, A226",
               "singlecelltranscri_5af90777": "RCC, A226", "freshsamples_6c87755e": "sibling frozensamples_34f5307e (same patients)"}
ORGAN = {"BRCA": {"breast"}, "COAD": {"colon", "large intestine", "intestine"}, "LUAD": {"lung"}, "KIRC": {"kidney"},
         "KIRP": {"kidney"}, "PRAD": {"prostate gland"}, "LIHC": {"liver"}}
A388_CANCER = {"breast": "BRCA", "colorectal": "COAD", "liver": "LIHC", "prostate": "PRAD", "thyroid": "THCA", "bladder": "BLCA",
               "esophageal": "ESCA", "gastric": "STAD", "pancreatic": "PAAD", "ovarian": "OV", "cervical": "CESC",
               "head and neck": "HNSC", "kidney": "kidney_unspecified", "cancer (unspecified)": "unspecified",
               "ampullary": "ampullary", "skin": "SKCM"}
GENERIC_ID = re.compile(r"^(p|pt|patient|case|sample|s|n|t|pat|donor)?\d{1,3}[a-z]?$")


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def norm_id(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def family(c):
    return FAMILY.get(c, c)


def census_cancer(disease):
    d = str(disease).lower()
    if "metastatic" in d:
        return None
    if "non-small" in d:
        return "lung_other"
    if "nonpapillary" in d:
        return "KIRC"
    for key, c in (("breast", "BRCA"), ("lung adenocarcinoma", "LUAD"), ("small cell lung", "SCLC"), ("squamous cell lung", "LUSC"),
                   ("lung", "lung_other"), ("colon", "COAD"), ("colorectal", "COAD"), ("papillary renal", "KIRP"),
                   ("renal cell", "KIRC"), ("prostat", "PRAD"), ("hepatocellular", "LIHC"), ("liver cancer", "LIHC"),
                   ("ovarian", "OV"), ("gastric", "STAD"), ("oropharynx", "HNSC"), ("bladder", "BLCA"), ("urothelial", "BLCA"),
                   ("endometri", "UCEC"), ("thyroid", "THCA")):
        if key in d:
            return c
    return None


def tidy(G, cohort_prefix=None):
    """Canonical gene symbols (duplicates summed) of a genome-wide pseudobulk frame."""
    genes = [c for c in G.columns if c not in META and c not in ("sample", "median_umi", "units", "metadata")]
    X = G[genes].copy()
    X.columns = [canonical(str(g)) for g in X.columns]
    X = X.T.groupby(level=0).sum().T.astype(np.float32)
    return pd.concat([G[[c for c in META if c in G]].reset_index(drop=True), X.reset_index(drop=True)], axis=1)


def pool(G):
    """One row per donor x compartment (summing counts, n_cells and depth over samples / sites / assays)."""
    genes = [c for c in G.columns if c not in META]
    agg = {c: "sum" for c in ["n_cells", "depth"] + genes}
    st = G.groupby(["donor", "compartment"]).stratum.first()
    P = G.groupby(["donor", "compartment"], sort=False).agg(agg).reset_index()
    P["stratum"] = [st[(d, c)] for d, c in zip(P.donor, P.compartment)]
    P[genes] = P[genes].astype(np.float32)
    return P[META + genes]


def finish(G, source, cancer_of):
    """Rows of the tumour malignant (renamed) and stromal compartments with >= 50 cells; donors with a malignant and >= 1
    anchor. Returns the genome-wide frame (or None) and the donor table."""
    G = G[G.compartment.isin(["tumour:epithelial", MALT] + ANCH) & (G.n_cells >= MIN_CELLS)].copy()
    G["compartment"] = G.compartment.where(G.compartment != "tumour:epithelial", MALT)
    G = pool(G)
    na = G[G.compartment.isin(ANCH)].groupby("donor").compartment.nunique()
    ok = sorted(d for d in set(G[G.compartment == MALT].donor) if na.get(d, 0) >= 1 and cancer_of.get(d) is not None)
    G = G[G.donor.isin(ok)].reset_index(drop=True)
    don = pd.DataFrame({"donor": ok, "source": source, "cancer": [cancer_of[d] for d in ok]})
    return G, don


def build_geo(g):
    d = pd.read_csv(GEO / g / "donors.csv")
    G = tidy(pd.read_parquet(GEO / g / "pseudobulk.parquet"))
    cancer = {f"{g}|{p}": GEO_CANCER[g] for p in d.patient.astype(str)}
    if g == "GSE125449":
        cancer = {f"{g}|{p}": ("LIHC" if s == "HCC" else None) for p, s in zip(d.patient.astype(str), d.subtype)}
    return finish(G, g, cancer)


def build_census(s):
    G = pd.read_parquet(CEN / f"{s}.parquet")
    G = G[G.compartment.str.startswith("tumour:")]
    dis = G.groupby("donor").disease.agg(lambda x: x.value_counts().index[0])
    cancer = {d: census_cancer(v) for d, v in dis.items()}
    org = G.groupby("donor").tissue_general.agg(lambda x: set(x))
    keep = []
    for d, c in cancer.items():
        keep.append(d)
    G = G[G.donor.isin(keep)]
    # primary-organ rows for matched cancers (metastasis sites dropped)
    c_row = G.donor.map(cancer)
    ok = np.array([(c not in ORGAN) or (t in ORGAN[c]) for c, t in zip(c_row, G.tissue_general)])
    G = G[ok]
    G = G.rename(columns={"stratum": "stratum"})
    meta_drop = [c for c in ("disease", "tissue", "tissue_general", "suspension", "tcga", "n_cells_sample", "median_umi") if c in G]
    G = tidy(G.drop(columns=meta_drop))
    return finish(G, s, cancer)


def build_a388(args):
    studies, inv = study_list(args)
    res = Parallel(n_jobs=args.genome_jobs)(delayed(B.cohort_table)(s, None, 1) for s in studies)
    err = {s: r for s, r in zip(studies, res) if isinstance(r, str)}
    tabs = {s: r[0] for s, r in zip(studies, res) if not isinstance(r, str) and len(r[0])}
    allp = {s: r[1] for s, r in zip(studies, res) if not isinstance(r, str)}
    mal = {s: r[2] for s, r in zip(studies, res) if not isinstance(r, str)}
    removed, drec, _ = B.dedup(tabs, allp, inv, set(inv.index[inv.in_H138_1.astype(bool)]))
    for s in tabs:
        tabs[s] = tabs[s][[(s, d) not in removed for d in tabs[s].donor]]
    tabs = {s: t for s, t in tabs.items() if len(t)}
    keeps = {s: set(zip(t.donor, t.compartment)) for s, t in tabs.items()}
    CACHE.mkdir(parents=True, exist_ok=True)
    gres = Parallel(n_jobs=args.genome_jobs)(delayed(B.genome_cohort)(s, mal[s], sorted(t.donor.unique()), keeps[s], str(CACHE))
                                             for s, t in tabs.items())
    out, don, dropped = {}, [], [dict(source=s, stage="table", error=e) for s, e in err.items()]
    for s, g, e in gres:
        if g is None:
            dropped.append(dict(source=s, stage="genome", error=e))
            continue
        g = g.copy()
        if "n_cells" not in g:
            g = g.merge(tabs[s][["donor", "compartment", "n_cells"]], on=["donor", "compartment"], how="left")
        g["donor"] = s + "|" + g.donor.astype(str)
        g["stratum"] = s + "|" + g.stratum.astype(str)
        cancer = inv.cancer[s]
        if cancer == "lung":
            pass
        out[s] = g
        for d in g.donor.unique():
            don.append(dict(donor=d, source=s, cancer=cancer, cohort_in_H138_1=bool(inv.in_H138_1[s])))
    don = pd.DataFrame(don)
    don["cancer"] = [("LUAD" if s in B.LUAD_COHORTS else "lung_mixed") if c == "lung" else A388_CANCER.get(c, c)
                     for c, s in zip(don.cancer, don.source)]
    return out, don, drec, dropped


def study_list(args):
    inv = pd.read_csv(B.FEAS).set_index("study")
    inv = inv[(inv.n_malignant > 0) & ((inv.n_fibroblast + inv.n_endothelial + inv.n_mural) > 0)]
    studies = list(inv.index)
    if args.smoke:
        studies = [s for s in ("SRP416180", "GSE222315", "GSE281736", "GSE164789") if s in studies]
        inv = inv.loc[studies]
    return studies, inv


# ---- counts, matched cancers ---------------------------------------------------------------------------------------

def critical_rho(n, alpha=0.05):
    """Exact one-sided critical Spearman rho (smallest rho with P(rho >= r) <= alpha) over all n! rankings."""
    from scipy.stats import rankdata  # noqa: F401
    x = np.arange(n)
    chunks = []
    for p in itertools.permutations(range(n)):
        chunks.append(np.sum((x - np.array(p)) ** 2))
    d2 = np.array(chunks)
    rho = 1 - 6 * d2 / (n * (n * n - 1))
    vals = np.sort(np.unique(np.round(rho, 12)))[::-1]
    best = np.nan
    for r in vals:
        if np.mean(rho >= r - 1e-12) > alpha + 1e-12:
            break
        best = r
    return float(best)


def counts_tables(donors, G_by_donor):
    rows = []
    for c, d in donors.groupby("cancer"):
        ids = set(d.donor)
        cl = {k: len([x for x in ids if G_by_donor.get((x, f"tumour:{k}"), False)]) for k in STROMA}
        n_anchor = pd.Series({x: sum(bool(G_by_donor.get((x, f"tumour:{k}"), False)) for k in STROMA) for x in ids})
        rows.append(dict(cancer=c, patients=len(ids), cohorts=d.source.nunique(), n_fibroblast=cl["fibroblast"],
                         n_endothelial=cl["endothelial"], n_mural=cl["mural"], patients_ge2_anchors=int((n_anchor >= 2).sum()),
                         patients_3_classes=int((n_anchor >= 3).sum()), tcga=MATCHED_TCGA.get(c, "")))
    return pd.DataFrame(rows).sort_values("patients", ascending=False)


# ---- P401 -------------------------------------------------------------------------------------------------------------

def fit_classes(Gc, classes, wf):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    t = cs.Dataset(Gc, anchors=classes, waveform=wf, mixed_strata=True).fit(n_boot=0, gate=False, targets=classes).table
    return {k: float(t.loc[k, "coupling"]) for k in classes}


def comp_chunk(Gc, classes, ws, rd):
    out = []
    for i, w in ws:
        try:
            out.append(dict(set=i, readout=rd, **fit_classes(Gc, classes, w)))
        except Exception as e:  # noqa: BLE001
            out.append(dict(set=i, readout=rd, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def clock_task(table_c, classes, rd, n_perm):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    wf = READOUTS[rd]
    ds = cs.Dataset(table_c, anchors=classes, waveform=wf, mixed_strata=True)
    try:
        g = ds.gate(n_perm=n_perm, seed=SEED)
        gate = dict(gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase, gate_clock_share=g.clock_share,
                    gate_share_null_mean=g.share_null_mean, gate_passes=bool(g.passes), gate_error="")
    except Exception as e:  # noqa: BLE001
        gate = dict(gate_passes=False, gate_error=f"{type(e).__name__}: {e}"[:200])
    try:
        c = fit_classes(table_c, classes, wf)
    except Exception as e:  # noqa: BLE001
        c, gate["fit_error"] = {}, f"{type(e).__name__}: {e}"[:200]
    return rd, gate, c


def stage_p401(cancers, frames, donors, args, out):
    n_set = 6 if args.smoke else 200
    n_perm = 100 if args.smoke else 1000
    from clocksync.paired import pseudo_waveforms
    rows, crow = [], []
    for c in cancers:
        ids = set(donors[donors.cancer == c].donor)
        parts = [f[f.donor.isin(ids) & f.compartment.isin(ANCH)] for f in frames.values()]
        parts = [p for p in parts if len(p)]
        Gc = pd.concat(parts, join="inner", ignore_index=True)  # genes measured in every cohort of the cancer
        genes = [g for g in Gc.columns if g not in META]
        classes = [a for a in ANCH if Gc[Gc.compartment == a].donor.nunique() >= 5]
        base = dict(cancer=c, patients=len(ids), cohorts=donors[donors.cancer == c].source.nunique(), genes_all_cohorts=len(genes),
                    classes="+".join(k.split(":")[1] for k in classes))
        if len(classes) < 3 or not set(B.BASE.genes) <= set(genes):
            for rd in READOUTS:
                rows.append(dict(base, readout=rd, p401_pass=False, status="not estimable: < 3 stromal classes with >= 5 patients"
                                 if len(classes) < 3 else "not estimable: clock genes missing"))
            log(f"{c}: not estimable ({base['classes']})")
            continue
        Gc = Gc[Gc.compartment.isin(classes)].reset_index(drop=True)
        tab = Gc[META + list(BASE.genes)]
        ref = Gc
        frac = ref[genes].div(ref.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        ws = list(enumerate(pseudo_waveforms(lm, n_set, np.random.default_rng([SEED, 0]), BASE, "matched")))
        used = sorted({g for _, w in ws for g in w.genes} | set(BASE.genes))
        Gu = Gc[META + used]
        clk = Parallel(n_jobs=min(2, args.jobs))(delayed(clock_task)(tab, classes, rd, n_perm) for rd in READOUTS)
        tasks = []
        for rd in READOUTS:
            w_rd = [(i, w if rd == "14" else C.restrict(w, BASE, CORE)) for i, w in ws]
            for k in range(args.jobs * 2):
                if w_rd[k::args.jobs * 2]:
                    tasks.append((Gu, classes, w_rd[k::args.jobs * 2], rd))
        comp = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(comp_chunk)(*t) for t in tasks) for r in x])
        for rd, gate, cc in clk:
            m = comp[comp.readout == rd]
            zs, ps = [], []
            for k in classes:
                v = m[k].dropna().to_numpy() if k in m else np.array([])
                if k not in cc or len(v) < min(10, n_set):
                    continue
                z = (cc[k] - v.mean()) / v.std(ddof=1)
                zs.append(z)
                crow.append(dict(cancer=c, readout=rd, cls=k.split(":")[1], coupling=cc[k], matched_mean=float(v.mean()),
                                 matched_sd=float(v.std(ddof=1)), matched_median=float(np.median(v)),
                                 matched_q95=float(np.quantile(v, .95)), z=z, p_one_sided=float(norm.sf(z)),
                                 share_matched_ge=float(np.mean(v >= cc[k])), n_matched_ok=len(v)))
            Z = float(np.sum(zs) / np.sqrt(len(zs))) if zs else np.nan
            p = float(norm.sf(Z)) if zs else np.nan
            rows.append(dict(base, readout=rd, n_classes_scored=len(zs), z_pooled=Z, p_pooled=p, **gate,
                             p401_pass=bool(gate.get("gate_passes", False) and np.isfinite(p) and p < .05),
                             status="ok" if zs else "no class fitted"))
        log(f"{c}: P401 done ({base['classes']})")
        pd.DataFrame(rows).to_csv(out / "p401.csv", index=False)
    P = pd.DataFrame(rows)
    P.to_csv(out / "p401.csv", index=False)
    pd.DataFrame(crow).to_csv(out / "p401_classes.csv", index=False)
    return P


# ---- purity -------------------------------------------------------------------------------------------------------------

def purity_table(matched):
    m = pd.read_csv(PURITY, sep="\t")
    m["case"] = m["array"].astype(str).str[:12]
    m = m[m["call status"] == "called"]
    p = m.drop_duplicates("case").set_index("case").purity
    rows = []
    for c in matched:
        f = TCGA_CACHE / f"tcga_tumours_{MATCHED_TCGA[c]}.parquet"
        if not f.exists():
            rows.append(dict(cancer=c, tcga=MATCHED_TCGA[c], error="no TCGA cache"))
            continue
        b = pd.read_parquet(f, columns=["barcode"]).barcode.astype(str).str[:12]
        v = p.reindex(b).dropna()
        rows.append(dict(cancer=c, tcga=MATCHED_TCGA[c], primaries=len(b), with_purity=len(v), mean_purity=float(v.mean()),
                         median_purity=float(v.median()), sd_purity=float(v.std())))
    return pd.DataFrame(rows)


# ---- main ---------------------------------------------------------------------------------------------------------------

def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    frames, dons, dropped = {}, [], []
    # (a) A388 cohorts
    f388, d388, drec, dr = build_a388(args)
    frames.update(f388)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)  # free the workers' memory before the genome-wide GEO / Census reads
    dons.append(d388.assign(priority=0))
    dropped += dr
    log(f"A388: {len(f388)} cohorts, {d388.donor.nunique()} patients (>= 1 anchor)")
    # (b) GEO builds
    geos = [g for g in GEO_CANCER if (GEO / g / "pseudobulk.parquet").exists()]
    if args.smoke:
        geos = geos[:2]
    for g in geos:
        try:
            G, d = build_geo(g)
            frames[g] = G
            dons.append(d.assign(priority=1))
            log(f"{g}: {d.donor.nunique()} patients")
        except Exception as e:  # noqa: BLE001
            dropped.append(dict(source=g, stage="geo", error=f"{type(e).__name__}: {e}"[:300]))
    # (c) Census
    cen = sorted(p.name[:-8] for p in CEN.glob("*.parquet") if p.name[:-8] not in DROP_CENSUS and not p.name.endswith("_merged.parquet"))
    if args.smoke:
        cen = [s for s in cen if s.startswith("singlecellatlas")][:1]
    excl = []
    for s in cen:
        try:
            G, d = build_census(s)
            if len(d):
                frames[s] = G
                dons.append(d.assign(priority=2))
            log(f"{s}: {d.donor.nunique()} patients ({', '.join(sorted(d.cancer.unique()))})")
        except Exception as e:  # noqa: BLE001
            dropped.append(dict(source=s, stage="census", error=f"{type(e).__name__}: {e}"[:300]))
    for s, why in DROP_CENSUS.items():
        excl.append(dict(source=s, reason=why))
    # a source lacking a clock gene (all 14 are needed by both readouts) is dropped from every test and listed
    for s in list(frames):
        miss = [g for g in BASE.genes if g not in frames[s].columns]
        if miss:
            dropped.append(dict(source=s, stage="clock genes", error=f"clock gene(s) not measured: {', '.join(miss)}"))
            del frames[s]
    dons = [d[d.source.isin(frames)] for d in dons]
    pd.DataFrame(excl).to_csv(out / "excluded.csv", index=False)
    pd.DataFrame(dropped, columns=["source", "stage", "error"]).to_csv(out / "dropped.csv", index=False)
    donors = pd.concat(dons, ignore_index=True)
    donors["id"] = donors.donor.str.split("|", n=1).str[1]
    donors["nid"] = donors.id.map(norm_id)
    donors["family"] = donors.cancer.map(family)
    # deduplication by id within a cancer family (priority A388 > GEO > Census) and by source batch
    rec = []
    drop = set()
    ga = pd.read_csv(CEN / "globalatlas_de5416ef_donors.csv")
    pal = set(ga[ga.src_batch.astype(str).isin(["pal_2021", "wu_2021"])].donor)
    for d in donors.donor:
        if d in pal:
            drop.add(d)
            rec.append(dict(removed=d, kept="", cancer="BRCA", reason="Global Atlas donor from Pal / Wu 2021 (src_batch): on-disk copy kept"))
    # the same, for a patient present in more than one source
    ok = donors[~donors.donor.isin(drop)]
    for (fam, nid), g in ok.groupby(["family", "nid"]):
        if len(g) < 2 or GENERIC_ID.match(nid) or len(nid) < 5:
            continue
        g = g.sort_values(["priority", "donor"])
        for _, r in g.iloc[1:].iterrows():
            if r.source == g.iloc[0].source:
                continue
            drop.add(r.donor)
            rec.append(dict(removed=r.donor, kept=g.iloc[0].donor, cancer=fam, reason="same normalised patient id in another source"))
    generic = ok[ok.duplicated(["family", "nid"], keep=False) & ok.nid.map(lambda x: bool(GENERIC_ID.match(x)) or len(x) < 5)]
    rec.append(dict(removed="", kept="", cancer="", reason=f"{generic.donor.nunique()} donors share a generic id (e.g. 'P1') with another source "
                                                        f"of the same cancer family: not treated as overlap (not checkable by id)"))
    pd.DataFrame(rec, columns=["removed", "kept", "cancer", "reason"]).to_csv(out / "dedup.csv", index=False)
    log(f"dedup: {len(drop)} patients removed")
    donors = donors[~donors.donor.isin(drop)].reset_index(drop=True)
    for s in list(frames):
        frames[s] = frames[s][~frames[s].donor.isin(drop)]
    # counts
    has = {(d, c): True for f in frames.values() for d, c in zip(f.donor, f.compartment)}
    counts = counts_tables(donors, has)
    counts.to_csv(out / "counts.csv", index=False)
    cs_rows = [dict(source=s, cancer=c, patients=g.donor.nunique()) for (s, c), g in donors.groupby(["source", "cancer"])]
    pd.DataFrame(cs_rows).to_csv(out / "counts_source.csv", index=False)
    pd.DataFrame([dict(source=s, path=str(CACHE / f"genome_{s}.parquet") if s in f388 else
                       str(GEO / s / "pseudobulk.parquet") if s in GEO_CANCER else str(CEN / f"{s}.parquet"))
                  for s in frames]).to_csv(out / "genome_paths.csv", index=False)
    matched = [c for c in counts.cancer if c in MATCHED_TCGA and counts.set_index("cancer").patients[c] >= 5]
    nm = len(matched)
    mt = pd.DataFrame([dict(cancer=c, tcga=MATCHED_TCGA[c], patients=int(counts.set_index("cancer").patients[c]), matched=c in matched)
                       for c in MATCHED_TCGA if c in set(counts.cancer)])
    mt["n_matched"] = nm
    mt["critical_rho_one_sided_05"] = critical_rho(nm) if 3 <= nm <= 10 else np.nan
    mt.to_csv(out / "matched.csv", index=False)
    pd.DataFrame([dict(n=n, critical_rho=critical_rho(n)) for n in range(5, 10)]).to_csv(out / "critical_rho.csv", index=False)
    log("counts\n" + counts.to_string(index=False))
    # P401 for every cancer with >= 5 patients
    cancers = [c for c in counts.cancer if counts.set_index("cancer").patients[c] >= 5]
    P = stage_p401(cancers, frames, donors, args, out)
    purity_table(matched).to_csv(out / "purity.csv", index=False)
    dec = []
    for rd in READOUTS:
        p = P[(P.readout == rd) & P.cancer.isin(matched)]
        ps = int(p.p401_pass.sum())
        dec.append(dict(readout=rd, matched_cancers=nm, matched_passing_p401=ps, needed=7, decidable=bool(ps >= 7),
                        passing="+".join(p[p.p401_pass].cancer)))
    pd.DataFrame(dec).to_csv(out / "decidability.csv", index=False)
    pd.set_option("display.width", 250)
    for n in ("matched", "p401", "p401_classes", "purity", "decidability"):
        log(n + "\n" + pd.read_csv(out / f"{n}.csv").round(4).to_string(index=False))
    log(f"power done {(time.time() - t0) / 60:.1f} min")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power",), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    stage_power(args)


if __name__ == "__main__":
    main()
