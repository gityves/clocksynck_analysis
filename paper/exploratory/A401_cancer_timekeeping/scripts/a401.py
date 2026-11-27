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

sys.path.insert(0, str(EXP / "A138_cancer_estimability/scripts"))
sys.path.insert(0, str(PAPER / "results/figure_5/scripts/cancer"))
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
ORIGINAL = {"BRCA": "BRCA", "COAD": "COAD", "LUAD": "LUAD", "LIHC": "LIHC", "PRAD": "PRAD", "KIRC": "KIRC", "KIRP": "KIRP",
            "THCA": "THCA", "BLCA": "BLCA", "UCEC": "UCEC"}
ADDED = {"STAD": "STAD", "ESCA": "ESCA", "LUSC": "LUSC", "HNSC": "HNSC"}  # note 1: bulk pairing computed here (blind)
MATCHED_TCGA = dict(ORIGINAL)  # extended in stage_power with the added projects that pass the TCGA gate with >= 10 pairs
FAMILY = {"LUAD": "lung", "lung_mixed": "lung", "LUSC": "lung", "lung_other": "lung", "SCLC": "lung"}
BASSEZ = DATA / "bassez2021_3ca" / "pseudobulk.parquet"
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


def tidy(G, cohort_prefix=None, only=None):
    """Canonical gene symbols (duplicates summed) of a genome-wide pseudobulk frame; `only`: keep only these canonical genes."""
    genes = [c for c in G.columns if c not in META and c not in ("sample", "median_umi", "units", "metadata")]
    if only is not None:
        genes = [c for c in genes if canonical(str(c)) in only]
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


def finish(G, source, cancer_of, all_comps=False):
    """Rows of the tumour malignant (renamed) and stromal compartments with >= 50 cells; donors with a malignant and >= 1
    anchor. Returns the genome-wide frame (or None) and the donor table."""
    keep = G.compartment.str.startswith("tumour:") & ~G.compartment.str.endswith(":unassigned") if all_comps \
        else G.compartment.isin(["tumour:epithelial", MALT] + ANCH)
    G = G[keep & (G.n_cells >= MIN_CELLS)].copy()
    G["compartment"] = G.compartment.where(G.compartment != "tumour:epithelial", MALT)
    G = pool(G)
    na = G[G.compartment.isin(ANCH)].groupby("donor").compartment.nunique()
    ok = sorted(d for d in set(G[G.compartment == MALT].donor) if na.get(d, 0) >= 1 and cancer_of.get(d) is not None)
    G = G[G.donor.isin(ok)].reset_index(drop=True)
    don = pd.DataFrame({"donor": ok, "source": source, "cancer": [cancer_of[d] for d in ok]})
    return G, don


def build_geo(g, all_comps=False, only=None):
    d = pd.read_csv(GEO / g / "donors.csv")
    G = tidy(pd.read_parquet(GEO / g / "pseudobulk.parquet"), only=only)
    cancer = {f"{g}|{p}": GEO_CANCER[g] for p in d.patient.astype(str)}
    if g == "GSE125449":
        cancer = {f"{g}|{p}": ("LIHC" if s == "HCC" else None) for p, s in zip(d.patient.astype(str), d.subtype)}
    return finish(G, g, cancer, all_comps)


def build_census(s, all_comps=False, only=None):
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
    G = tidy(G.drop(columns=meta_drop), only=only)
    return finish(G, s, cancer, all_comps)


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


# ---- bulk side of the added projects (note 1 (a)) -------------------------------------------------------------------------

def bulk_pairing(out, added):
    """A138 H138.4 statistic and gate of the normals for the added TCGA projects, by the A138 path (`_shared.tcga_builder`
    build, `pair_frames`, `q3_lib.reference_gate` rule: n_gate 200, rng [52, 0, i]; statistic `corr_stat`; pairing null
    `perm_null`, 10,000, rng [138, 4, 0, i]; p as `score.perm_p`). BRCA (i = 1) is recomputed as the reproduction check
    against `A138 .../score/per_unit.csv`. Added projects take i = 10..13 so that the original streams are unchanged."""
    import q3_lib as q
    from clocksync.paired import cohort_states, reference_gate
    import _shared as sh
    ref = pd.read_csv(EXP / "A138_cancer_estimability/output/score/per_unit.csv").query("hypothesis == 'H138.4'").set_index("unit")
    b = sh.tcga_builder()
    rows = []
    for c, i in [("BRCA", 1)] + [(p, 10 + j) for j, p in enumerate(ADDED) if p in added]:
        samples, _ = b.build([c])
        case, rf = sh.pair_frames(samples, c)
        r = dict(project=c, kind="reproduction" if c == "BRCA" else "added", pairs=len(case))
        if len(case) < 10:
            rows.append(dict(r, error="< 10 pairs"))
            continue
        yr = np.log(np.maximum(rf[q.GENES].astype(float), 1e-300)).sub(np.log(rf.depth.to_numpy(float)), axis=0).to_numpy()
        g = reference_gate(yr, q.WF, 200, np.random.default_rng([52, 0, i]))
        bt, vt, bn, vn = cohort_states(case, rf, True, q.WF)
        za, zb = q.centred_complex(bt), q.centred_complex(bn)
        obs = q.corr_stat(za, vt, zb, vn)
        null = q.perm_null(za, vt, zb, vn, 10_000, np.random.default_rng([138, 4, 0, i]))
        x = null[np.isfinite(null)]
        r.update(gate_p=g["p"], gated=bool(g["p"] < 0.05), statistic=obs, null_mean=float(np.nanmean(null)),
                 null_q95=float(np.quantile(x, .95)), p=float((1 + np.sum(x >= obs)) / (len(x) + 1)))
        if c == "BRCA":
            r.update(published_statistic=ref.statistic[c], published_null_mean=ref.null_mean[c], published_p=ref.p[c],
                     published_gate_p=ref.gate_p[c], reproduced=bool(np.isclose(obs, ref.statistic[c], atol=1e-9)))
        rows.append(r)
        log(f"bulk pairing {c}: {len(case)} pairs, gate p {r['gate_p']:.4f}, statistic {obs:.4f}, p {r['p']:.4f}")
    t = pd.DataFrame(rows)
    t.to_csv(out / "bulk_pairing_added.csv", index=False)
    return t


def a83_share(out, added):
    """A83 (a83 part b, `compute_tumour_gate.py`): clock-plane share of the 14 centred log fractions of ALL primary tumours
    of a project (`reference_gate`, 200 permuted-acrophase planes); BRCA recomputed as the reproduction check. The share is
    deterministic; the null (p, share_null) uses its own stream because A83's rng ran on after the paired loop."""
    import _shared as sh
    from clocksync.paired import log_expr, reference_gate
    ref = pd.read_csv(PAPER / "results/figure_5/output/cancer/a83_all_tumours.csv").set_index("project")
    b = sh.tcga_builder()
    rows = []
    for k, c in enumerate(["BRCA"] + [p for p in ADDED if p in added]):
        tum = b.tumours(c)
        y = log_expr(tum[sh.GENES].astype(float), tum.depth.to_numpy(float)).to_numpy()
        g = reference_gate(y, sh.WF, 200, np.random.default_rng([401, 83, k]))
        r = dict(project=c, kind="reproduction" if c == "BRCA" else "added", tumours=len(tum), share=g["share"],
                 share_null=g["share_null_mean"], p=g["p"])
        if c == "BRCA":
            r.update(published_share=ref.share[c], published_tumours=int(ref.tumours[c]), published_p=ref.p[c],
                     reproduced=bool(np.isclose(g["share"], ref.share[c], atol=1e-9) and len(tum) == ref.tumours[c]))
        rows.append(r)
        log(f"A83 share {c}: {len(tum)} tumours, share {g['share']:.4f}")
    t = pd.DataFrame(rows)
    t.to_csv(out / "a83_added.csv", index=False)
    return t


def build_bassez(all_comps=False, only=None):
    """BIOKEY (Bassez 2021, 3CA copy with ARNTL): pre-treatment samples only, donor `bassez2021|<sample>`."""
    G = pd.read_parquet(BASSEZ)
    G = G[G.timepoint == "Pre"]
    G = G.drop(columns=[c for c in ("timepoint", "expansion", "cohort", "BC_type") if c in G])
    G = tidy(G, only=only)
    don = sorted(set(G.donor))
    return finish(G, "bassez2021", {d: "BRCA" for d in don}, all_comps)


# ---- two-class P401 (note 1 (b)) -----------------------------------------------------------------------------------------

def rc_units(frame_c, classes):
    """Per cohort (source): anchors-only table and the patients having both classes (>= 4)."""
    out = []
    a, b = classes
    for s, f in frame_c.items():
        da = set(f[f.compartment == a].donor)
        db = set(f[f.compartment == b].donor)
        pts = sorted(da & db)
        if len(pts) >= 4:
            out.append((s, f, pts))
    return out


def rc_pooled(units, classes, wf, cols):
    """Complex-pooled |sum_k sum_i z_a conj(z_b)| / sqrt(sum |z_a|^2 sum |z_b|^2) (A387 Part 1 form) over cohorts; states
    on one projector per cohort (`a149_lib.pair_states`, centred within stratum, then within cohort)."""
    import a149_lib as L
    num, ea, eb, n = 0j, 0.0, 0.0, 0
    for s, f, pts in units:
        pb = f[f.compartment.isin(classes)][META + cols]
        p = L.pair_states(pb, classes[0], classes[1], wf=wf, pts=pts)
        if p is None:
            continue
        za, zb = p["za"] - p["za"].mean(), p["zb"] - p["zb"].mean()
        num += np.sum(za * np.conj(zb))
        ea += float(np.sum(np.abs(za) ** 2))
        eb += float(np.sum(np.abs(zb) ** 2))
        n += len(za)
    return (abs(num) / np.sqrt(ea * eb) if ea > 0 and eb > 0 else np.nan), n


def rc_chunk(units, classes, ws, rd):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in ws:
        try:
            r, n = rc_pooled(units, classes, w, list(w.genes))
            out.append(dict(set=i, readout=rd, R=r, n=n))
        except Exception as e:  # noqa: BLE001
            out.append(dict(set=i, readout=rd, error=f"{type(e).__name__}: {e}"[:200]))
    return out


# ---- counts, matched cancers ---------------------------------------------------------------------------------------

def critical_rho(n, alpha=0.05):
    """Exact one-sided critical Spearman rho (smallest rho with P(rho >= r) <= alpha) over all n! rankings."""
    if n > 9:  # exact enumeration is infeasible: Monte Carlo over 4,000,000 random rankings (seed 401), upper-tail rule as below
        rng = np.random.default_rng([SEED, n])
        x = np.arange(n)
        d2 = np.concatenate([np.sum((x[None, :] - np.argsort(rng.random((200_000, n)), axis=1)) ** 2, axis=1) for _ in range(20)])
        rho = 1 - 6 * d2 / (n * (n * n - 1))
        vals = np.sort(np.unique(np.round(rho, 12)))[::-1]
        best = np.nan
        for r in vals:
            if np.mean(rho >= r - 1e-12) > alpha + 1e-12:
                break
            best = r
        return float(best)
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


def two_class(c, base, classes, frames, ids, tab, Gu, ws, n_perm, args, rows, crow):
    """P401 of a cancer with exactly two stromal classes: anchors-only gate (two compartments) and the complex-pooled pairwise
    coupling of the classes against the same statistic of the 200 matched sets (z, one-sided p)."""
    for rd in READOUTS:
        wf = READOUTS[rd]
        try:
            g = cs.Dataset(tab, anchors=classes, waveform=wf, mixed_strata=True).gate(n_perm=n_perm, seed=SEED)
            gate = dict(gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase, gate_clock_share=g.clock_share,
                        gate_share_null_mean=g.share_null_mean, gate_passes=bool(g.passes), gate_error="")
        except Exception as e:  # noqa: BLE001
            gate = dict(gate_passes=False, gate_error=f"{type(e).__name__}: {e}"[:200])
        fc = {s: f[f.donor.isin(ids) & f.compartment.isin(classes)] for s, f in frames.items()}
        fc = {s: f for s, f in fc.items() if len(f)}
        clock_units = rc_units({s: f[META + list(BASE.genes)] for s, f in fc.items()}, classes)
        match_units = rc_units({s: f[[x for x in f.columns if x in set(META) | set(Gu.columns)]] for s, f in fc.items()}, classes)
        if not clock_units:
            rows.append(dict(base, readout=rd, form="pairwise", **gate, p401_pass=False, status="not estimable: no cohort with >= 4 pairs"))
            continue
        r_clock, n_pairs = rc_pooled(clock_units, classes, wf, list(wf.genes))
        w_rd = [(i, w if rd == "14" else C.restrict(w, BASE, CORE)) for i, w in ws]
        tasks = [(match_units, classes, w_rd[k::args.jobs * 2], rd) for k in range(args.jobs * 2) if w_rd[k::args.jobs * 2]]
        comp = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(rc_chunk)(*t) for t in tasks) for r in x])
        v = comp.R.dropna().to_numpy() if "R" in comp else np.array([])
        z = (r_clock - v.mean()) / v.std(ddof=1) if len(v) >= min(10, len(ws)) else np.nan
        p = float(norm.sf(z))
        crow.append(dict(cancer=c, readout=rd, cls="+".join(k.split(":")[1] for k in classes), form="pairwise", coupling=r_clock,
                         matched_mean=float(v.mean()), matched_sd=float(v.std(ddof=1)), matched_median=float(np.median(v)),
                         matched_q95=float(np.quantile(v, .95)), z=z, p_one_sided=p, share_matched_ge=float(np.mean(v >= r_clock)),
                         n_matched_ok=len(v), n_pairs=n_pairs))
        rows.append(dict(base, readout=rd, form="pairwise", n_classes_scored=2, z_pooled=z, p_pooled=p, **gate,
                         p401_pass=bool(gate.get("gate_passes", False) and np.isfinite(p) and p < .05), status="ok"))
    log(f"{c}: P401 (pairwise) done ({base['classes']})")


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
        if len(classes) < 2 or not set(B.BASE.genes) <= set(genes):
            for rd in READOUTS:
                rows.append(dict(base, readout=rd, p401_pass=False, status="not estimable: < 2 stromal classes with >= 5 patients"
                                 if len(classes) < 2 else "not estimable: clock genes missing"))
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
        if len(classes) == 2:
            two_class(c, base, classes, frames, ids, tab, Gu, ws, n_perm, args, rows, crow)
            pd.DataFrame(rows).to_csv(out / "p401.csv", index=False)
            continue
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
    import _shared as sh
    b = sh.tcga_builder()
    m = pd.read_csv(PURITY, sep="\t")
    m["case"] = m["array"].astype(str).str[:12]
    m = m[m["call status"] == "called"]
    p = m.drop_duplicates("case").set_index("case").purity
    rows = []
    for c in matched:
        smp = b.samples_of(MATCHED_TCGA[c])
        t = smp[smp.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
        v = p.reindex(t.patient).dropna()
        rows.append(dict(cancer=c, tcga=MATCHED_TCGA[c], primaries=len(t), with_purity=len(v), mean_purity=float(v.mean()),
                         median_purity=float(v.median()), sd_purity=float(v.std())))
    return pd.DataFrame(rows)


# ---- score stage ---------------------------------------------------------------------------------------------------------

A387 = EXP / "A387_tumour_full_model/output/score"
A387_NAME = {"COAD": "colorectal", "LUAD": "lung", "lung_mixed": "lung", "LUSC": "lung", "LIHC": "liver", "PRAD": "prostate"}
ALL_TCGA = {**ORIGINAL, **ADDED}
N_BOOT, N_PERM_FIT = 500, 1000


def fit_one(tab, classes, rd, targets, ref, n_boot, n_perm, readout):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    ds = cs.Dataset(tab, anchors=classes, waveform=READOUTS[rd], mixed_strata=True)
    kw = dict(seed=SEED, n_perm=n_perm, gate=False, reference=ref, attenuation="df", level=0.90, anchor_readout=readout, targets=targets)
    try:
        t = ds.fit(n_boot=n_boot, **kw).table
        t["interval_fallback"] = False
    except ValueError as e:  # a bootstrap replicate without an estimate for the reference: point estimates only, no interval
        t = ds.fit(n_boot=0, **kw).table
        t["interval_fallback"] = f"no intervals ({str(e)[:80]})"
    return t


def mal_chunk(Gu, classes, ws, rd):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in ws:
        try:
            t = cs.Dataset(Gu, anchors=classes, waveform=w, mixed_strata=True).fit(n_boot=0, gate=False, targets=[MALT]).table
            out.append(dict(set=i, readout=rd, coupling=float(t.loc[MALT, "coupling"])))
        except Exception as e:  # noqa: BLE001
            out.append(dict(set=i, readout=rd, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def spearman_exact(x, y):
    """Spearman rho and the exact one-sided permutation p (all n! orderings of y; ties as average ranks)."""
    from scipy.stats import rankdata
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan
    rx, ry = rankdata(x), rankdata(y)
    cx, cy = rx - rx.mean(), ry - ry.mean()
    den = np.sqrt((cx ** 2).sum() * (cy ** 2).sum())
    rho = float((cx * cy).sum() / den)
    perms = np.array(list(itertools.permutations(range(n))))
    r = (cx[None, :] * cy[perms]).sum(1) / den
    return rho, float(np.mean(r >= rho - 1e-12))


def stage_score(args):
    out = OUT / ("score_smoke" if args.smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    pw = OUT / "power"
    t0 = time.time()
    n_boot, n_perm, n_set = (20, 50, 6) if args.smoke else (N_BOOT, N_PERM_FIT, 200)
    frames, donors, _ = pool_sources(args, out)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    counts = pd.read_csv(pw / "counts.csv")
    if not args.smoke:
        now = donors.groupby("cancer").donor.nunique()
        assert (now.reindex(counts.cancer).to_numpy() == counts.patients.to_numpy()).all(), "pool differs from the power stage"
        log("pool reproduces the power-stage counts")
    if args.smoke:  # code check: malignant rows' donor labels permuted within stratum
        rng = np.random.default_rng(1)
        for f in frames.values():
            for st, g in f[f.compartment == MALT].groupby("stratum"):
                f.loc[g.index, "donor"] = g.donor.to_numpy()[rng.permutation(len(g))]
    p401 = pd.read_csv(pw / "p401.csv")
    p401c = pd.read_csv(pw / "p401_classes.csv")
    matched = pd.read_csv(pw / "matched.csv")
    purity = pd.read_csv(pw / "purity.csv").set_index("cancer")
    a387p = pd.read_csv(A387 / "part2_tests.csv")
    a387c = pd.read_csv(A387 / "per_cancer_tests.csv")
    from clocksync.paired import pseudo_waveforms
    cancers = [c for c in counts.cancer if counts.set_index("cancer").patients[c] >= 5]
    if args.smoke:
        cancers = [c for c in cancers if c in set(donors.cancer)]
    rows, comp_all = [], []
    for c in cancers:
        ids = set(donors[donors.cancer == c].donor)
        parts = [f[f.donor.isin(ids) & f.compartment.isin([MALT] + ANCH)] for f in frames.values()]
        parts = [p for p in parts if len(p)]
        Gc = pd.concat(parts, join="inner", ignore_index=True)
        genes = [g for g in Gc.columns if g not in META]
        classes = [a for a in ANCH if Gc[(Gc.compartment == a)].donor.nunique() >= 5]
        base = dict(cancer=c, patients=len(ids), cohorts=donors[donors.cancer == c].source.nunique(),
                    classes="+".join(k.split(":")[1] for k in classes))
        if len(classes) < 2 or not set(BASE.genes) <= set(genes):
            for rd in READOUTS:
                rows.append(dict(base, readout=rd, status="no host: < 2 stromal classes with >= 5 patients"))
            log(f"{c}: no host")
            continue
        Gc = Gc[Gc.compartment.isin([MALT] + classes)].reset_index(drop=True)
        tab = Gc[META + list(BASE.genes)]
        npat = {k: Gc[Gc.compartment == k].donor.nunique() for k in classes}
        ref = "tumour:fibroblast" if "tumour:fibroblast" in classes else max(classes, key=lambda k: npat[k])
        mode = "loo" if len(classes) >= 3 else "full"
        # composition from the single-cell pseudobulks of this cancer (depth = UMI)
        dep = Gc.pivot_table(index="donor", columns="compartment", values="depth", aggfunc="sum")
        s_sc = (dep[MALT] / dep[[MALT] + classes].sum(axis=1)).mean()
        w = {k: dep[k].sum() / sum(dep[j].sum() for j in classes) for k in classes}
        # matched sets (same draw as P401: the cancer's pooled anchors, genes in every cohort)
        anc = Gc[Gc.compartment.isin(classes)]
        frac = anc[genes].div(anc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(anc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        ws = list(enumerate(pseudo_waveforms(lm, n_set, np.random.default_rng([SEED, 0]), BASE, "matched")))
        used = sorted({g for _, wv in ws for g in wv.genes} | set(BASE.genes))
        Gu = Gc[META + used]
        tasks = [("fit_c", rd, dict(targets=[MALT], ref=None, n_boot=n_boot, n_perm=n_perm, readout="loo" if mode == "loo" else "full"))
                 for rd in READOUTS]
        tasks += [("fit_ref", rd, dict(targets=[MALT] + [k for k in classes], ref=ref, n_boot=n_boot, n_perm=0, readout=mode))
                  for rd in READOUTS]
        res = Parallel(n_jobs=min(args.jobs, len(tasks)))(
            delayed(fit_one)(tab, classes, rd, kw["targets"], kw["ref"], kw["n_boot"], kw["n_perm"], kw["readout"]) for _, rd, kw in tasks)
        fits = {(t[0], t[1]): r for t, r in zip(tasks, res)}
        mt = []
        for rd in READOUTS:
            w_rd = [(i, wv if rd == "14" else C.restrict(wv, BASE, CORE)) for i, wv in ws]
            for k in range(args.jobs * 2):
                if w_rd[k::args.jobs * 2]:
                    mt.append((Gu, classes, w_rd[k::args.jobs * 2], rd))
        comp = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(mal_chunk)(*t) for t in mt) for r in x])
        comp.assign(cancer=c).to_csv(out / f"matched_{c}.csv", index=False)
        for rd in READOUTS:
            fc, fr = fits[("fit_c", rd)].loc[MALT], fits[("fit_ref", rd)]
            m = comp[comp.readout == rd].coupling.dropna().to_numpy()
            cm = float(fc.coupling)
            z = (cm - m.mean()) / m.std(ddof=1)
            e = max(cm - float(np.median(m)), 0.0)
            amp = {k: float(np.exp(fr.loc[k, "log_amp_rel"])) for k in classes}
            r = dict(base, readout=rd, status="ok", interval_fallback=str(fc.get("interval_fallback", False)) + "|" + str(fr.interval_fallback.iloc[0]),
                     amplitude_reference=ref.split(":")[1], anchor_readout_amp=mode,
                     c_mal=cm, c_lo90=fc.coupling_lo, c_hi90=fc.coupling_hi, p_c=fc.p_coupling, boundary_type=fc.boundary_type,
                     host_reliability=fc.host_reliability, n_mal=int(fc.n), matched_mean=float(m.mean()), matched_sd=float(m.std(ddof=1)),
                     matched_median=float(np.median(m)), matched_q95=float(np.quantile(m, .95)), z_mal=z, e_mal=e,
                     share_matched_ge=float(np.mean(m >= cm)), n_matched_ok=len(m),
                     A_mal=float(np.exp(fr.loc[MALT, "log_amp_rel"])), A_mal_lo=float(np.exp(fr.loc[MALT, "log_amp_lo"])),
                     A_mal_hi=float(np.exp(fr.loc[MALT, "log_amp_hi"])), offset_mal_h=fr.loc[MALT, "offset_h"],
                     offset_mal_lo=fr.loc[MALT, "offset_lo"], offset_mal_hi=fr.loc[MALT, "offset_hi"],
                     **{f"A_{k.split(':')[1]}": amp[k] for k in classes},
                     **{f"offset_{k.split(':')[1]}_h": float(fr.loc[k, "offset_h"]) for k in classes},
                     sc_umi_share_mal=float(s_sc), **{f"w_{k.split(':')[1]}": float(w[k]) for k in classes})
            # composition and c_str
            if c in purity.index:
                s_mal, s_src = float(purity.mean_purity[c]), "TCGA ABSOLUTE purity"
            else:
                s_mal, s_src = float(s_sc), "single-cell malignant UMI share (non-TCGA cancer)"
            nm = A387_NAME.get(c)

            def cstr(kind):
                if nm is not None:
                    q = a387c[(a387c.readout.astype(str) == rd) & (a387c.kind == kind) & (a387c.cancer == nm)]
                    if len(q):
                        return float(q.coupling.iloc[0]), f"A387 per cancer ({nm})"
                q = a387p[(a387p.readout.astype(str) == rd) & (a387p.kind == kind)]
                return float(q.coupling.iloc[0]), "A387 pooled"
            cs_ = {}
            for k in classes:
                kind = "endothelial" if k.endswith("endothelial") else "fibroblast"
                cs_[k] = cstr(kind)
            c_chain = float(np.mean([v[0] for v in cs_.values()]))
            sk = {k: (1 - s_mal) * w[k] for k in classes}
            num = s_mal * r["A_mal"] * np.sqrt(e * c_chain) + sum(sk[k] * amp[k] * np.sqrt(cs_[k][0]) for k in classes)
            den = np.sqrt((s_mal * r["A_mal"]) ** 2 + sum((sk[k] * amp[k]) ** 2 for k in classes))
            r.update(s_mal=s_mal, s_mal_source=s_src, c_str_chain=c_chain, P=float(num / den),
                     **{f"s_{k.split(':')[1]}": sk[k] for k in classes},
                     **{f"c_str_{k.split(':')[1]}": cs_[k][0] for k in classes},
                     **{f"c_str_{k.split(':')[1]}_source": cs_[k][1] + (" (mural uses the fibroblast value)" if k.endswith("mural") else "")
                        for k in classes})
            rows.append(r)
        log(f"{c}: scored (reference {ref.split(':')[1]}, {mode})")
        pd.DataFrame(rows).to_csv(out / "per_cancer.csv", index=False)
    R = pd.DataFrame(rows)
    # bulk and A83 values
    a138 = pd.read_csv(EXP / "A138_cancer_estimability/output/score/per_unit.csv").query("hypothesis == 'H138.4'").set_index("unit").statistic
    added = pd.read_csv(pw / "bulk_pairing_added.csv").query("kind == 'added'").set_index("project").statistic
    bulk = pd.concat([a138, added])
    a83 = pd.concat([pd.read_csv(PAPER / "results/figure_5/output/cancer/a83_all_tumours.csv").set_index("project").share,
                     pd.read_csv(pw / "a83_added.csv").query("kind == 'added'").set_index("project").share])
    R["tcga"] = R.cancer.map(ALL_TCGA)
    R["bulk_pairing"] = R.tcga.map(bulk)
    R["a83_share"] = R.tcga.map(a83)
    pf = p401.set_index(["cancer", "readout"]).p401_pass
    R["p401_pass"] = [bool(pf.get((c, str(rd)), False)) for c, rd in zip(R.cancer, R.readout)]
    R["matched"] = R.cancer.isin(set(matched[matched.matched].cancer))
    R.to_csv(out / "per_cancer.csv", index=False)
    # D401 table: stroma couplings within the tumour from the power stage
    pc = p401c.assign(readout=p401c.readout.astype(str)).pivot_table(index=["cancer", "readout"], columns="cls", values="coupling")
    pc.columns = [f"tumour_stroma_coupling_{x}" for x in pc.columns]
    D = R.merge(pc.reset_index().assign(readout=lambda d: d.readout.astype(str)), on=["cancer", "readout"], how="left")
    D.to_csv(out / "d401.csv", index=False)
    # H401.1 / H401.2 and descriptive variants
    rowsH, loo = [], []
    crit = pd.read_csv(pw / "critical_rho.csv").set_index("n").critical_rho

    def sets(rd):
        s14 = R[(R.readout == "14") & R.p401_pass & R.matched & R.P.notna()].cancer.tolist()
        sc = R[(R.readout == "core") & R.p401_pass & R.matched & R.P.notna()].cancer.tolist()
        return s14 if rd == "14" else sc, s14
    for rd in READOUTS:
        own, s14 = sets(rd)
        for label, cs_set in (("passing P401 on this readout", own), ("the 14-gene passing set", s14)):
            d = R[(R.readout == rd) & R.cancer.isin(cs_set)]
            for hyp, xv, yv in (("H401.1 P vs bulk pairing", "P", "bulk_pairing"), ("H401.2 e_mal vs A83 share", "e_mal", "a83_share"),
                                ("descriptive c_mal vs bulk pairing", "c_mal", "bulk_pairing"),
                                ("descriptive e_mal vs bulk pairing", "e_mal", "bulk_pairing"),
                                ("descriptive z_mal vs bulk pairing", "z_mal", "bulk_pairing")):
                rho, p = spearman_exact(d[xv], d[yv])
                n = len(d)
                rowsH.append(dict(readout=rd, set=label, hypothesis=hyp, n=n, cancers="+".join(d.cancer), rho=rho, p_exact_one_sided=p,
                                  critical_rho=float(crit.get(n, np.nan)),
                                  decides=bool(rd == "14" and label.startswith("passing") and hyp.startswith("H401")),
                                  passes=bool(p < 0.05) if np.isfinite(p) else False))
        d14 = R[(R.readout == rd) & R.cancer.isin(s14)]
        if rd == "14" or True:
            for omit in d14.cancer:
                d = d14[d14.cancer != omit]
                for hyp, xv, yv in (("H401.1 P vs bulk pairing", "P", "bulk_pairing"), ("H401.2 e_mal vs A83 share", "e_mal", "a83_share")):
                    rho, p = spearman_exact(d[xv], d[yv])
                    loo.append(dict(readout=rd, omitted=omit, hypothesis=hyp, n=len(d), rho=rho, p_exact_one_sided=p))
    H = pd.DataFrame(rowsH)
    H.to_csv(out / "h401.csv", index=False)
    pd.DataFrame(loo).to_csv(out / "h401_leave_one_cancer_out.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 60)
    log("H401\n" + H.round(4).to_string(index=False))
    log("D401\n" + R[["cancer", "readout", "patients", "c_mal", "c_lo90", "c_hi90", "z_mal", "e_mal", "A_mal", "P", "bulk_pairing", "a83_share", "s_mal",
                      "p401_pass"]].round(3).to_string(index=False))
    log("leave-one-out\n" + pd.DataFrame(loo).round(4).to_string(index=False))
    log(f"score done {(time.time() - t0) / 60:.1f} min")


# ---- main ---------------------------------------------------------------------------------------------------------------

def pool_sources(args, out):
    """Sources, deduplication and genome-wide frames of the power stage (same code path; writes excluded / dropped / dedup)."""
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
    if (BASSEZ.exists() and not args.smoke) or (BASSEZ.exists() and args.smoke):
        try:
            G, d = build_bassez()
            frames["bassez2021"] = G
            dons.append(d.assign(priority=1))
            log(f"bassez2021 (Pre): {d.donor.nunique()} patients")
        except Exception as e:  # noqa: BLE001
            dropped.append(dict(source="bassez2021", stage="build", error=f"{type(e).__name__}: {e}"[:300]))
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
    return frames, donors, f388


def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    added = list(ADDED)[:1] if args.smoke else list(ADDED)
    bulk = bulk_pairing(out, added)
    a83_share(out, added)
    for _, r in bulk[bulk.kind == "added"].iterrows():
        if r.get("gated") is True and r.pairs >= 10:
            MATCHED_TCGA[r.project] = ADDED[r.project]
    log(f"matched TCGA set: {sorted(MATCHED_TCGA)}")
    frames, donors, f388 = pool_sources(args, out)
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
    mt = pd.DataFrame([dict(cancer=c, tcga=MATCHED_TCGA[c], patients=int(counts.set_index("cancer").patients[c]), matched=c in matched,
                            added_project=c in ADDED) for c in MATCHED_TCGA if c in set(counts.cancer)])
    mt["n_matched"] = nm
    mt["critical_rho_one_sided_05"] = critical_rho(nm) if 3 <= nm <= 14 else np.nan
    mt.to_csv(out / "matched.csv", index=False)
    pd.DataFrame([dict(n=n, critical_rho=critical_rho(n)) for n in range(5, 15)]).to_csv(out / "critical_rho.csv", index=False)
    log("counts\n" + counts.to_string(index=False))
    # P401 for every cancer with >= 5 patients
    cancers = [c for c in counts.cancer if counts.set_index("cancer").patients[c] >= 5]
    P = stage_p401(cancers, frames, donors, args, out)
    purity_table(matched).to_csv(out / "purity.csv", index=False)
    dec = []
    for rd in READOUTS:
        p = P[(P.readout == rd) & P.cancer.isin(matched)]
        ps = int(p.p401_pass.sum())
        dec.append(dict(readout=rd, role="primary (decides)" if rd == "14" else "beside (descriptive unless >= 7)",
                        matched_cancers=nm, matched_passing_p401=ps, needed=7, decidable=bool(ps >= 7),
                        passing="+".join(p[p.p401_pass].cancer), failing="+".join(p[~p.p401_pass].cancer)))
    pd.DataFrame(dec).to_csv(out / "decidability.csv", index=False)
    pd.set_option("display.width", 250)
    for n in ("matched", "p401", "p401_classes", "purity", "decidability"):
        log(n + "\n" + pd.read_csv(out / f"{n}.csv").round(4).to_string(index=False))
    log(f"power done {(time.time() - t0) / 60:.1f} min")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    stage_score(args) if args.stage == "score" else stage_power(args)


if __name__ == "__main__":
    main()
