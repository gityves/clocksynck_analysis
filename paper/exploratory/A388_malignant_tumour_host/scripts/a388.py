"""A388: malignant cells on their own tumour's stromal host, pooled over cancer types. Power stage.

Registered in paper/docs/registered_experiments.md (A388). Tables, labels and rules of the feasibility inventory
(`feasibility.py`, `breast_atlas.py`): A149 rules (>= 50 cells, sample median UMI >= 1,500, presence rule), malignant =
copy-number call where >= 5 patients are called, otherwise tumour epithelium; Wu / Pal author labels (Pal stroma by the
`cancer_common` marker rule). One `clocksync.Dataset` over all cohorts (donor and stratum prefixed by cohort,
`mixed_strata`), anchors = `tumour:{fibroblast, endothelial, mural}`, target = `tumour:malignant` (the malignant
compartment of every cohort under one name). Unit = patient with a malignant pseudobulk and >= 2 anchors.

BLINDNESS (power stage): no statistic of the malignant target on the real data is computed, stored or printed. Used:
counts; `Dataset.subset(compartments=anchors).gate` (anchors only, also per class); the blind simulator
`clocksync.simulate.like(ds, target, loading=z)` (the target is never regressed on the host); companion fits of
pseudo-waveforms (comparators).

  1 pooled table, deduplication (`dedup.csv`), labels new vs H138.1 / cancer / A83 class          -> counts.csv
  2 gates on the anchors (G388): pooled, new, keeping-class donors, at-null-class donors           -> gates.csv
  3 H388.1 (new patients) and H388.1b (all) power: simulated at |z| 0.5 / 1 and a grid of c        -> power_h1.csv
  4 H388.2 power: keeping vs at-null by simulating each class at its real design. The power proxy is the difference of
    the two classes' simulated couplings (the common-spread contrast d* is not simulated); null = both at c 0.02
    (c = 0 has no finite noise), power = share of differences above the null q95                    -> power_h2.csv
  5 companions on the real data: genome-wide pseudobulks, 200 matched + 200 narrow responsive sets, coupling of
    `tumour:malignant` for the new patients and all patients                                        -> companions.csv,
                                                                                       companion_summary.csv, dropped.csv
  6 decidability (MDE <= 0.25)                                                                   -> decidability.csv

New vs H138.1 is a PATIENT label: a patient is old if its cohort is an H138.1 cohort or it is the same patient as one
of an H138.1 cohort (deduplication), new otherwise.

    python a388.py --stage power [--jobs 4] [--genome-jobs 2] [--smoke]       (long: paper/common/run.sh)
    python a388.py --stage breast      (builds only the cached Wu / Pal genome-wide tables)
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
for p in ("A227_tme_compartment_pairing/scripts", "A149_paired_tumour_normal/scripts", "A151_fibroblast_pairing/scripts",
          "A226_stroma_pairing_confirm/scripts", "A336_two_readouts/scripts", "A170_sc_breast_er/scripts"):
    sys.path.insert(0, str(EXP / p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import a336_common as C  # noqa: E402
import clocksync as cs  # noqa: E402
import compute_a226 as K  # noqa: E402
import compute_a227 as M  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

K.install()

OUT = HERE.parents[1] / "output"
CACHE = OUT / "cache"
FEAS = OUT / "feasibility" / "inventory.csv"
BASE = C.WF
CORE = C.core(BASE)
READOUTS = {"14": BASE, "core": CORE}
ANCHORS = ["tumour:fibroblast", "tumour:endothelial", "tumour:mural"]
MALT = "tumour:malignant"
SEED = 388
N_SET, N_SIM, N_PERM = 200, 200, 1000
C_GRID = (0.02, 0.1, 0.2, 0.3, 0.4, 0.5)
C_KEEP = (0.2, 0.3, 0.4)
Z_GRID = (0.5, 1.0)
MAX_MDE = 0.25
KEEPING = ("breast", "lung", "liver", "prostate", "thyroid", "kidney")
AT_NULL = ("colorectal", "ovarian", "bladder", "endometrial")
BREAST = {"GSE176078": "Wu", "GSE161529": "Pal"}
MIN_CELLS, MIN_UMI = 50, 1500
META = ["donor", "compartment", "stratum", "n_cells", "depth"]  # table meta columns


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def a83(cancer):
    return "keeping" if cancer in KEEPING else "at_null" if cancer in AT_NULL else "other"


# ---- Wu and Pal genome-wide pseudobulks (author labels; cached) ------------------------------------------------------

def _genome_frame(rows, cohort):
    """rows: [(patient, compartment, n_cells, median_umi, Series of counts)] -> table in the genome_pair format."""
    X = pd.DataFrame({i: r[4] for i, r in enumerate(rows)}).T.fillna(0.0)
    X.columns = [canonical(str(g)) for g in X.columns]
    X = X.T.groupby(level=0).sum().T.reset_index(drop=True)
    meta = pd.DataFrame({"donor": [r[0] for r in rows], "compartment": [r[1] for r in rows], "stratum": cohort,
                         "n_cells": [r[2] for r in rows], "depth": X.sum(axis=1).to_numpy(), "median_umi": [r[3] for r in rows]})
    return pd.concat([meta, X.astype(np.float32)], axis=1)


def build_wu():
    import prep_a170 as P
    meta = pd.read_csv(P.WU / "metadata.csv", index_col=0)
    bc = pd.read_csv(P.WU / "count_matrix_barcodes.tsv", header=None)[0].to_numpy()
    genes = pd.read_csv(P.WU / "count_matrix_genes.tsv", header=None)[0].to_numpy()
    meta = meta.reindex(bc)
    assert meta["orig.ident"].notna().all()
    pat = meta["orig.ident"].to_numpy()
    cmap = {"Cancer Epithelial": MALT, "CAFs": "tumour:fibroblast", "Endothelial": "tumour:endothelial", "PVL": "tumour:mural"}
    cls = meta.celltype_major.map(cmap).to_numpy()
    units = [(p, c) for p in sorted(set(pat)) for c in cmap.values()]
    ix = {u: j for j, u in enumerate(units)}
    u_of = np.array([ix.get((p, c), -1) for p, c in zip(pat, cls)])
    U, G = len(units), len(genes)
    acc, n = np.zeros(G * U), 0
    for ch in pd.read_csv(P.WU / "count_matrix_sparse.mtx", sep=" ", skiprows=2, header=None, chunksize=20_000_000,
                          dtype={0: np.int64, 1: np.int64, 2: np.float64}):
        g, c, v = ch[0].to_numpy() - 1, ch[1].to_numpy() - 1, ch[2].to_numpy()
        k = u_of[c] >= 0
        acc += np.bincount(g[k] * U + u_of[c][k], weights=v[k], minlength=G * U)
        n += len(ch)
        log(f"Wu: {n:,} entries")
    acc = acc.reshape(G, U)
    ncell = pd.Series(u_of[u_of >= 0]).value_counts()
    med = meta.groupby("orig.ident").nCount_RNA.median()
    rows = [(p, c, int(ncell.get(j, 0)), float(med[p]), pd.Series(acc[:, j], index=genes)) for j, (p, c) in enumerate(units)]
    return _genome_frame(rows, "Wu")


def build_pal():
    import zipfile

    import prep_a170 as P
    import rds_parse  # noqa: F401  (prep_a170 imports it; path check)
    from cancer_common.compartments import marker_compartments
    import scipy.sparse as sp
    z = zipfile.ZipFile(P.PAL_ZIP)
    rows = []
    for obj, (er, tum) in P.PAL_OBJ.items():
        md = P.R.read_rds(P.D / f"GSE161529/figshare17058077/SeuratObject_{obj}.rds")["meta.data"]
        vals = dict(zip(list(md["names"].data), md.data))
        cells = pd.DataFrame({"group": P.fac(vals["group"]), "cluster": P.fac(vals["seurat_clusters"]).astype(int)},
                             index=list(md["row.names"].data))
        for g, d in cells.groupby("group"):
            bc, sym, gi, ci, v = P.read10x(z, g.replace("_", "-"))
            pos = pd.Series(np.arange(len(bc)), index=[f"{g}_{b}" for b in bc])[d.index].to_numpy()
            X = sp.csr_matrix((v, (ci, gi)), shape=(len(bc), len(sym)))[pos]
            umi = np.asarray(X.sum(1)).ravel()
            comp = np.asarray(marker_compartments(X, sym, umi)).astype(object)
            tumc = d.cluster.isin(tum).to_numpy()
            comp[tumc] = "malignant"
            for k, name in (("malignant", MALT), ("fibroblast", "tumour:fibroblast"), ("endothelial", "tumour:endothelial"),
                            ("mural", "tumour:mural")):
                i = np.flatnonzero(comp == k)
                s = np.asarray(X[i].sum(0)).ravel() if len(i) else np.zeros(len(sym))
                rows.append((g, name, len(i), float(np.median(umi)), pd.Series(s, index=sym).groupby(level=0).sum()))
            log(f"Pal {obj} {g}: {len(d)} cells")
    return _genome_frame(rows, "Pal")


def breast_genome(study):
    f = CACHE / f"breast_{BREAST[study]}_genome.parquet"
    if f.exists():
        return pd.read_parquet(f)
    CACHE.mkdir(parents=True, exist_ok=True)
    G = build_wu() if BREAST[study] == "Wu" else build_pal()
    G.to_parquet(f)
    return G


def qualify(G):
    """Breast atlas rows passing the A149 rules (>= 50 cells, sample median UMI >= 1,500)."""
    return G[(G.n_cells >= MIN_CELLS) & (G.median_umi >= MIN_UMI)].drop(columns="median_umi")


# ---- cohort tables ---------------------------------------------------------------------------------------------------

def cohort_table(study, mal_label=None):
    """Pooled-format table of one cohort (clock genes only; target renamed `tumour:malignant`; patients with a malignant
    pseudobulk and >= 2 anchors) with the list of all its patients (for deduplication). Returns (table, all patients,
    malignant label) or an error string."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        if study in BREAST:
            pb, mal = qualify(breast_genome(study)), MALT
            allp = sorted(breast_genome(study).donor.unique())
        else:
            if study in K.COHORTS:
                pb = K.table(study)
                mal = L.MAL if L.MAL in set(pb.compartment) else L.EPI
            else:
                pb, mal = L.table(study)
            allp = sorted(pb.donor.astype(str).unique())
        pb = pb[pb.compartment.isin([mal] + ANCHORS)].copy()
        pb["donor"] = pb.donor.astype(str)
        pb["compartment"] = pb.compartment.where(pb.compartment != mal, MALT)
        n_anchor = pb[pb.compartment.isin(ANCHORS)].groupby("donor").compartment.nunique()
        has_mal = set(pb[pb.compartment == MALT].donor)
        keep = sorted(d for d in has_mal if n_anchor.get(d, 0) >= 2)
        pb = pb[pb.donor.isin(keep)]
        pb = pb[(pb.compartment == MALT) | pb.compartment.isin(ANCHORS)]
        pb["stratum"] = pb.stratum.astype(str) if "stratum" in pb else ""
        pb = pb[META + list(BASE.genes)].reset_index(drop=True)
        return pb, allp, mal
    except Exception as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"[:300]


def study_list(smoke):
    inv = pd.read_csv(FEAS)
    inv = inv[(inv.n_ge2_anchors > 0)]
    inv = inv.set_index("study")
    studies = list(inv.index)
    if smoke:
        studies = [s for s in ("SRP416180", "GSE222315", "GSE315534", "GSE281736", "GSE176078") if s in studies]
    return studies, inv


# ---- deduplication ---------------------------------------------------------------------------------------------------

def dedup(tabs, allp, inv, h138):
    """Patients shared between an scBaseCount cohort (SRP / ERP) and the GEO cohort of the same study: the inventory's GEO
    series (`gse`) of the scBaseCount cohort is matched to the cohorts in the list; patient labels (case-insensitive)
    identify the shared patients. A shared patient is kept once, from the GEO / author copy if that copy qualifies (>= 2
    anchors and a malignant pseudobulk), else from the scBaseCount copy. Returns (set of removed (study, patient),
    dedup records, set of retained patients that are also an H138.1 patient)."""
    iv = pd.read_csv(L.ST.INVENTORY, low_memory=False)
    qual = {s: set(t.donor) for s, t in tabs.items()}
    removed, rec, h138_dup = set(), [], set()
    norm = lambda x: str(x).strip().lower()  # noqa: E731
    for s in tabs:
        if s.startswith("GSE") or s in BREAST:
            continue
        gses = sorted(set(iv[iv.study == s].gse.dropna().astype(str)))
        geo = [g for g in gses if g in tabs]
        if gses and not geo:
            extra = (f"; the Wu ({'GSE176078'}) and Pal ({'GSE161529'}) atlases have other series, so no overlap with them"
                     if inv.cancer[s] == "breast" else "")
            rec.append(dict(kept=s, removed="", patient="", status="no overlap",
                            note=f"its GEO series ({'; '.join(gses)}) has no qualifying patient in the pooled list, so there is "
                                 f"no second copy{extra}"))
            continue
        if not geo:
            same = [o for o in tabs if o != s and inv.cancer[o] == inv.cancer[s]]
            coll = sum(len({norm(d) for d in allp[s]} & {norm(d) for d in allp[o]}) > 0 for o in same if o in allp)
            rec.append(dict(kept=s, removed="", patient="", status="not checkable",
                            note=f"no GEO series in the sample inventory: overlap with the "
                                 f"{len(same)} same-cancer cohorts can only be tested by patient label ({coll} cohorts share "
                                 f"a label, uninformative for generic ids); no overlap established"))
            continue
        for g in geo:
            if g not in tabs:
                rec.append(dict(kept=s, removed="", patient="", status="no overlap", note=f"GEO series {g} is not in the pooled list"))
                continue
            la = {norm(d): d for d in allp[g]}
            lb = {norm(d): d for d in allp[s]}
            shared = sorted(set(la) & set(lb))
            rec.append(dict(kept=g, removed=s, patient="", status="pair", note=f"accession {g} of {s}; {len(shared)} shared patient labels "
                                                                              f"of {len(la)} (GEO) / {len(lb)} (scBaseCount)"))
            for p in shared:
                gq, sq = la[p] in qual[g], lb[p] in qual[s]
                if gq and sq:
                    removed.add((s, lb[p]))
                    rec.append(dict(kept=g, removed=s, patient=la[p], status="removed scBaseCount copy",
                                    note="both copies qualify; GEO copy kept"))
                    if inv.in_H138_1.get(s, False):
                        h138_dup.add((g, la[p]))
                elif sq and not gq:
                    removed.add((g, la[p])) if la[p] in qual[g] else None
                    rec.append(dict(kept=s, removed=g, patient=lb[p], status="kept scBaseCount copy",
                                    note="GEO copy does not qualify; scBaseCount copy kept"))
                else:
                    rec.append(dict(kept="", removed="", patient=la[p], status="not qualifying", note="shared but qualifies in neither"))
    return removed, pd.DataFrame(rec), h138_dup


# ---- genome-wide pseudobulks of the anchors and the malignant compartment ---------------------------------------------

def genome_cohort(study, mal, pts, keep):
    """Genome-wide pseudobulks (donor x compartment) of the malignant and anchor compartments of `pts` of one cohort, with the
    labels of the table (copy-number calls for `tumour:malignant` cohorts, as `a149_lib.genome_epithelium`; stromal classes
    by the marker rule, as `compute_a227.genome_pair`); rows restricted to `keep` = {(patient, compartment)}."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    K.install()  # routes the A226 cohorts of `cancer_common.studies.units` (a side effect lost in loky workers)
    f = CACHE / f"genome_{study}.parquet"
    if f.exists():
        G = pd.read_parquet(f)
        if {(d, c) for d, c in zip(G.donor, G.compartment)} == keep:
            return study, G, ""
    try:
        if study in BREAST:
            G = breast_genome(study)
            G = G[[(d, c) in keep for d, c in zip(G.donor, G.compartment)]].drop(columns="median_umi")
        else:
            us, _ = L.ST.units(study)
            us = [u for u in L.ST.majority_chemistry(us) if u[1] in set(pts)]
            calls = None
            if mal == L.MAL:
                calls = L.new_calls(study) if study in L.NEW else L.ST.read_calls(study)
            sums, depth, strat, genes = {}, {}, {}, None
            for sample, patient, kind, load, chem, _ in us:
                if kind != "tumour":
                    continue
                X, _g, umi, comp, barcodes, sym = load()
                lab = L.B._labels(comp, kind, sample, barcodes, calls)
                if genes is None:
                    genes = pd.Index(sym)
                for c in [mal] + ANCHORS:
                    i = np.flatnonzero(lab == c)
                    if len(i):
                        k = (patient, c)
                        sums[k] = sums.get(k, 0) + np.asarray(X[i].sum(0)).ravel()
                        depth[k] = depth.get(k, 0) + float(umi[i].sum())
                        strat[k] = chem
            keys = sorted(sums)
            G = pd.DataFrame(np.vstack([sums[k] for k in keys]).astype(np.float32), columns=genes)
            G.columns = [canonical(str(g)) for g in G.columns]
            G = G.T.groupby(level=0).sum().T
            meta = pd.DataFrame({"donor": [k[0] for k in keys], "compartment": [k[1] for k in keys],
                                 "stratum": [strat[k] for k in keys], "depth": [depth[k] for k in keys]})
            G = pd.concat([meta, G.reset_index(drop=True)], axis=1)
            G["compartment"] = G.compartment.where(G.compartment != mal, MALT)
            G = G[[(d, c) in keep for d, c in zip(G.donor, G.compartment)]]
        G = G.reset_index(drop=True)
        if not len(G):
            return study, None, "no rows"
        miss = keep - {(d, c) for d, c in zip(G.donor, G.compartment)}
        if miss:
            return study, None, f"{len(miss)} table rows have no genome row"
        CACHE.mkdir(parents=True, exist_ok=True)
        G.to_parquet(f)
        return study, G, ""
    except Exception as e:  # noqa: BLE001
        return study, None, f"{type(e).__name__}: {e}"[:300]


# ---- datasets, gates, simulation ---------------------------------------------------------------------------------------

def dataset(table, wf, name=""):
    return cs.Dataset(table, anchors=ANCHORS, waveform=wf, mixed_strata=True, name=name)


def gate_row(ds, label, rd, n_perm, donors=None):
    sub = ds.subset(donors=donors, compartments=ANCHORS) if donors is not None else ds.subset(compartments=ANCHORS)
    g = sub.gate(n_perm=n_perm, seed=SEED)
    return dict(set=label, readout=rd, n_donors=g.n_donors, n_pairs=g.n_pairs, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                clock_share=g.clock_share, share_null_mean=g.share_null_mean, passes=bool(g.passes))


def sim_task(table, rd, targets_c, z, seed, n_sim):
    """Draws of the coupling of simulated targets at the real design (blind: `like(..., loading=z)`). targets_c:
    {target: [c, ...]}; returns {(target, c): array}."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    ds = dataset(table, READOUTS[rd])
    out = {}
    for ti, (t, cs_) in enumerate(targets_c.items()):
        sim = cs.simulate.like(ds, t, loading=z)
        rng = np.random.default_rng([*seed, ti])
        for c in cs_:
            out[(t, c)] = np.array([sim.draw(c, rng)["coupling"] for _ in range(n_sim)], float)
    return out


def mde(power_by_c, thr_c=C_GRID[0]):
    ok = [c for c, p in power_by_c.items() if c > thr_c and p >= 0.8]
    return min(ok) if ok else np.nan


def stage_sims(table, donors, args, out):
    """Power of H388.1 (new patients), H388.1b (all) and H388.2 (keeping vs at-null)."""
    n_sim = 20 if args.smoke else N_SIM
    ds_new = table[table.donor.isin(donors[donors.new].donor)]
    sets = {"H388.1": ds_new, "H388.1b": table}
    tasks = [(h, rd, z) for h in sets for rd in READOUTS for z in Z_GRID]
    res = Parallel(n_jobs=args.jobs, return_as="generator")(
        delayed(sim_task)(sets[h], rd, {MALT: list(C_GRID)}, z, (SEED, 1, i), n_sim) for i, (h, rd, z) in enumerate(tasks))
    rows = []
    for (h, rd, z), r in zip(tasks, res):
        t = sets[h]
        thr = float(np.nanquantile(r[(MALT, C_GRID[0])], .95))
        pw = {c: float(np.nanmean(r[(MALT, c)] > thr)) for c in C_GRID}
        m = mde(pw)
        pat = t[t.compartment == MALT].donor.nunique()
        for c in C_GRID:
            rows.append(dict(hypothesis=h, readout=rd, z_abs=z, c_true=c, c_median=float(np.nanmedian(r[(MALT, c)])),
                             null_threshold=thr, power=pw[c], mde_c=m, patients=pat,
                             cohorts=t.donor.str.split("|").str[0].nunique(), n_sim=n_sim))
        log(f"power {h} {rd} |z| {z}: MDE {m}")
    P1 = pd.DataFrame(rows)
    P1.to_csv(out / "power_h1.csv", index=False)
    # H388.2
    cls = donors.set_index("donor").a83
    t2 = table[table.donor.map(cls).isin(["keeping", "at_null"])].copy()
    lab = t2.donor.map(cls)
    t2["compartment"] = np.where(t2.compartment == MALT, MALT + "|" + lab, t2.compartment)
    tk, tn = MALT + "|keeping", MALT + "|at_null"
    tasks = [(rd, z) for rd in READOUTS for z in Z_GRID]
    res = Parallel(n_jobs=args.jobs, return_as="generator")(
        delayed(sim_task)(t2, rd, {tk: [C_GRID[0], *C_KEEP], tn: [C_GRID[0]]}, z, (SEED, 2, i), n_sim)
        for i, (rd, z) in enumerate(tasks))
    rows = []
    for (rd, z), r in zip(tasks, res):
        null = r[(tk, C_GRID[0])] - r[(tn, C_GRID[0])]
        thr = float(np.nanquantile(null, .95))
        pw = {c: float(np.nanmean((r[(tk, c)] - r[(tn, C_GRID[0])]) > thr)) for c in C_KEEP}
        m = mde(pw)
        for c in C_KEEP:
            d = r[(tk, c)] - r[(tn, C_GRID[0])]
            rows.append(dict(hypothesis="H388.2", readout=rd, z_abs=z, c_keep=c, c_null=C_GRID[0], diff_median=float(np.nanmedian(d)),
                             null_diff_q95=thr, null_diff_sd=float(np.nanstd(null)), power=pw[c], mde_c_keep=m,
                             mde_diff=round(m - C_GRID[0], 4) if np.isfinite(m) else np.nan,
                             patients_keeping=t2[t2.compartment == tk].donor.nunique(),
                             patients_at_null=t2[t2.compartment == tn].donor.nunique(),
                             cohorts_keeping=t2[t2.compartment == tk].donor.str.split("|").str[0].nunique(),
                             cohorts_at_null=t2[t2.compartment == tn].donor.str.split("|").str[0].nunique(), n_sim=n_sim))
        log(f"power H388.2 {rd} |z| {z}: MDE c_keep {m}")
    P2 = pd.DataFrame(rows)
    P2.to_csv(out / "power_h2.csv", index=False)
    return P1, P2


# ---- companions -------------------------------------------------------------------------------------------------------

def companion_chunk(G, ws, rd, name):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in ws:
        try:
            t = dataset(G, w).fit(n_boot=0, gate=False, targets=[MALT]).table
            out.append(dict(cohort_set=name, readout=rd, set=i, coupling=float(t.loc[MALT, "coupling"])))
        except Exception as e:  # noqa: BLE001
            out.append(dict(cohort_set=name, readout=rd, set=i, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def stage_companions(G, table, donors, args, out):
    from clocksync.paired import pseudo_waveforms
    n_set = 6 if args.smoke else N_SET
    genes = [c for c in G.columns if c not in META and c != "median_umi"]
    genes = [g for g in genes if G[g].notna().all()]  # genes measured in every cohort
    ref = G[G.compartment.isin(ANCHORS)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    sets = [(f, w) for code, f in enumerate(("matched", "responsive"))
            for w in pseudo_waveforms(lm, n_set, np.random.default_rng([SEED, code]), BASE, f)]
    used = sorted({g for _, w in sets for g in w.genes} | set(BASE.genes))
    Gu = G[META + used]
    log(f"companions: {len(sets)} sets, {len(used)} genes, {len(Gu)} rows, {len(genes)} genes in every cohort")
    new_donors = set(donors[donors.new].donor)
    jobs = []
    for name, Gs in (("new", Gu[Gu.donor.isin(new_donors)]), ("all", Gu)):
        for rd in READOUTS:
            ws = [(i, w if rd == "14" else C.restrict(w, BASE, CORE)) for i, (_, w) in enumerate(sets)]
            for k in range(args.jobs * 4):
                if ws[k::args.jobs * 4]:
                    jobs.append((Gs, ws[k::args.jobs * 4], rd, name))
    res = Parallel(n_jobs=args.jobs, return_as="generator")(delayed(companion_chunk)(*j) for j in jobs)
    rows, n = [], 0
    for r in res:
        rows += r
        n += 1
        if n % 8 == 0:
            log(f"companions: {n}/{len(jobs)} chunks")
    comp = pd.DataFrame(rows)
    comp["family"] = [sets[i][0] for i in comp.set]
    comp.to_csv(out / "companions.csv", index=False)
    if "error" in comp:
        log("companion fit errors:", int(comp.error.notna().sum()))
    s = comp.groupby(["cohort_set", "readout", "family"]).coupling.agg(
        n_ok=lambda x: int(x.notna().sum()), median=lambda x: float(np.nanmedian(x)), q95=lambda x: float(np.nanquantile(x, .95)),
        q99=lambda x: float(np.nanquantile(x, .99))).reset_index()
    mal_d = Gu[Gu.compartment == MALT].donor
    s["patients"] = [int(mal_d.isin(new_donors).sum()) if n == "new" else int(mal_d.nunique()) for n in s.cohort_set]
    s.to_csv(out / "companion_summary.csv", index=False)
    return s


# ---- main ---------------------------------------------------------------------------------------------------------------

def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    studies, inv = study_list(args.smoke)
    h138 = set(inv.index[inv.in_H138_1.astype(bool)])
    log(f"{len(studies)} cohorts with >= 1 patient (feasibility): {studies}")
    res = Parallel(n_jobs=args.genome_jobs)(delayed(cohort_table)(s) for s in studies)
    dropped = [dict(study=s, stage="table", error=r) for s, r in zip(studies, res) if isinstance(r, str)]
    tabs = {s: r[0] for s, r in zip(studies, res) if not isinstance(r, str) and len(r[0])}
    allp = {s: r[1] for s, r in zip(studies, res) if not isinstance(r, str)}
    mal = {s: r[2] for s, r in zip(studies, res) if not isinstance(r, str)}
    log(f"tables: {len(tabs)} cohorts, {sum(t.donor.nunique() for t in tabs.values())} patients before deduplication")
    removed, drec, h138_dup = dedup(tabs, allp, inv, h138)
    drec.to_csv(out / "dedup.csv", index=False)
    for s in tabs:
        tabs[s] = tabs[s][[(s, d) not in removed for d in tabs[s].donor]]
    tabs = {s: t for s, t in tabs.items() if len(t)}
    log(f"after deduplication: {sum(t.donor.nunique() for t in tabs.values())} patients ({len(removed)} removed)")
    # genome-wide tables (companions) first: a cohort that cannot be built is dropped from every test
    keeps = {s: set(zip(t.donor, t.compartment)) for s, t in tabs.items()}
    gres = Parallel(n_jobs=args.genome_jobs)(
        delayed(genome_cohort)(s, mal[s], sorted(t.donor.unique()), keeps[s]) for s, t in tabs.items())
    Gs = {s: g for s, g, e in gres if g is not None}
    dropped += [dict(study=s, stage="genome", error=e) for s, g, e in gres if g is None]
    log(f"genome tables: {len(Gs)} built, {len(dropped)} dropped in total")
    tabs = {s: t for s, t in tabs.items() if s in Gs}
    pd.DataFrame(dropped, columns=["study", "stage", "error"]).to_csv(out / "dropped.csv", index=False)
    # pooled table, donor labels
    T = []
    for s, t in tabs.items():
        t = t.copy()
        t["donor"] = s + "|" + t.donor
        t["stratum"] = s + "|" + t.stratum
        T.append(t)
    table = pd.concat(T, ignore_index=True)
    dl = []
    for s, t in tabs.items():
        for d in t.donor.unique():
            dl.append(dict(donor=f"{s}|{d}", cohort=s, patient=d, cancer=inv.cancer[s], a83=a83(inv.cancer[s]),
                           cohort_in_H138_1=s in h138, dup_of_H138_1=(s, d) in h138_dup))
    donors = pd.DataFrame(dl)
    donors["new"] = ~(donors.cohort_in_H138_1 | donors.dup_of_H138_1)
    # counts
    cnt = []
    for by in (["cohort", "cancer", "a83", "new"], ["cancer", "a83", "new"], ["a83", "new"], ["a83"], ["new"], []):
        g = donors.groupby(by) if by else [((), donors)]
        for k, d in g:
            k = k if isinstance(k, tuple) else (k,)
            cnt.append(dict(by="+".join(by) or "total", cohort=d.cohort.iloc[0] if "cohort" in by else "",
                            cancer=d.cancer.iloc[0] if "cancer" in by else "", a83=d.a83.iloc[0] if "a83" in by else "",
                            new=bool(d.new.iloc[0]) if "new" in by else "", patients=len(d), cohorts=d.cohort.nunique()))
    pd.DataFrame(cnt).to_csv(out / "counts.csv", index=False)
    # G388
    new_d = list(donors[donors.new].donor)
    gates = []
    for rd, wf in READOUTS.items():
        ds = dataset(table, wf, f"A388_{rd}")
        for label, dd in (("pooled", None), ("new", new_d), ("keeping", list(donors[donors.a83 == "keeping"].donor)),
                          ("at_null", list(donors[donors.a83 == "at_null"].donor))):
            gates.append(gate_row(ds, label, rd, 100 if args.smoke else N_PERM, dd))
            log(f"gate {label} {rd}: p_pairing {gates[-1]['p_pairing']:.4f}, p_acrophase {gates[-1]['p_acrophase']:.4f}")
    pd.DataFrame(gates).to_csv(out / "gates.csv", index=False)
    P1, P2 = stage_sims(table, donors, args, out)
    # companions
    G = pd.concat([Gs[s].assign(donor=s + "|" + Gs[s].donor.astype(str)) for s in tabs], ignore_index=True)
    G = G.drop(columns=[c for c in ("stratum", "n_cells") if c in G]).merge(
        table[["donor", "compartment", "stratum", "n_cells"]], on=["donor", "compartment"])
    S = stage_companions(G, table, donors, args, out)
    # decidability
    dec = []
    for (h, rd, z), g in P1.groupby(["hypothesis", "readout", "z_abs"]):
        m = g.mde_c.iloc[0]
        dec.append(dict(hypothesis=h, readout=rd, z_abs=z, mde=m, patients=g.patients.iloc[0], decidable=bool(m <= MAX_MDE)))
    for (rd, z), g in P2.groupby(["readout", "z_abs"]):
        m = g.mde_diff.iloc[0]
        dec.append(dict(hypothesis="H388.2", readout=rd, z_abs=z, mde=m, patients=int(g.patients_keeping.iloc[0] + g.patients_at_null.iloc[0]),
                        decidable=bool(m <= MAX_MDE)))
    D = pd.DataFrame(dec)
    D.to_csv(out / "decidability.csv", index=False)
    pd.set_option("display.width", 250)
    for name, t in (("counts", pd.read_csv(out / "counts.csv")), ("gates", pd.DataFrame(gates)), ("decidability", D),
                    ("companions", S)):
        log(f"{name}\n" + t.round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "breast"), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    if args.stage == "breast":
        for s in BREAST:
            breast_genome(s)
            log("built", s)
    else:
        stage_power(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
