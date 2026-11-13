"""Cell classes of the single-cell fetal adrenal cohorts (A244): per-cell marker counts, the per-cell rule, class
pseudobulks per fetus, foreign-RNA fractions and the design table.

Same cohorts, fetuses, samples, ages and cells as `A242/scripts/prep_adrenal.py` (its parsers and constants are
imported): kameneva, wang, sasaki, kildisiute, dong, han. Stages (each is run once, in this order):

  markers      per cell: total UMI count and the counts of the registered marker genes (cortex, medulla, stroma,
               blood). Nothing else is read for the rule: the 14 clock genes and pseudo-clock sets never enter it.
               Also the published cell annotation where one is on disk (Kildisiute, Dong). Writes
               `$SCCC_DATA/natexp/fetal_adrenal_open/a244_cells.parquet` (one row per cell).
  hist         marker-fraction histograms (`../output/marker_hist.csv`) from which the thresholds in `RULE` were set.
  classes      applies `RULE`; class sizes per fetus (`../output/class_sizes.csv`), agreement with the published
               annotation (`../output/agreement_*.csv`).
  pseudobulk   genome-wide sums per fetus and class (all genes: the clock genes enter only as counts), checked against
               the whole-adrenal pseudobulks of A242; `a244_class_pseudobulk.parquet` (every class and fetus) and
               `a244_class_samples.csv` (classes with >= MIN_CELLS cells).
  design       foreign-RNA fractions (`../output/foreign_fetus.csv`, `foreign_cohort.csv`), the design table
               (`../output/design_classes.csv`) and the median clock-gene counts of the class pseudobulks
               (`../output/clock_detection_classes.csv`).

NOT computed: any clock statistic (coherence, pairing, phase, companion). The median clock-gene count of a class
pseudobulk is a count, as in the detection table of A242.

    python -I prep_a244.py --stage markers|hist|classes|pseudobulk|design
"""
import argparse
import gzip
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "A242_fetal_adrenal_replication" / "scripts"))
import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import prep_adrenal as a242  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC, DEST, OUT = a242.SRC, a242.DEST, HERE.parents[1] / "output"
CELLS, PB, SHEET = DEST / "a244_cells.parquet", DEST / "a244_class_pseudobulk.parquet", DEST / "a244_class_samples.csv"
MARKERS = {"cortex": ["STAR", "CYP11A1", "CYP17A1", "CYP21A2", "CYP11B1", "NR5A1"],
           "medulla": ["CHGA", "CHGB", "TH", "DBH", "PHOX2B", "STMN2", "ELAVL4", "SOX10", "PLP1"],
           "stroma": ["COL1A1", "COL1A2", "COL3A1", "PDGFRB", "PECAM1", "CDH5"],
           "blood": ["HBG1", "HBG2", "HBA1", "HBB", "PTPRC"]}
ALL_MARKERS = [g for v in MARKERS.values() for g in v]
CLASSES = ("cortex", "medulla", "stroma")
MIN_CELLS = 200  # cells of a class and fetus for a pseudobulk (A190)
# Per-cell rule (thresholds set on the marker fractions alone: `--stage hist`; see the design addendum):
# a class is met when its marker count >= MIN_COUNT and its marker fraction (marker UMIs / total UMIs) >= THRESH.
RULE = dict(min_count=3, cortex=0.005, medulla=0.005, stroma=0.005, blood=0.1)
COHORTS = ["kameneva", "wang", "sasaki", "kildisiute", "dong", "han"]


def _canon():
    cache = {}
    return lambda names: [cache.setdefault(g, canonical(g)) for g in names]


def _wide(h, symbol_col=None):
    """Gene x cell text table -> (csc genes x cells float32, symbols, barcodes), float32 to save memory."""
    first = pd.read_csv(h, sep="\t", nrows=0).columns.tolist()
    h.seek(0) if hasattr(h, "seek") else None
    types = defaultdict(lambda: np.float32, {first[0]: str, **({symbol_col: str} if symbol_col else {})})
    t = pd.read_csv(h, sep="\t", index_col=0, dtype=types)
    sym = t[symbol_col].to_numpy(str) if symbol_col else t.index.to_numpy(str)
    t = t.drop(columns=[symbol_col]) if symbol_col else t
    return sp.csc_matrix(t.to_numpy(np.float32)), sym, t.columns.to_numpy(str)


def _tenx(h, min_umi=0):
    with h5py.File(h, "r") as f:
        m = f["matrix"]
        x = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
        names, bc = m["features/name"][:].astype(str), m["barcodes"][:].astype(str)
    keep = np.asarray(x.sum(0)).ravel() >= max(min_umi, 1)
    return x[:, keep].astype(np.float32), names, bc[keep]


def samples(cohorts=COHORTS):
    """Yield (cohort, fetus, sample, X genes x cells, gene names, barcodes): the cells of A242, sample by sample."""
    if "kameneva" in cohorts:
        with tarfile.open(SRC / "kameneva2021" / "raw" / "GSE147821_RAW.tar") as tar:
            for name in a242.KAMENEVA:
                x, n, b = _tenx(tar.extractfile(f"{name}.raw_feature_bc_matrix.h5"), a242.MIN_UMI)
                yield "kameneva", name.split("_", 1)[1], name, x, n, b
    if "dong" in cohorts:
        with tarfile.open(SRC / "dong2020" / "raw" / "GSE137804_RAW.tar") as tar:
            for name in a242.DONG:
                with gzip.open(tar.extractfile(f"{name}_gene_cell_exprs_table.xls.gz"), "rt") as h:
                    x, n, b = _wide(h, "Symbol")
                yield "dong", name.split("_")[1], name, x, n, b
    if "wang" in cohorts:
        for name, (fetus, _) in a242.WANG.items():
            x, n, b = _tenx(SRC / "wang2022" / "raw" / f"{name}_filtered_feature_bc_matrix.h5")
            yield "wang", fetus, name, x, n, b
    if "han" in cohorts:
        for name, age in a242.HAN.items():
            with gzip.open(SRC / "han2020" / "raw" / f"{name}_dge.txt.gz", "rt") as h:
                x, n, b = _wide(h)
            yield "han", f"{age}w", name, x, n, b
    if "sasaki" in cohorts:
        for name, age in a242.SASAKI.items():
            x, n, b = _tenx(SRC / f"{name}_filtered_feature_bc_matrix.h5")
            yield "sasaki", f"{age}wpf", name.split("/")[-1], x, n, b
    if "kildisiute" in cohorts:
        a = ad.read_h5ad(SRC / "kildisiute2021" / "raw" / "adr_all_cellxgene.h5ad")
        age = a.obs.GestationalAge.astype(str).to_numpy()
        for label in a.obs.GestationalAge.astype(str).unique():
            for s in a.obs.SampleName.astype(str)[age == label].unique():
                i = np.flatnonzero((age == label) & (a.obs.SampleName.astype(str).to_numpy() == s))
                x = sp.csr_matrix(a.X[i]).expm1()  # counts recovered as in A242
                x = sp.diags(a.obs.nCount_RNA.to_numpy(float)[i] / np.asarray(x.sum(1)).ravel()) @ x
                yield "kildisiute", label, s, sp.csc_matrix(x.T).astype(np.float32), a.var_names.to_numpy(str), \
                    a.obs_names.to_numpy(str)[i]


def stage_markers():
    canon, rows = _canon(), []
    for cohort, fetus, sample, x, names, bc in samples():
        sym = np.array(canon(names))
        m = np.column_stack([np.asarray(x[np.flatnonzero(sym == g)].sum(0)).ravel() for g in ALL_MARKERS])
        missing = [g for g in ALL_MARKERS if not (sym == g).any()]
        t = pd.DataFrame(m.astype(np.float32), columns=ALL_MARKERS)
        t.insert(0, "total", np.asarray(x.sum(0)).ravel().astype(np.float32))
        t.insert(0, "barcode", bc)
        t.insert(0, "i", np.arange(len(bc)))
        t.insert(0, "sample", sample)
        t.insert(0, "fetus", f"{cohort}:{fetus}")
        t.insert(0, "cohort", cohort)
        rows.append(t)
        print(f"{cohort}:{fetus:12s} {sample:40s} cells {len(bc):6d} markers missing {missing}", flush=True)
    c = pd.concat(rows, ignore_index=True)
    c["published"] = ""
    k = ad.read_h5ad(SRC / "kildisiute2021" / "raw" / "adr_all_cellxgene.h5ad", backed="r").obs.Annotation.astype(str)
    m = c.cohort == "kildisiute"
    c.loc[m, "published"] = k.reindex(c.barcode[m]).to_numpy()
    d = pd.read_csv(SRC / "dong2020" / "raw" / "GSE137804_Adrenal_gland_annotation.csv.gz")
    d["key"] = d["sample"] + ":" + d.cell_id.str.split(".").str[0]
    key = c["sample"].str.split("_").str[1] + ":" + c.barcode.str.split(".").str[0]
    m = c.cohort == "dong"
    c.loc[m, "published"] = d.set_index("key").annotation.reindex(key[m]).to_numpy()
    c["published"] = c.published.fillna("")
    c.to_parquet(CELLS)
    print("cells", len(c), "published annotation:", c.groupby("cohort").published.apply(lambda s: (s != "").mean()).round(3).to_dict(),
          flush=True)


def fractions(c):
    """Marker fraction (class marker UMIs / total UMIs) and marker count of each class."""
    cnt = pd.DataFrame({k: c[v].sum(axis=1) for k, v in MARKERS.items()})
    return cnt.div(c.total, axis=0), cnt


def classify(c, rule=RULE):
    """Class of each cell: the one of cortex / medulla / stroma that it meets alone and not blood; '' otherwise
    (two classes, none, or blood). Reads only the marker columns and the total."""
    f, n = fractions(c)
    met = pd.DataFrame({k: (n[k] >= rule["min_count"]) & (f[k] >= rule[k]) for k in MARKERS})
    one = met[list(CLASSES)].sum(axis=1) == 1
    out = met[list(CLASSES)].idxmax(axis=1).where(one & ~met.blood, "")
    return out, met


def stage_hist():
    c = pd.read_parquet(CELLS)
    f, n = fractions(c)
    edges = np.concatenate([[0], np.logspace(-4, 0, 41)])
    rows = []
    for k in MARKERS:
        for g, idx in [("all", c.index), *c.groupby("cohort").groups.items()]:
            h = np.histogram(f.loc[idx, k], edges)[0]
            for e0, e1, v in zip(edges[:-1], edges[1:], h):
                rows.append(dict(marker=k, group=g, lo=e0, hi=e1, cells=int(v)))
    pd.DataFrame(rows).to_csv(OUT / "marker_hist.csv", index=False)
    print("cells", len(c))


def stage_classes():
    c = pd.read_parquet(CELLS)
    c["cls"], met = classify(c)
    c["blood"] = met.blood
    c["n_met"] = met[list(CLASSES)].sum(axis=1)
    sizes = c.groupby(["cohort", "fetus", "cls"]).size().unstack(fill_value=0)
    sizes = sizes.rename(columns={"": "left_out"})
    sizes["cells"] = c.groupby(["cohort", "fetus"]).size()
    sizes.to_csv(OUT / "class_sizes.csv")
    print(sizes.to_string())
    for coh, g in c[c.published != ""].groupby("cohort"):
        t = pd.crosstab(g.published, g.cls.replace("", "(left out)"))
        t.to_csv(OUT / f"agreement_{coh}.csv")
        print(coh, t.to_string(), sep="\n")


def stage_pseudobulk():
    c = pd.read_parquet(CELLS)
    c["cls"], _ = classify(c)
    canon, acc, whole = _canon(), {}, {}
    for cohort, fetus, sample, x, names, bc in samples():
        cc = c[(c["sample"] == sample) & (c.fetus == f"{cohort}:{fetus}")].sort_values("i")
        assert len(cc) == x.shape[1] and (cc.barcode.to_numpy() == bc).all(), sample
        for k, mask in [("all", np.ones(len(cc), bool)), *((k, (cc.cls == k).to_numpy()) for k in CLASSES)]:
            if not mask.any():
                continue
            s = a242.by_symbol(np.asarray(x[:, mask].sum(1)).ravel(), names)
            key = (f"{cohort}:{fetus}", k)
            acc[key] = (acc[key][0].add(s, fill_value=0), acc[key][1] + int(mask.sum())) if key in acc else (s, int(mask.sum()))
        print(f"{cohort}:{fetus} {sample}", flush=True)
    ref = pd.read_parquet(DEST / "pseudobulk.parquet")
    for fetus in {k[0] for k in acc}:  # the classes and the cells left out add up to the A242 whole-adrenal pseudobulk
        s = acc[(fetus, "all")][0]
        r = ref.loc[fetus].dropna()
        assert np.allclose(s.reindex(r.index).to_numpy(), r.to_numpy(), rtol=1e-3, atol=0.5), fetus
    keep = {k: v for k, v in acc.items() if k[1] != "all" and v[1] >= MIN_CELLS}
    sheet = pd.DataFrame([dict(id=f"{f}|{k}", donor=f, compartment=k, cohort=f.split(":")[0], cells=n,
                               depth=float(s.sum())) for (f, k), (s, n) in keep.items()])
    allcells = pd.DataFrame([dict(fetus=f, cls=k, cells=n, depth=float(s.sum())) for (f, k), (s, n) in acc.items()])
    allcells.to_csv(OUT / "class_pseudobulk_all_sizes.csv", index=False)
    pb = pd.DataFrame({f"{f}|{k}": s for (f, k), (s, _) in acc.items() if k != "all"}).T.sort_index()
    pb.columns = pb.columns.astype(str)
    pb.to_parquet(PB)
    s1 = pd.read_csv(DEST / "samples.csv").set_index("id")
    sheet["age"] = sheet.donor.map(s1.age)
    sheet.to_csv(SHEET, index=False)
    print(sheet.groupby(["cohort", "compartment"]).agg(fetuses=("donor", "size"), cells=("cells", "median")).to_string())


OTHER = {"cortex": "medulla", "medulla": "cortex"}


def foreign(pb, sheet, genes):
    """RNA of the other class in the cortex and medulla pseudobulks of a fetus (the A190 limit; as
    `score_a241.contamination`). a = (other-class marker counts / depth) of the pseudobulk, relative to the same
    fraction in the other class's pseudobulk of the fetus (any size; its cells are reported). The expected foreign
    clock-gene counts are a x depth x (clock fraction of the other class); `corrected` has them subtracted (floored at
    0, rounded). Returns (table per fetus and class, corrected counts of the clock genes)."""
    rows, corr = [], pb.loc[sheet.id, genes].copy()
    for _, r in sheet[sheet.compartment.isin(OTHER)].iterrows():
        ref = pb.loc[f"{r.donor}|{OTHER[r.compartment]}"] if f"{r.donor}|{OTHER[r.compartment]}" in pb.index else None
        if ref is None:
            continue
        x, mk = pb.loc[r.id], MARKERS[OTHER[r.compartment]]
        a = (x[mk].sum() / x.sum()) / (ref[mk].sum() / ref.sum())
        exp = a * x.sum() * ref[genes] / ref.sum()
        corr.loc[r.id] = np.round(np.clip(x[genes].to_numpy(float) - exp.to_numpy(), 0, None))
        rows.append(dict(id=r.id, cohort=r.cohort, donor=r.donor, compartment=r.compartment, a=a,
                         ref_cells=int(pb.attrs["cells"][(r.donor, OTHER[r.compartment])]),
                         ref_ge_min=bool(f"{r.donor}|{OTHER[r.compartment]}" in set(sheet.id)),
                         share_median=float(np.median(exp / x[genes].astype(float))),
                         foreign_counts=float(exp.sum()), clock_counts=float(x[genes].sum()),
                         **{f"exp_{g}": exp[g] for g in genes}))
    return pd.DataFrame(rows), corr


def stage_design():
    import clocksync as cs
    genes = list(cs.gtex().genes)
    pb, sheet = pd.read_parquet(PB), pd.read_csv(SHEET)
    sizes = pd.read_csv(OUT / "class_pseudobulk_all_sizes.csv")
    pb.attrs["cells"] = {(f, k): n for f, k, n in zip(sizes.fetus, sizes.cls, sizes.cells)}
    f, _ = foreign(pb, sheet, genes)
    f.drop(columns=[c for c in f if c.startswith("exp_")]).to_csv(OUT / "foreign_fetus.csv", index=False)
    rows = []
    for (coh, k), g in f.groupby(["cohort", "compartment"]):
        e, x = g[[c for c in g if c.startswith("exp_")]].sum(), pb.loc[g.id, genes].sum()
        share = pd.Series(e.to_numpy() / x.to_numpy(), index=genes)
        rows.append(dict(cohort=coh, compartment=k, n_fetuses=len(g), n_ref_ge_min=int(g.ref_ge_min.sum()),
                         a_median=g.a.median(), a_min=g.a.min(), a_max=g.a.max(), share_median_over_genes=share.median(),
                         share_median_over_fetuses=g.share_median.median(), share_max_gene=share.max(),
                         share_max_gene_name=share.idxmax()))
    for k, g in f.groupby("compartment"):  # the whole class (the confounding rule: median foreign share > 0.5)
        e, x = g[[c for c in g if c.startswith("exp_")]].sum(), pb.loc[g.id, genes].sum()
        share = pd.Series(e.to_numpy() / x.to_numpy(), index=genes)
        rows.append(dict(cohort="all", compartment=k, n_fetuses=len(g), n_ref_ge_min=int(g.ref_ge_min.sum()),
                         a_median=g.a.median(), a_min=g.a.min(), a_max=g.a.max(), share_median_over_genes=share.median(),
                         share_median_over_fetuses=g.share_median.median(), share_max_gene=share.max(),
                         share_max_gene_name=share.idxmax()))
    pd.DataFrame(rows).to_csv(OUT / "foreign_cohort.csv", index=False)
    x = pb.loc[sheet.id, genes]
    det = x.groupby([sheet.cohort.to_numpy(), sheet.compartment.to_numpy()]).median().round(1)
    det.index.names = ["cohort", "compartment"]
    det.to_csv(OUT / "clock_detection_classes.csv")
    d = sheet.assign(clock_total=x.sum(axis=1).to_numpy(), clock_median=x.median(axis=1).to_numpy()).groupby(
        ["compartment", "cohort"], sort=False).agg(
        n_fetuses=("id", "size"), ages=("age", lambda v: " ".join(f"{a:.1f}" for a in sorted(v))), cells_median=("cells", "median"),
        cells_min=("cells", "min"), depth_median=("depth", "median"), clock_gene_median_count=("clock_median", "median"),
        clock_total_median=("clock_total", "median")).reset_index()
    a = f.groupby(["compartment", "cohort"]).a.median().rename("a_median").reset_index()
    d = d.merge(a, how="left")
    d.to_csv(OUT / "design_classes.csv", index=False)
    pd.set_option("display.width", 250)
    print(d.round(3).to_string(index=False), pd.read_csv(OUT / "foreign_cohort.csv").round(4).to_string(index=False), sep="\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True, choices=["markers", "hist", "classes", "pseudobulk", "design"])
    globals()[f"stage_{ap.parse_args().stage}"]()
