"""Whole-organ pseudobulks of open human fetal heart, kidney, lung and liver (A243): one row per fetus, organ, cohort.

Sources (`$SCCC_DATA/natexp/<folder>/raw/`, each with a PROVENANCE.md; processed matrices only; downloads by
`download_a243.sh`):
  cao2020              GSE156793 (on disk): sci-RNA-seq3 organ pseudobulks (genome-wide), organs with >= 1,000 cells
  cardoso_moreira2019  E-MTAB-6814 (on disk): bulk RNA-seq, prenatal libraries with a resolved individual (as A237);
                       counts recovered from the released CPM (CPM / smallest non-zero CPM of the library)
  roost2015            GSE66302 (on disk): DeepSAGE tag counts; atrium and ventricle of one fetus summed
  han2020_organs       GSE134355: Microwell-seq dge tables (UMI), heart 2, kidney 4, lung 2 (liver not used: E5/E3)
  encode_fetal_organs  ENCODE UW polyA+ (heart 2, kidney 4) and CSHL total RNA (heart 2, lung 2, liver 2): RSEM
                       expected counts
  cui2019              GSE106118: plate scRNA-seq (UMI) of 20 fetal hearts (all cells of a heart summed)
  knight_schrijver2023 GSE216019: 10x, apex + base libraries of 6 fetal hearts summed
  hochane2019          GSE114530: 10x, 5 whole fetal kidneys (raw matrices: barcodes with >= 500 UMIs)
  lindstrom_gse*       GSE124472, GSE112570, GSE139280, GSE102596: 10x, dissected kidney cortex of 4 fetuses
  sountoulidis2023     GSE215895: 10x, whole fetal lungs (libraries of one donor summed)
  renesme2025          GSE280880: snRNA-seq UMI matrix of 9 fetal lungs

Not used: Renesme 2025 (GSE280880): the UMI matrix has 24,441 rows without gene names or identifiers (E2: genes
not identifiable); Quach (GSE264407): Seurat .rds only (not parsed); He 2022 (E-MTAB-11278): CD326-enriched samples.

Unit: a fetus (one whole-organ pseudobulk per fetus and cohort). Samples of one source with the same stated age are
one fetus and are summed unless the source names the individuals (Cui: heart names; Han: sample names; Roost,
Knight-Schrijver, Sountoulidis: donor codes). Age in weeks as stated (days / 7 where given in days); the convention is
recorded per source. Symbols through `clocksync.genes.canonical`; several rows of one symbol are summed; NaN = the
symbol is absent from the source.

Computed: sums of counts, cell numbers, depths, the median count of each clock gene per organ and cohort (E4
detection). NOT computed: any clock statistic.

    python -I prep_a243.py        (downloads are parsed, never executed)

Outputs: `$SCCC_DATA/natexp/fetal_organs_pooled/{pseudobulk.parquet, samples.csv}`;
`../output/{design_cohorts.csv, clock_detection.csv}`.
"""
import gzip
import io
import re
import sys
import tarfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.io as sio  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "natexp"
DEST = SRC / "fetal_organs_pooled"
OUT = HERE.parents[1] / "output"
GENES = list(cs.gtex().genes)
MIN_UMI = 500  # raw 10x matrices: barcodes kept
MIN_CELLS_CAO = 1000
E4_MIN = 20
ORGANS = ["heart", "kidney", "lung", "liver"]
ENS = {}  # Ensembl gene id (no version) -> symbol


def by_symbol(values, symbols) -> pd.Series:
    """Counts summed per canonical symbol."""
    s = pd.Series(np.asarray(values, float), index=[canonical(g) for g in symbols])
    return s.groupby(level=0).sum()


def from_ensembl(values, ids) -> pd.Series:
    sym = pd.Series([ENS.get(str(i).split(".")[0]) for i in ids])
    ok = sym.notna().to_numpy()
    return by_symbol(np.asarray(values, float)[ok], sym[ok])


def mtx(matrix_fh, genes_fh, barcodes_fh=None, min_umi: int = 0):
    """(counts per symbol, cells) of a 10x matrix.mtx(.gz) with features/genes tsv (genes x barcodes)."""
    x = sp.csc_matrix(sio.mmread(matrix_fh))
    g = pd.read_csv(genes_fh, sep="\t", header=None)
    ids, names = g.iloc[:, 0].astype(str), g.iloc[:, 1].astype(str)
    ENS.update(zip(ids.str.split(".").str[0], names))
    if x.shape[1] > 30000:  # raw matrix (empty droplets included)
        min_umi = max(min_umi, MIN_UMI)
    keep = np.asarray(x.sum(0)).ravel() >= max(min_umi, 1)
    return by_symbol(np.asarray(x[:, keep].sum(1)).ravel(), names), int(keep.sum())


def tar_member(tar, name):
    return gzip.open(tar.extractfile(name), "rb")


def wide(path, sep="\t", min_umi: int = 0):
    """(counts per symbol, cells) of a gene x cell text table (gz)."""
    t = pd.read_csv(path, sep=sep, index_col=0)
    keep = t.sum(axis=0).to_numpy() >= max(min_umi, 1)
    return by_symbol(t.loc[:, keep].sum(axis=1).to_numpy(), t.index), int(keep.sum())


class Collector:
    def __init__(self):
        self.rows, self.pb = [], {}

    def add(self, cohort, organ, fetus, age, stated, convention, assay, accession, counts, cells=np.nan, sample="",
            note=""):
        key = f"{cohort}:{organ}:{fetus}"
        self.pb[key] = counts if key not in self.pb else self.pb[key].add(counts, fill_value=0)
        self.rows.append(dict(id=key, cohort=cohort, organ=organ, fetus=str(fetus), age_weeks=float(age),
                              age_stated=stated, convention=convention, assay=assay, accession=accession,
                              sample=sample, cells=cells, depth=float(counts.sum()), note=note))
        print(f"{key:34s} {sample[:40]:40s} cells {cells} depth {counts.sum():.3g}", flush=True)


def cao(c):
    g = ad.read_h5ad(SRC / "cao2020" / "genome" / "pseudobulk_organ.h5ad")
    meta = pd.read_csv(SRC / "cao2020" / "pseudobulk_organ" / "samples.csv").set_index(["donor_id", "organ"])
    x = g.X.toarray() if sp.issparse(g.X) else np.asarray(g.X)
    for i, (donor, comp, n) in enumerate(zip(g.obs.donor, g.obs.compartment, g.obs.n_cells)):
        if comp.lower() not in ORGANS or n < MIN_CELLS_CAO:
            continue
        m = meta.loc[(donor, comp)]
        # the encode_uw / Cao pairs of E6: heart H27098 (120 d M), kidney H27798 (85 d F)
        note = "possibly the fetus of ENCODE UW " + {"H27098": "ENCSR047LLJ (heart 120 d M)",
                                                    "H27798": "ENCSR896QPD (kidney 85 d F)"}.get(donor, "") \
            if donor in ("H27098", "H27798") else ""
        c.add("cao2020", comp.lower(), donor, m.ga_days / 7, f"{m.ga_days:.0f} days", "post conception (ga_pcw)",
              "sci-RNA-seq3", "GSE156793", by_symbol(x[i], g.var_names), int(n), f"{comp}_{donor}", note)


def cardoso(c):
    d = SRC / "cardoso_moreira2019"
    s = pd.read_csv(d / "individuals.csv").set_index("sample")
    s = s[s.prenatal & (s.link != "open") & s.organ.isin(ORGANS)]
    major = s.groupby("individual").sex.agg(lambda v: v.mode().iloc[0])
    s = s[s.sex == s.individual.map(major)]  # as A237: a library of the other sex than its individual is dropped
    cpm = pd.read_parquet(d / "expr_cpm.parquet")[s.index]
    sym = pd.read_csv(d / "genes.csv").set_index("gene_id").symbol
    ENS.update({k: v for k, v in sym.items() if isinstance(v, str)})
    x = cpm.to_numpy(float)
    x = x / np.where(x > 0, x, np.inf).min(0)  # counts: CPM / smallest non-zero CPM of the library
    for j, name in enumerate(s.index):
        r = s.loc[name]
        c.add("cardoso_moreira", r.organ, f"i{r.individual:.0f}", r.age_wpc, f"{r.age_wpc:g} wpc", "post conception",
              "bulk RNA-seq", "E-MTAB-6814", from_ensembl(np.rint(x[:, j]), cpm.index), sample=name,
              note=f"sex {r.sex}; counts recovered from CPM")


def roost(c):
    t = pd.read_csv(SRC / "roost2015" / "raw" / "GSE66302_Roost_et_al_Raw_Reads.txt.gz", sep="\t", index_col=0)
    gm = SRC / "roost2015_organs" / "raw" / "GSE66302_series_matrix.txt.gz"
    meta = {}  # sample title -> gestational age
    titles = ages = None
    with gzip.open(gm, "rt") as h:
        for line in h:
            if line.startswith("!Sample_title"):
                titles = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_characteristics_ch1") and "gestational age: " in line:
                ages = [float(re.sub(r"[^\d.]", "", v.split(": ")[1]) or "nan") for v in line.rstrip("\n").split("\t")[1:]]
    meta = dict(zip(titles, ages))
    for col in t.columns:
        m = re.match(r"(\w+) (heart [AV]|heart|kidney|lung|liver)_(\d+)$", col)
        if not m:
            continue
        fetus, part, _ = m.groups()
        organ = part.split()[0]
        c.add("roost", organ, fetus, meta[col], f"{meta[col]:g} weeks (gestational age)", "gestational",
              "DeepSAGE tags", "GSE66302", from_ensembl(t[col], t.index), sample=col,
              note="chamber summed" if organ == "heart" else "")


HAN = {"GSM4008686_Fetal-Heart1": ("heart", "Heart1", 12, "F"), "GSM4008687_Fetal-Heart2": ("heart", "Heart2", 11, "F"),
       "GSM4008693_Fetal-Kidney3": ("kidney", "Kidney3", 13, "M"), "GSM4008694_Fetal-Kidney4": ("kidney", "Kidney4", 11, "M"),
       "GSM4008695_Fetal-Kidney5": ("kidney", "Kidney5", 12, "M"), "GSM4008696_Fetal-Kidney6": ("kidney", "Kidney6", 11, "F"),
       "GSM4008699_Fetal-Lung1": ("lung", "Lung1", 12, "F"), "GSM4008700_Fetal-Lung2": ("lung", "Lung2", 11, "F")}


def han(c):
    for name, (organ, fetus, age, sex) in HAN.items():
        cnt, n = wide(SRC / "han2020_organs" / "raw" / f"{name}_dge.txt.gz")
        c.add("han_hcl", organ, fetus, age, f"{age}-week-old", "not stated", "Microwell-seq", "GSE134355", cnt, n,
              name, f"sex {sex}")


# file, cohort, organ, donor, age stated (days or weeks), sex, assay, experiment
ENCODE = [("ENCFF646WRS", "encode_uw", "heart", "ENCDO605ZCF", 91, "day", "F", "ENCSR863BUL"),
          ("ENCFF292LKT", "encode_uw", "heart", "ENCDO249ZDZ", 120, "day", "M", "ENCSR047LLJ"),
          ("ENCFF564JYO", "encode_uw", "kidney", "ENCDO165PKL", 85, "day", "F", "ENCSR896QPD"),
          ("ENCFF193UXG", "encode_uw", "kidney", "ENCDO443KWV", 87, "day", "M", "ENCSR907KDH"),
          ("ENCFF012IBV", "encode_uw", "kidney", "ENCDO718LHH", 105, "day", "F", "ENCSR212AMA"),
          ("ENCFF598STD", "encode_uw", "kidney", "ENCDO012PVK", 108, "day", "F", "ENCSR495UXA"),
          ("ENCFF589BBG", "encode_cshl", "heart", "ENCDO069AAA", 28, "week", "F", "ENCSR000AEZ"),
          ("ENCFF181TFL", "encode_cshl", "heart", "ENCDO070AAA", 19, "week", "F", "ENCSR000AEZ"),
          ("ENCFF675DVW", "encode_cshl", "lung", "ENCDO064AAA", 20, "week", "F", "ENCSR000AFC"),
          ("ENCFF678NWU", "encode_cshl", "lung", "ENCDO067AAA", 24, "week", "F", "ENCSR000AFC"),
          ("ENCFF780YUS", "encode_cshl", "liver", "ENCDO064AAA", 20, "week", "F", "ENCSR000AFB"),
          ("ENCFF846SOK", "encode_cshl", "liver", "ENCDO065AAA", 22, "week", "M", "ENCSR000AFB")]
E6_NOTE = {"ENCSR047LLJ": "same age and sex as Cao heart H27098: possibly one fetus, kept distinct (as A242)",
           "ENCSR896QPD": "same age and sex as Cao kidney H27798: possibly one fetus, kept distinct (as A242)"}


def encode(c):
    for f, cohort, organ, donor, age, unit, sex, exp in ENCODE:
        t = pd.read_csv(SRC / "encode_fetal_organs" / "raw" / f"{f}.tsv", sep="\t")
        c.add(cohort, organ, donor, age / 7 if unit == "day" else age, f"{age} {unit}", "not stated (ENCODE)",
              "polyA+ RNA-seq" if cohort == "encode_uw" else "total RNA-seq", exp,
              from_ensembl(t.expected_count, t.gene_id), sample=f, note="; ".join(
                  x for x in (f"sex {sex}; RSEM expected counts", E6_NOTE.get(exp, "")) if x))


def cui(c):
    d = SRC / "cui2019" / "raw"
    t = pd.read_csv(d / "GSE106118_UMI_count_merge.txt.gz", sep="\t", index_col=0)
    heart = t.columns.str.replace(r"_[A-Z]+\.\d+$", "", regex=True)  # HE17W_2_LA.5 -> HE17W_2
    for h in pd.unique(heart):
        age = float(re.match(r"HE(\d+)W", h).group(1))
        cols = heart == h
        c.add("cui", "heart", h, age, f"{age:g} W", "not stated", "plate scRNA-seq (UMI)", "GSE106118",
              by_symbol(t.loc[:, cols].sum(axis=1).to_numpy(), t.index), int(cols.sum()), h,
              "cells picked from chambers/valves, all summed")


KS_AGE = {"F1": 58, "F2": 61, "F3": 64, "F4": 66, "F6": 69, "F7": 80}


def knight_schrijver(c):
    with tarfile.open(SRC / "knight_schrijver2023" / "raw" / "GSE216019_RAW.tar") as tar:
        names = [m.name for m in tar.getmembers() if m.name.endswith("_matrix.mtx.gz")]
        for n in sorted(names):
            fetus = re.search(r"_(F\d)_", n).group(1)
            stem = n[:-len("matrix.mtx.gz")]
            cnt, cells = mtx(tar_member(tar, n), tar_member(tar, stem + "features.tsv.gz"), min_umi=0)
            c.add("knight_schrijver", "heart", fetus, KS_AGE[fetus] / 7, f"{KS_AGE[fetus]} days", "not stated",
                  "10x 3' v3", "GSE216019", cnt, cells, n[:-len("_matrix.mtx.gz")], "apex + base summed")


def hochane(c):
    ages = {"w16": 16, "wk09": 9, "wk11": 11, "wk13": 13, "wk18": 18}
    with tarfile.open(SRC / "hochane2019" / "raw" / "GSE114530_RAW.tar") as tar:
        for n in sorted(m.name for m in tar.getmembers() if m.name.endswith("_matrix.mtx.gz")):
            tag = re.search(r"_(w16|wk\d+)_matrix", n).group(1)
            stem = n[:-len("matrix.mtx.gz")]
            cnt, cells = mtx(tar_member(tar, n), tar_member(tar, stem + "genes.tsv.gz"), min_umi=0)
            c.add("hochane", "kidney", tag, ages[tag], f"Week {ages[tag]}", "not stated", "10x 3' v2", "GSE114530",
                  cnt, cells, n[:-len("_matrix.mtx.gz")], f"barcodes >= {MIN_UMI} UMIs")


def lindstrom(c):
    """Four fetuses (14, 15, 16, 17 weeks): dissected kidney cortex, one cohort. The 17-week fetus of GSE124472 and
    GSE112570 is one fetus by E6 (same age, same laboratory): `finish` keeps the deeper source."""
    note = "dissected kidney cortex"
    d = SRC / "lindstrom_gse124472" / "raw"
    with tarfile.open(d / "GSE124472_RAW.tar") as tar:  # nested archives, read in memory (matrix, genes, barcodes)
        for gsm, tag, age in (("GSM3534656_H17w_Z1", "Z1", 17), ("GSM3534657_H17w_Z2", "Z2", 17),
                              ("GSM3534658_H15w_Z1", "Z1", 15), ("GSM3534659_H15w_Z2", "Z2", 15)):
            with tarfile.open(fileobj=tar.extractfile(f"{gsm}_raw_counts.tar.gz"), mode="r:gz") as inner:
                cnt, cells = mtx(inner.extractfile("matrix.mtx"), inner.extractfile("genes.tsv"))
            c.add("lindstrom_combes", "kidney", f"{age}w:GSE124472", age, f"Week {age}", "not stated", "10x 3'",
                  "GSE124472", cnt, cells, gsm, f"{note} (outer + inner zone summed)")
    with tarfile.open(SRC / "lindstrom_gse112570" / "raw" / "GSE112570_RAW.tar") as tar:
        for gsm in ("GSM3073088_H17w_K1", "GSM3073089_H17w_K2"):
            cnt, cells = mtx(tar_member(tar, f"{gsm}_matrix.mtx.gz"), tar_member(tar, f"{gsm}_genes.tsv.gz"))
            c.add("lindstrom_combes", "kidney", "17w:GSE112570", 17, "Week 17", "not stated", "10x 3'", "GSE112570",
                  cnt, cells, gsm, f"{note} (two kidneys of one stated age summed)")
    with tarfile.open(SRC / "lindstrom_gse139280" / "raw" / "GSE139280_RAW.tar") as tar:
        # kidney 1 (fresh, 0 h) only: kidney 2 of the same fetus was stored 16 h at 4 C (a storage test)
        z = zipfile.ZipFile(io.BytesIO(gzip.open(tar.extractfile("GSM4135994_count_table_H14w_0h.zip.gz")).read()))
        cnt, cells = mtx(z.open("matrix.mtx"), z.open("genes.tsv"))
        c.add("lindstrom_combes", "kidney", "14w:GSE139280", 14, "Week 14", "not stated", "10x 3'", "GSE139280", cnt,
              cells, "GSM4135994_H14w_0h", f"{note}; the 16 h stored kidney 2 not used")
    with tarfile.open(SRC / "lindstrom_gse102596" / "raw" / "GSE102596_RAW.tar") as tar:
        t = pd.read_csv(tar.extractfile("GSM2741551_count-table-human16w.tsv.gz"), sep="\t", index_col=0,
                        compression="gzip")
        c.add("lindstrom_combes", "kidney", "16w:GSE102596", 16, "week 16", "not stated", "10x 3'", "GSE102596",
              by_symbol(t.sum(axis=1).to_numpy(), t.index), t.shape[1], "GSM2741551_human16w", note)


def tenx_h5(fh):
    """(counts per symbol, cells) of a filtered 10x h5 (file object)."""
    with h5py.File(fh, "r") as f:
        m = f["matrix"]
        x = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
        names, ids = m["features/name"][:].astype(str), m["features/id"][:].astype(str)
    ENS.update(zip(pd.Series(ids).str.split(".").str[0], names))
    return by_symbol(np.asarray(x.sum(1)).ravel(), names), x.shape[1]


def sountoulidis(c):
    """Libraries `10X<donor>w<k>` of one donor are one fetus (donor numbers name the individuals); age from the GEO
    sample titles ("lung,PCW8.5, rep1 [10X140w1]")."""
    age = {}
    with gzip.open(SRC / "sountoulidis2023" / "raw" / "GSE215895_series_matrix.txt.gz", "rt") as h:
        for line in h:
            if line.startswith("!Sample_title"):
                for v in line.rstrip("\n").split("\t")[1:]:
                    m = re.search(r"PCW([\d.]+).*\[10[Xx](\d+)w(\d+)\]", v)
                    age[(m.group(2), m.group(3))] = float(m.group(1))
    with tarfile.open(SRC / "sountoulidis2023" / "raw" / "GSE215895_RAW.tar") as tar:
        for n in sorted(m.name for m in tar.getmembers() if m.name.endswith("_counts.h5")):
            donor, k = re.search(r"10[Xx](\d+)w(\d+)_counts", n).groups()
            cnt, cells = tenx_h5(tar.extractfile(n))
            c.add("sountoulidis", "lung", donor, age[(donor, k)], f"PCW{age[(donor, k)]:g}", "post conception",
                  "10x 3'", "GSE215895", cnt, cells, n.split("_")[0] + f"_{donor}w{k}", "libraries of one donor summed")


def finish(c, out: bool = True):
    s = pd.DataFrame(c.rows)
    # E6, Lindstrom/Combes 17-week fetus (GSE124472 and GSE112570): the deeper source is kept
    pair = s[s.id.isin(["lindstrom_combes:kidney:17w:GSE124472", "lindstrom_combes:kidney:17w:GSE112570"])]
    if pair.id.nunique() == 2:
        depth = pair.groupby("id").depth.sum()
        drop = depth.idxmin()
        print(f"E6 17-week Lindstrom/Combes fetus: kept {depth.idxmax()}, dropped {drop}", flush=True)
        s = s[s.id != drop]
        c.pb.pop(drop)
        s.loc[s.id == depth.idxmax(), "note"] += f"; E6: the same 17-week fetus as {drop.split(':')[-1]} (shallower), counted once"
    assert (s.groupby("id").age_weeks.nunique() == 1).all(), "samples of one fetus differ in age"
    sheet = s.groupby("id", sort=False).agg(
        cohort=("cohort", "first"), organ=("organ", "first"), fetus=("fetus", "first"), age_weeks=("age_weeks", "first"),
        age_stated=("age_stated", "first"), convention=("convention", "first"), assay=("assay", "first"),
        accession=("accession", "first"), n_samples=("sample", "size"), samples=("sample", "+".join),
        cells=("cells", "sum"), depth=("depth", "sum"), note=("note", lambda v: "; ".join(sorted(set(x for x in v if x)))))
    sheet.loc[sheet.cells == 0, "cells"] = np.nan  # bulk
    sheet = sheet.reset_index()
    counts = pd.DataFrame(c.pb).T.loc[sheet.id]
    counts.columns = counts.columns.astype(str)
    det = counts[GENES].groupby([sheet.organ.to_numpy(), sheet.cohort.to_numpy()], sort=False).median().round(1)
    det.index.names = ["organ", "cohort"]
    det["E4_median_of_medians"] = det[GENES].median(axis=1)
    det["E4_pass"] = det.E4_median_of_medians >= E4_MIN
    nfet = sheet.groupby(["organ", "cohort"], sort=False).id.size()
    det["E5_pass"] = (nfet.reindex(det.index) >= 2).to_numpy()
    det["n_fetuses"] = nfet.reindex(det.index).to_numpy()
    sheet = sheet.merge(det[["E4_pass", "E5_pass"]].reset_index(), on=["organ", "cohort"])
    sheet["eligible"] = sheet.E4_pass & sheet.E5_pass
    sheet["window"] = np.where(sheet.age_weeks >= 12, "from 12 weeks", "before 12 weeks")
    design = sheet.groupby(["organ", "cohort"], sort=False).agg(
        accession=("accession", lambda v: "+".join(sorted(set(v)))), assay=("assay", "first"),
        age_convention=("convention", lambda v: "; ".join(sorted(set(v)))), n_fetuses=("id", "size"),
        n_before_12w=("window", lambda v: int((v == "before 12 weeks").sum())),
        n_from_12w=("window", lambda v: int((v == "from 12 weeks").sum())),
        age_min=("age_weeks", "min"), age_max=("age_weeks", "max"), median_cells=("cells", "median"),
        median_depth=("depth", "median")).reset_index()
    design = design.merge(det[["E4_median_of_medians", "E4_pass", "E5_pass"]].reset_index(), on=["organ", "cohort"])
    design["organ"] = pd.Categorical(design.organ, ORGANS)
    design = design.sort_values(["organ", "cohort"]).round({"age_min": 2, "age_max": 2, "median_depth": 0})
    if out:
        DEST.mkdir(exist_ok=True)
        OUT.mkdir(exist_ok=True)
        counts.to_parquet(DEST / "pseudobulk.parquet")
        sheet.to_csv(DEST / "samples.csv", index=False)
        design.to_csv(OUT / "design_cohorts.csv", index=False)
        det.reset_index().to_csv(OUT / "clock_detection.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(design.to_string(index=False), det.drop(columns=["n_fetuses"]).to_string(), f"symbols {counts.shape[1]}",
          sep="\n", flush=True)
    return sheet, counts, design, det


STEPS = dict(cardoso=cardoso, cao=cao, han=han, encode=encode, knight_schrijver=knight_schrijver, hochane=hochane,
             cui=cui, roost=roost, lindstrom=lindstrom, sountoulidis=sountoulidis)


def main():
    steps = sys.argv[1:] or list(STEPS)  # a subset prints the tables without writing the outputs
    c = Collector()
    for k in steps:
        STEPS[k](c)
    finish(c, out=not sys.argv[1:])


if __name__ == "__main__":
    main()
