"""Whole-adrenal pseudobulks of the open human fetal adrenal transcriptomes (A242): one row per fetus and cohort.

Sources (`$SCCC_DATA/natexp/<folder>/raw/`, each with a PROVENANCE.md; processed matrices only):
  encode_fetal_adrenal  ENCODE (UW): long-read RNA-seq of 7 fetuses (transcript quantifications, summed per gene),
                        polyA plus RNA-seq of 4 of them and total RNA-seq of an eighth (RSEM expected counts)
  kameneva2021          GSE147821, 10x raw matrices of 11 adrenal samples of 11 fetuses (barcodes with >= 500 UMIs)
  dong2020              GSE137804, 10x gene x cell tables of 4 fetal adrenals
  wang2022              GSE167860, 10x filtered matrices of 10 fetal adrenals
  han2020               GSE134355, Microwell-seq dge tables of 3 fetal adrenal samples
  cheng2022, gse291218  GSE193672 and GSE291218 (one laboratory), 10x filtered matrices; three GSE291218 files hold the
                        barcodes of a Cheng 2022 sample and are dropped (`DUPLICATE`)
  roost2015             GSE66302, DeepSAGE tag counts of 3 fetal adrenals
  kildisiute2021        cellxgene object of 7 adrenal glands (log-normalised; counts recovered as expm1(x) rescaled
                        to the cell's nCount_RNA)
  schroder2022          GSE196700, targeted (amplicon) transcriptome counts of 13 fetal adrenals of unstated age

Unit: a fetus. Samples of one source with the same stated age are one fetus and are summed, unless the source names
them as different individuals (Kameneva: 11 fetuses; Wang: numbered fetuses; Cheng: embryo numbers). Age: weeks as
stated by the source (ENCODE: days / 7). Every cell of a sample is used (whole tissue, as the organ-level pseudobulk
of A189). Symbols through `clocksync.genes.canonical`; several rows of one symbol are summed.

Computed: sums of counts, cell numbers, depths, the median count of the clock genes per cohort (detection).
NOT computed: any clock statistic.

    python -I prep_adrenal.py        (downloads are parsed, never executed)

Outputs: `$SCCC_DATA/natexp/fetal_adrenal_open/{pseudobulk.parquet, samples.csv}` (fetus x symbol; NaN = symbol absent
from the source); `../output/{design_cohorts.csv, clock_detection.csv}`.
"""
import gzip
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "natexp"
DEST = SRC / "fetal_adrenal_open"
OUT = HERE.parents[1] / "output"
GENES = list(cs.gtex().genes)
MIN_UMI = 500  # Kameneva raw matrices: barcodes kept

# cohort, fetus, age (weeks as stated), convention, file
ENCODE_LR = [("H23365", 108, "ENCFF799IIR"), ("H24409", 113, "ENCFF066HKE"), ("H23887", 108, "ENCFF039THD"),
             ("H24272", 85, "ENCFF867GXK"), ("H23327", 85, "ENCFF743CJH"), ("H24279", 101, "ENCFF773OPR"),
             ("H23769", 108, "ENCFF059IFS")]
ENCODE_SR = [("encode_pa", "H23887", 108, "ENCFF295BLS"), ("encode_pa", "H24272", 85, "ENCFF773UYG"),
             ("encode_pa", "H24279", 101, "ENCFF151YRN"), ("encode_pa", "H23769", 108, "ENCFF358WTK"),
             ("encode_total", "H22662", 96, "ENCFF517LIG")]
KAMENEVA = {"GSM4446535_10X_19_001": 8, "GSM4446536_10X_19_063": 9, "GSM4446537_10X_19_088": 6,
            "GSM4446538_10X_19_123": 14, "GSM4446539_10X_19_124": 12, "GSM4446540_10X_19_125": 8,
            "GSM4446541_10X_20_005": 9, "GSM4446542_10X_20_006": 11, "GSM4446543_10X_20_007": 9,
            "GSM4734601_10X_20_016": 8, "GSM4734603_10X_20_035": 12}
DONG = {"GSM4088785_F2": 8 + 6 / 7, "GSM4088786_F7": 8 + 3 / 7, "GSM4088787_F106": 14 + 4 / 7,
        "GSM4088788_F107": 13 + 2 / 7}
WANG = {"GSM5114030_Sample1": ("7GW_male1", 7), "GSM5114031_Sample2": ("7GW_male2", 7),
        "GSM5114032_Sample3": ("7GW_female1", 7), "GSM5114035_Sample6": ("9GW_male1", 9),
        "GSM5114036_Sample7": ("9GW_female1", 9), "GSM5114039_Sample10": ("10GW_male1", 10),
        "GSM5114040_Sample11": ("10GW_male2", 10), "GSM5114042_Sample13": ("13GW_male1", 13),
        "GSM5114043_Sample14": ("13GW_female1", 13), "GSM5114046_Sample17": ("14GW_female1", 14)}
HAN = {"GSM4008675_Fetal-Adrenal-Gland2": 12, "GSM4008676_Fetal-Adrenal-Gland3": 14,
       "GSM4008677_Fetal-Adrenal-Gland4": 12}
SASAKI = {"cheng2022/raw/GSM5815985_HS30H_adrenal": 5, "cheng2022/raw/GSM5815987_HS26H_adrenal": 6,
          "cheng2022/raw/GSM5815989_HS12H_adrenal": 8, "cheng2022/raw/GSM5815990_Hs46_Adrenal": 9,
          "cheng2022/raw/GSM5815992_Hs43_Fresh": 12, "gse291218/raw/GSM8831091_fetal_adrenal_17wpf": 17,
          "gse291218/raw/GSM8831095_fetal_adrenal_5wpf": 5, "gse291218/raw/GSM8831096_fetal_adrenal_19wpf": 19}
DUPLICATE = {"GSM8831100_fetal_adrenal_6wpf": "GSM5815987_HS26H_adrenal",  # same barcodes (checked in main)
             "GSM8831101_fetal_adrenal_8wpf_1_Adrenal_trimmed": "GSM5815990_Hs46_Adrenal",
             "GSM8831102_fetal_adrenal_8wpf_2": "GSM5815989_HS12H_adrenal"}
ROOST = {"CH1 adrenal_9": 9.4, "CC2 adrenal_18": 18, "CF2 adrenal_22": 22}
KILDISIUTE = {"w8": 8, "w8d6": 8 + 6 / 7, "w10d5": 10 + 5 / 7, "w11": 11, "w21": 21}
SCHRODER_BATCH = {**{f"FA{i}": "148" for i in range(1, 5)}, **{f"FA{i}": "150" for i in range(5, 12)},
                  "FA12": "153", "FA13": "153"}
ENS = {}  # Ensembl gene id (no version) -> symbol, filled from the sources that give both


def by_symbol(values, symbols) -> pd.Series:
    """Counts summed per canonical symbol."""
    s = pd.Series(np.asarray(values, float), index=[canonical(g) for g in symbols])
    return s.groupby(level=0).sum()


def tenx(h, min_umi: int = 0):
    """(counts per symbol, cells, barcodes) of a 10x h5 (file or file object)."""
    with h5py.File(h, "r") as f:
        m = f["matrix"]
        x = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
        names, ids = m["features/name"][:].astype(str), m["features/id"][:].astype(str)
        barcodes = m["barcodes"][:].astype(str)
    ENS.update(zip(ids, names))
    keep = np.asarray(x.sum(0)).ravel() >= max(min_umi, 1)
    return by_symbol(np.asarray(x[:, keep].sum(1)).ravel(), names), int(keep.sum()), set(barcodes[keep])


def wide(h, symbol_col: str | None = None):
    """(counts per symbol, cells) of a gene x cell text table."""
    t = pd.read_csv(h, sep="\t", index_col=0)
    if symbol_col:
        ENS.update(zip(t.index, t[symbol_col]))
        t = t.set_index(symbol_col)
    return by_symbol(t.sum(axis=1).to_numpy(), t.index), t.shape[1]


def encode_long(path):
    t = pd.read_csv(path, sep="\t")
    t = t[t.gene_novelty == "Known"]
    ENS.update(zip(t.annot_gene_id.str.split(".").str[0], t.annot_gene_name))
    return by_symbol(t.iloc[:, -1].to_numpy(), t.annot_gene_name)


def from_ensembl(values, ids) -> pd.Series:
    sym = pd.Series([ENS.get(str(i).split(".")[0]) for i in ids])
    ok = sym.notna().to_numpy()
    return by_symbol(np.asarray(values, float)[ok], sym[ok])


def main():
    DEST.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)
    rows, pb = [], {}

    def add(cohort, source, assay, fetus, age, convention, counts, cells=np.nan, sample=""):
        key = f"{cohort}:{fetus}"
        pb[key] = counts if key not in pb else pb[key].add(counts, fill_value=0)
        rows.append(dict(id=key, cohort=cohort, source=source, assay=assay, fetus=str(fetus), age=float(age),
                         convention=convention, sample=sample, n_cells=cells, depth=float(counts.sum())))
        print(f"{key:28s} {sample[:40]:40s} cells {cells} depth {counts.sum():.3g}", flush=True)

    enc = SRC / "encode_fetal_adrenal" / "raw"
    for fetus, days, f in ENCODE_LR:
        add("encode_lr", "ENCODE", "long-read RNA-seq", fetus, days / 7, "post conception", encode_long(enc / f"{f}.tsv"),
            sample=f)
    with tarfile.open(SRC / "kameneva2021" / "raw" / "GSE147821_RAW.tar") as tar:
        for name, age in KAMENEVA.items():
            c, n, _ = tenx(tar.extractfile(f"{name}.raw_feature_bc_matrix.h5"), MIN_UMI)
            add("kameneva", "GSE147821", "10x 3'", name.split("_", 1)[1], age, "post conception", c, n, name)
    with tarfile.open(SRC / "dong2020" / "raw" / "GSE137804_RAW.tar") as tar:
        for name, age in DONG.items():
            with gzip.open(tar.extractfile(f"{name}_gene_cell_exprs_table.xls.gz"), "rt") as h:
                c, n = wide(h, "Symbol")
            add("dong", "GSE137804", "10x 3'", name.split("_")[1], age, "post conception", c, n, name)
    for name, (fetus, age) in WANG.items():
        c, n, _ = tenx(SRC / "wang2022" / "raw" / f"{name}_filtered_feature_bc_matrix.h5")
        add("wang", "GSE167860", "10x 3'", fetus, age, "gestational", c, n, name)
    for name, age in HAN.items():
        c, n = wide(SRC / "han2020" / "raw" / f"{name}_dge.txt.gz")
        add("han", "GSE134355", "Microwell-seq", f"{age}w", age, "not stated", c, n, name)
    barcodes = {}
    for name, age in SASAKI.items():
        c, n, barcodes[name.split("/")[-1]] = tenx(SRC / f"{name}_filtered_feature_bc_matrix.h5")
        add("sasaki", name.split("/")[0], "10x 3'", f"{age}wpf", age, "post fertilisation", c, n, name.split("/")[-1])
    for dup, orig in DUPLICATE.items():
        b = tenx(SRC / "gse291218" / "raw" / f"{dup}_filtered_feature_bc_matrix.h5")[2]
        share = len(b & barcodes[orig]) / min(len(b), len(barcodes[orig]))
        assert share > 0.95, (dup, orig, share)
        print(f"duplicate {dup} = {orig}: shared barcodes {share:.3f}", flush=True)
    for cohort, fetus, days, f in ENCODE_SR:  # after the 10x sources: needs the Ensembl -> symbol map
        t = pd.read_csv(enc / f"{f}.tsv", sep="\t")
        add(cohort, "ENCODE", "polyA plus RNA-seq" if cohort == "encode_pa" else "total RNA-seq", fetus, days / 7,
            "post conception", from_ensembl(t.expected_count, t.gene_id), sample=f)
    t = pd.read_csv(SRC / "roost2015" / "raw" / "GSE66302_Roost_et_al_Raw_Reads.txt.gz", sep="\t", index_col=0)
    for name, age in ROOST.items():
        add("roost", "GSE66302", "DeepSAGE", name.split()[0], age, "gestational", from_ensembl(t[name], t.index),
            sample=name)
    a = ad.read_h5ad(SRC / "kildisiute2021" / "raw" / "adr_all_cellxgene.h5ad")
    x = sp.csr_matrix(a.X).expm1()
    x = sp.diags(a.obs.nCount_RNA.to_numpy(float) / np.asarray(x.sum(1)).ravel()) @ x
    for age_label, g in a.obs.reset_index(drop=True).groupby(a.obs.GestationalAge.astype(str).to_numpy()):
        add("kildisiute", "kildisiute2021", "10x 3' (SoupX-corrected)", age_label, KILDISIUTE[age_label],
            "post conception", by_symbol(np.asarray(x[g.index.to_numpy()].sum(0)).ravel(), a.var_names), len(g),
            "+".join(sorted(g.SampleName.astype(str).unique())))
    t = pd.read_csv(SRC / "schroder2022" / "raw" / "GSE196700_Read_counts_runs_combined.txt.gz", sep="\t")
    for name, batch in SCHRODER_BATCH.items():
        add("schroder", "GSE196700", "targeted transcriptome (amplicons)", name, np.nan, f"batch {batch}",
            by_symbol(t[name], t.Gene), sample=name)

    s = pd.DataFrame(rows)
    sheet = s.groupby("id", sort=False).agg(
        cohort=("cohort", "first"), source=("source", lambda v: "+".join(sorted(set(v)))), assay=("assay", "first"),
        fetus=("fetus", "first"), age=("age", "first"), convention=("convention", "first"),
        n_samples=("sample", "size"), samples=("sample", "+".join), n_cells=("n_cells", "sum"),
        depth=("depth", "sum")).reset_index()
    counts = pd.DataFrame(pb).T.loc[sheet.id]
    counts.columns = counts.columns.astype(str)
    counts.to_parquet(DEST / "pseudobulk.parquet")
    sheet.to_csv(DEST / "samples.csv", index=False)
    sheet["window"] = np.where(sheet.age.isna(), "unknown", np.where(sheet.age >= 12, "from 12 weeks", "before 12 weeks"))
    design = sheet.groupby(["cohort", "window"], sort=False).agg(
        assay=("assay", "first"), convention=("convention", "first"), n_fetuses=("id", "size"),
        ages=("age", lambda v: " ".join(f"{x:.1f}" for x in sorted(v.dropna()))),
        median_cells=("n_cells", "median"), median_depth=("depth", "median")).reset_index()
    design.to_csv(OUT / "design_cohorts.csv", index=False)
    det = counts[GENES].groupby(sheet.cohort.to_numpy(), sort=False).median().round(1)
    det.to_csv(OUT / "clock_detection.csv")
    pd.set_option("display.width", 250)
    print(design.to_string(index=False), det.to_string(), f"symbols {counts.shape[1]}", sep="\n", flush=True)


if __name__ == "__main__":
    main()
