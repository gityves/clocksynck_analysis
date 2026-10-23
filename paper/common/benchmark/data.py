"""Bulk expression for the cyclic-ordering benchmark: GTEx tissues with recorded time of death, and TCGA tumours
(recount3) as analysed by COFE (Ananthasubramaniam et al. 2025, S4 Text).

The same preprocessing is applied to both, following COFE's TCGA analysis: protein-coding genes (GENCODE v29) with
a mean of at least 64 reads per sample, log2 counts per million, then genes with a mean log2 CPM of at least 7.
"""

from __future__ import annotations

import gzip
from functools import cache

import anndata as ad
import numpy as np
import pandas as pd

from common.paths import DATA

GTEX_LABELLED = "/home/maxine/Documents/paychere/scRNAseq_database/phase_inference/baselines/data/gtex/gtex_v10_labelled.h5ad"
RECOUNT3 = DATA / "tcga" / "recount3"
WINDOW = (8.0, 18.0)  # hospital hours

# COFE's 11 adenocarcinomas and the histological type kept for each (S2 Table).
TCGA_TYPES = {
    "BRCA": "Infiltrating Ductal Carcinoma",
    "KIRC": "Kidney Clear Cell Renal Carcinoma",
    "PRAD": "Prostate Adenocarcinoma Acinar Type",
    "OV": "Serous Cystadenocarcinoma",
    "BLCA": "Muscle invasive urothelial carcinoma (pT2 or above)",
    "UCEC": "Endometrioid endometrial adenocarcinoma",
    "COAD": "Colon Adenocarcinoma",
    "THCA": "Thyroid Papillary Carcinoma - Classical/usual",
    "LIHC": "Hepatocellular Carcinoma",
    "LUAD": "Lung Adenocarcinoma- Not Otherwise Specified (NOS)",
    "KIRP": "Kidney Papillary Renal Cell Carcinoma",
}


@cache
def gencode_genes() -> pd.DataFrame:
    """GENCODE v29 genes of recount3: gene_id (versioned), gene_name, gene_type."""
    rows = []
    with gzip.open(RECOUNT3 / "annotation" / "human.gene_sums.G029.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            if col[2] != "gene":
                continue
            attr = dict(a.strip().split(" ", 1) for a in col[8].strip(";").split(";") if a.strip())
            rows.append({k: attr[k].strip('"') for k in ("gene_id", "gene_name", "gene_type")})
    return pd.DataFrame(rows).set_index("gene_id")


def log2_normalised(counts: pd.DataFrame, depth: float = 4e7) -> pd.DataFrame:
    """log2(1 + counts / size factor), DESeq2 median-of-ratios size factors scaled so that the mean normalised
    library is `depth` reads (recount3's TCGA scaling), which makes the COFE thresholds comparable across datasets."""
    pos = counts.loc[:, (counts > 0).all()]
    logc = np.log(pos)
    sf = np.exp((logc.sub(logc.mean(), axis=1)).median(axis=1))
    norm = counts.div(sf, axis=0)
    return np.log2(1 + norm * (depth / norm.sum(1).mean()))


def cofe_filter(counts: pd.DataFrame, min_reads: float = 64, min_log2: float = 7) -> pd.DataFrame:
    """Samples x genes read counts (gene symbols as columns) -> log2 normalised expression of the protein-coding
    genes with a mean of `min_reads` reads and a mean log2 expression of `min_log2`. COFE's S4 Text calls the
    second scale "log2 cpm"; its thresholds and its ~13,000 expressed genes per cancer only fit normalised counts
    (LIHC: 11,570 genes; with counts per million, 830)."""
    coding = set(gencode_genes().query("gene_type == 'protein_coding'").gene_name)
    counts = counts.loc[:, [g in coding for g in counts.columns]]
    counts = counts.loc[:, counts.mean() >= min_reads]
    lg = log2_normalised(counts)
    return lg.loc[:, lg.mean() >= min_log2]


# --- GTEx ---------------------------------------------------------------------------------------------------

@cache
def _gtex() -> ad.AnnData:
    a = ad.read_h5ad(GTEX_LABELLED, backed="r")
    return a[a.obs.DTHHRDY.isin([1, 2]).to_numpy()]  # fast deaths: time of death reflects the living clock


def gtex_tissue(tissue: str) -> tuple[pd.DataFrame, pd.Series]:
    """Read counts (donors x gene symbols, duplicates summed) and time of death (h, 0-24) of the fast deaths."""
    a = _gtex()
    a = a[(a.obs.tissue == tissue).to_numpy()].to_memory()
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    donors = a.obs.subject.astype(str).to_numpy()
    counts = pd.DataFrame(X, index=donors, columns=a.var_names).T.groupby(level=0).sum().T
    return counts, pd.Series(a.obs.tod.to_numpy(float) % 24, index=donors, name="tod")


def in_window(hours, window=WINDOW) -> np.ndarray:
    h = np.asarray(hours) % 24
    return (h >= window[0]) & (h < window[1])


# --- TCGA (recount3) ----------------------------------------------------------------------------------------

def _metadata(project: str, table: str) -> pd.DataFrame:
    return pd.read_csv(RECOUNT3 / project / f"tcga.{table}.{project}.MD.gz", sep="\t", low_memory=False)


def tcga_samples(project: str) -> pd.DataFrame:
    """One row per sample: external_id, barcode, patient, sample type, histology, AUC (bc_auc.all_reads_all_bases)."""
    md = _metadata(project, "tcga")
    qc = _metadata(project, "recount_qc")[["external_id", "bc_auc.all_reads_all_bases"]]
    md = md.merge(qc, on="external_id")
    return pd.DataFrame({
        "external_id": md.external_id, "barcode": md.tcga_barcode, "patient": md.tcga_barcode.str[:12],
        "sample_type": md.cgc_sample_sample_type, "histology": md.cgc_case_histological_diagnosis,
        "auc": md["bc_auc.all_reads_all_bases"],
    })


def tcga_counts(project: str, samples: pd.DataFrame) -> pd.DataFrame:
    """Read counts (samples x gene symbols, duplicates summed) for `samples` (rows of `tcga_samples`), scaled as
    recount3::transform_counts does by default: round(coverage * 4e7 / AUC)."""
    path = RECOUNT3 / project / f"tcga.gene_sums.{project}.G029.gz"
    cov = pd.read_csv(path, sep="\t", comment="#", index_col=0, usecols=lambda c: c == "gene_id" or c in set(samples.external_id))
    cov = cov[samples.external_id]
    counts = np.round(cov * (4e7 / samples.auc.to_numpy()))
    counts.index = gencode_genes().gene_name.reindex(counts.index).to_numpy()
    counts = counts.groupby(level=0).sum().T
    counts.index = samples.barcode.to_numpy()
    return counts


def tcga_tumours(project: str) -> pd.DataFrame:
    """COFE's selection: primary tumours of the project's histological type, one per patient (lowest barcode)."""
    s = tcga_samples(project)
    s = s[(s.sample_type == "Primary Tumor") & (s.histology == TCGA_TYPES[project])]
    return s.sort_values("barcode").drop_duplicates("patient")


def tcga_normals(project: str) -> pd.DataFrame:
    """Solid tissue normal samples, one per patient."""
    s = tcga_samples(project)
    return s[s.sample_type == "Solid Tissue Normal"].sort_values("barcode").drop_duplicates("patient")
