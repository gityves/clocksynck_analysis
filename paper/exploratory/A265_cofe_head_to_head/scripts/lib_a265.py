"""Data of A265 (COFE against clocksync on the TCGA tumour-normal pairs), built as COFE's TCGA analysis (Ananthasubramaniam
and Venkataramanan 2025, PLoS Biol 23:e3003196, S4 Text and S2 Table).

Per cancer: primary tumours of the histological type of S2 Table, one per patient (lowest barcode); recount3 scaled counts
(`tcga_paired.counts_of`); protein-coding genes (GENCODE v29) with a mean count >= 64 over the tumours; log2 cpm (library
size = sum of the protein-coding counts, pseudocount 1: S4 Text names DESeq2 without the call, see the registration);
then genes with mean log2 cpm >= 7. Solid-tissue normals: one per patient whose tumour is in that set (`tcga_paired`
rule), on the same genes and the same transform, never used to train COFE.
"""
from __future__ import annotations

import gzip
import sys
from functools import cache
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import data  # noqa: E402

B = data._builder("tcga_paired")
CANCERS = B.CANCERS  # the ten with normals; COFE's eleventh (OV) has none
HISTOLOGY = {"BRCA": "Infiltrating Ductal Carcinoma", "KIRC": "Kidney Clear Cell Renal Carcinoma",
             "PRAD": "Prostate Adenocarcinoma Acinar Type", "BLCA": "Muscle invasive urothelial carcinoma (pT2 or above)",
             "UCEC": "Endometrioid endometrial adenocarcinoma", "COAD": "Colon Adenocarcinoma",
             "THCA": "Thyroid Papillary Carcinoma - Classical/usual", "LIHC": "Hepatocellular Carcinoma",
             "LUAD": "Lung Adenocarcinoma- Not Otherwise Specified (NOS)", "KIRP": "Kidney Papillary Renal Cell Carcinoma"}
MIN_COUNT, MIN_LOG2CPM = 64.0, 7.0


@cache
def protein_coding() -> frozenset:
    names = set()
    with gzip.open(B.RECOUNT3 / "annotation" / "human.gene_sums.G029.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.split("\t", 8)
            if col[2] == "gene" and 'gene_type "protein_coding"' in col[8]:
                names.add(col[8].split('gene_name "', 1)[1].split('"', 1)[0])
    return frozenset(names)


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    return np.log2(counts.div(counts.sum(1), axis=0) * 1e6 + 1.0)


def cancer(c: str) -> dict:
    """tumour / normal counts (protein-coding genes), the COFE matrices (log2 cpm, filtered genes) and the pairs."""
    s = B.samples_of(c)
    t = s[(s.sample_type == "Primary Tumor") & (s.histology == HISTOLOGY[c])].sort_values("barcode").drop_duplicates("patient")
    n = s[s.sample_type == "Solid Tissue Normal"].sort_values("barcode").drop_duplicates("patient")
    n = n[n.patient.isin(set(t.patient))]
    ct, cn = B.counts_of(c, t), B.counts_of(c, n)
    dt, dn = ct.sum(1).to_numpy(float), cn.sum(1).to_numpy(float)  # depth over all genes (as `tcga_paired`)
    pc = [g for g in ct.columns if g in protein_coding()]
    ct, cn = ct[pc], cn[pc]
    keep = ct.columns[ct.mean() >= MIN_COUNT]
    lt, ln = log2cpm(ct[keep]), log2cpm(cn[keep])
    genes = lt.columns[lt.mean() >= MIN_LOG2CPM]
    lt.index, ln.index = t.patient.to_numpy(), n.patient.to_numpy()
    ct.index, cn.index = t.patient.to_numpy(), n.patient.to_numpy()
    return dict(counts_t=ct, counts_n=cn, depth_t=pd.Series(dt, index=ct.index), depth_n=pd.Series(dn, index=cn.index),
                cofe_t=lt[genes], cofe_n=ln[genes], pairs=list(n.patient))


def cofe_matrices(ct: pd.DataFrame, cn: pd.DataFrame):
    """COFE's filters and transform on (possibly simulated) protein-coding counts: (tumour, normal) log2 cpm."""
    keep = ct.columns[ct.mean() >= MIN_COUNT]
    lt, ln = log2cpm(ct[keep]), log2cpm(cn[keep])
    genes = lt.columns[lt.mean() >= MIN_LOG2CPM]
    return lt[genes], ln[genes]
