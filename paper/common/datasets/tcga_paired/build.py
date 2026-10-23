"""TCGA primary tumour vs solid-tissue normal of the same patient (ten cancers), recount3 gene counts.

Source: `$SCCC_DATA/tcga/recount3/<project>/tcga.gene_sums.<project>.G029.gz` with the recount3 metadata and the
GENCODE v29 annotation (recount3, human G029; PROVENANCE in `tcga/recount3/`).
Legacy builder replaced: `paired_damping/scripts/a52_tcga.py:pairs` with `sccc.benchmark.data.tcga_samples/tcga_counts`.
Used by: A52 (and a52_secondary), A80, A83, A84, A87, A97 (paired tumour / normal cohorts), A88.

Choices (as legacy): per patient one primary tumour and one solid-tissue normal (lowest barcode when several), kept
for patients with both; counts scaled as recount3 `transform_counts` does by default, round(coverage x 4e7 / AUC)
(`bc_auc.all_reads_all_bases`); gene symbols from GENCODE v29 with duplicated symbols summed; depth = the sum of
the scaled counts over all genes. unit = patient (12-character barcode), condition = tumour / normal, cohort = cancer,
`barcode` and the clinical fields (histology, pathologic stage, age at diagnosis, vital status, days to death / last
follow-up) are kept.
`gene_means()` returns the legacy expression-matched control input (mean log fraction of the normals, genes with a
mean count >= 1). Unpaired tumours (A83, A97 all primary tumours) are not in this table (`tumours()` below).
"""

from __future__ import annotations

import gzip
from functools import cache

import numpy as np
import pandas as pd

import sys
from pathlib import Path

from common.paths import DATA

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_bulk_common"))
import bulk_tables as T  # noqa: E402

RECOUNT3 = DATA / "tcga" / "recount3"
CANCERS = ["BLCA", "BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "PRAD", "THCA", "UCEC"]
LEVELS = ("normal", "tumour")
CLINICAL = ["histology", "stage", "age", "vital", "t_death", "t_fu"]  # patient-level TCGA clinical fields


@cache
def gene_names() -> pd.Series:
    """GENCODE v29 gene_id -> gene_name (recount3 annotation)."""
    ids, names = [], []
    with gzip.open(RECOUNT3 / "annotation" / "human.gene_sums.G029.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.split("\t", 8)
            if col[2] != "gene":
                continue
            attr = dict(a.strip().split(" ", 1) for a in col[8].strip().strip(";").split(";") if a.strip())
            ids.append(attr["gene_id"].strip('"'))
            names.append(attr["gene_name"].strip('"'))
    return pd.Series(names, index=ids)


def _md(project: str, table: str) -> pd.DataFrame:
    return pd.read_csv(RECOUNT3 / project / f"tcga.{table}.{project}.MD.gz", sep="\t", low_memory=False)


def samples_of(project: str) -> pd.DataFrame:
    """One row per sample: external_id, barcode, patient, sample_type, histology, auc."""
    md = _md(project, "tcga").merge(_md(project, "recount_qc")[["external_id", "bc_auc.all_reads_all_bases"]],
                                    on="external_id")
    clin = {"age": "cgc_case_age_at_diagnosis", "stage": "cgc_case_pathologic_stage", "vital": "cgc_case_vital_status",
            "t_death": "cgc_case_days_to_death", "t_fu": "cgc_case_days_to_last_follow_up"}
    return pd.DataFrame({"external_id": md.external_id, "barcode": md.tcga_barcode,
                         "patient": md.tcga_barcode.str[:12], "sample_type": md.cgc_sample_sample_type,
                         "histology": md.cgc_case_histological_diagnosis, "auc": md["bc_auc.all_reads_all_bases"],
                         **{k: md[v] if v in md else np.nan for k, v in clin.items()}})


def counts_of(project: str, samples: pd.DataFrame) -> pd.DataFrame:
    """Scaled counts (samples x gene symbols, duplicates summed) of the rows of `samples`, index = barcode."""
    ids = set(samples.external_id)
    cov = pd.read_csv(RECOUNT3 / project / f"tcga.gene_sums.{project}.G029.gz", sep="\t", comment="#",
                      index_col=0, usecols=lambda c: c == "gene_id" or c in ids)[samples.external_id]
    counts = np.round(cov * (4e7 / samples.auc.to_numpy()))
    counts.index = gene_names().reindex(counts.index).to_numpy()
    counts = counts.groupby(level=0).sum().T
    counts.index = samples.barcode.to_numpy()
    return counts


def _pairs(project: str):
    s = samples_of(project)
    t = s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
    n = s[s.sample_type == "Solid Tissue Normal"].sort_values("barcode").drop_duplicates("patient")
    pts = sorted(set(t.patient) & set(n.patient))
    t = t.set_index("patient").loc[pts].reset_index()
    n = n.set_index("patient").loc[pts].reset_index()
    return t, n, counts_of(project, t), counts_of(project, n)


def _table(project, meta, counts, condition):
    d = pd.DataFrame({"unit": meta.patient.to_numpy(), "condition": condition, "cohort": project,
                      "value_type": "counts", "depth": counts.sum(1).to_numpy(), "barcode": meta.barcode.to_numpy(),
                      **{c: meta[c].to_numpy() for c in CLINICAL}})
    cols = T.clock_columns(counts.columns)
    return pd.concat([d, pd.DataFrame(counts[cols].to_numpy(float), columns=T.GENES)], axis=1)


def build(cancers=CANCERS) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(samples, units): tumour and normal rows of the paired patients of each cancer."""
    parts = []
    for c in cancers:
        t, n, ct, cn = _pairs(c)
        parts += [_table(c, t, ct, "tumour"), _table(c, n, cn, "normal")]
    samples = pd.concat(parts, ignore_index=True)
    return samples, T.units_table(samples, ("histology", "stage", "age"))


def gene_means(cancers=CANCERS) -> dict[str, pd.Series]:
    """cancer -> mean log fraction of all genes in the paired normals (genes with mean count >= 1)."""
    out = {}
    for c in cancers:
        _, _, _, cn = _pairs(c)
        lm = np.log(cn.div(cn.sum(1), axis=0).mean() + 1e-12)
        out[c] = lm[cn.mean() >= 1]
    return out


def tumours(project: str) -> pd.DataFrame:
    """All primary tumours (one per patient, lowest barcode, any histology) of a project incl. OV, unpaired."""
    s = samples_of(project)
    t = s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
    return _table(project, t, counts_of(project, t), "tumour")
