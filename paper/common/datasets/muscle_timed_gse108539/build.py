"""Timed human skeletal muscle, Geneva hourly-feeding cohort: 10 volunteers x 6 biopsies, 4-hourly over 24 h (GSE108539).

Source: `$SCCC_DATA/muscle_timed/GSE108539/` (Perrin et al. 2018, eLife 7:e34114; same biopsies as Loizides-Mangold 2017).
Vastus lateralis, total RNA-seq (HiSeq, hg38). Controlled laboratory routine, hourly isocaloric meal replacement while awake
(lights on 07:00-22:00). The processed matrix is log2 edgeR RPKM per gene, exon (mRNA) and intron (pre-mRNA) layers, Ensembl ids:
no counts, so the table is an array-type table (`value_type` = log2_rpkm), as `blood_timed_gse39445`. ARNTL is absent from the
filtered matrix (13 of 14 panel genes; column left NaN).
Used by: A113 (WP4; timed validation, estimators of the registration), scout of A130.

Table: `unit` = library, `block` = subject (GEO number; the unit of within-subject centring), `clock_h` (clock time of the
biopsy, h, 12 16 20 0 4 8), `order` (biopsy index 0-5 from 12:00), `layer` (exon | intron), `value_type`, the 14 genes.
Biopsy times are recorded clock times of a laboratory routine, not hospital times. Donors table (index = subject): `n_biopsies`.
Ages and sex are public only as a cohort aggregate (supplement Table S1: 9 M / 1 F, 22-54 y); the GEO subject numbers cannot be
linked to it, so there is no per-subject age.
"""

from __future__ import annotations

import re

import pandas as pd

import clocksync as cs
from common.paths import DATA

D = DATA / "muscle_timed" / "GSE108539"
GENES_REF = DATA / "reference" / "gene_positions_gencode46.tsv"  # Ensembl id, symbol, ...


def _symbols() -> pd.Series:
    g = pd.read_csv(GENES_REF, sep="\t", header=None, names=["id", "symbol", "chr", "start", "end", "type"])
    return g.assign(id=g.id.str.split(".").str[0]).drop_duplicates("id").set_index("id").symbol


def build(layer: str = "exon") -> tuple[pd.DataFrame, pd.DataFrame]:
    """(libraries, subjects): one row per biopsy with the 14 genes (log2 RPKM of `layer`), and the subject table."""
    x = pd.read_csv(D / "GSE108539_norm_reads_count.txt.gz", sep="\t", index_col=0)
    x.index = x.index.str.split(".").str[0]
    sym = _symbols().reindex(x.index)
    genes = list(cs.gtex().genes)
    x = x[sym.isin(genes).to_numpy()].set_axis(sym[sym.isin(genes)].to_numpy(), axis=0)
    x = x.loc[:, [c.lower().startswith(layer) for c in x.columns]].reindex(genes)
    rows = [re.match(r"(?:Exon|Intron)_(\d+)_S(\d+)_(\d+)\.\d+$", c).groups() for c in x.columns]
    t = pd.DataFrame(dict(unit=x.columns, block=[r[0] for r in rows], clock_h=[float(r[2]) for r in rows],
                          layer=layer, value_type="log2_rpkm", depth=float("nan")))
    t["order"] = ((t.clock_h - 12) % 24 / 4).round().astype(int)
    t = pd.concat([t, x.T.reset_index(drop=True)], axis=1)
    subjects = t.groupby("block").size().rename("n_biopsies").to_frame()
    return t, subjects
