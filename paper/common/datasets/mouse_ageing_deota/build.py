"""Mouse restricted feeding with known time: 22 tissues, ad libitum vs 9 h time-restricted feeding (Deota et al. 2023).

Source: `$SCCC_DATA/mouse_ageing/Deota2023/GSE190389_geo/GSE190389_raw.csv.gz` (GSE190389, raw counts, 2 mice per arm
and ZT, ZT 0-22, latin-1 encoded; tissues of the same mouse are inferred from the library names).
Legacy builder replaced: `desync/validation/mouse_ageing/scripts/a49_deota.py:load` with `common.clock_table`.
Used by: A49 (H49.1, H49.2).

Choices (as legacy): library names `<tissue>_<ZT>_<arm>_<mouse>_<run>`; the Excel-mangled space of "Kidney cortex"
is repaired; runs of the same tissue x mouse are averaged and rounded half to even (np.rint); donor =
`{arm}_ZT{zt}_m{mouse}` (one mouse, all its tissues share the time), compartment = tissue, depth = sum of all
(rounded) counts; stratum = arm (as A49 sets it). Tissue QC and the unit rule (`drop_units`: >= 3 tissues and >= 2 brain
anchors) are analysis steps, provided as functions. Donors table: `zt`, `arm`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import clocksync as cs
from common.paths import DATA

SRC = DATA / "mouse_ageing" / "Deota2023" / "GSE190389_geo" / "GSE190389_raw.csv.gz"
BRAIN = ["SCN", "Arcuate", "DMH", "Hypothalamus", "DG", "Amygdala"]  # anchors of A49
MODULE = ["Liver", "eWAT", "iWAT", "BAT", "Pancreas", "Muscle", "Stomach", "Duodenum", "Jejunum", "Ileum"]
COMPANION = ["Heart", "Lung", "Spleen", "Kidney cortex", "Kidney medulla"]
EXPR_MIN = 5.0


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    b = pd.read_csv(SRC, index_col=0, encoding="latin-1")
    b.columns = [c.replace("æ", " ") for c in b.columns]
    rec = []
    for c in b.columns:
        tis, zt, arm, mouse, run = c.rsplit("_", 4)
        rec.append(dict(col=c, compartment=tis, zt=int(zt[2:]), arm=arm, donor=f"{arm}_ZT{int(zt[2:])}_m{mouse}"))
    m = pd.DataFrame(rec)
    cols, rows = [], []
    for _, g in m.groupby(["compartment", "donor"]):
        cols.append(np.rint(b[g.col].mean(1).to_numpy()).astype(np.int64))
        rows.append(g.iloc[0][["compartment", "zt", "arm", "donor"]])
    cnt = pd.DataFrame(np.array(cols), columns=b.index)
    meta = pd.DataFrame(rows).reset_index(drop=True)
    meta["stratum"] = meta.arm
    pb = cs.from_bulk(cnt, meta, donor="donor", compartment="compartment", stratum="stratum")
    donors = meta.drop_duplicates("donor").set_index("donor")[["zt", "arm"]]
    return pb, donors


def tissue_qc(pb: pd.DataFrame) -> pd.DataFrame:
    med = pb.groupby("compartment")[list(cs.gtex().genes)].median()
    q = pd.DataFrame(dict(n_units=pb.groupby("compartment").donor.nunique(), min_gene_median=med.min(1),
                          weakest_gene=med.idxmin(1)))
    q["kept"] = q.min_gene_median >= EXPR_MIN
    return q


def drop_units(pb: pd.DataFrame, anchors, min_tissues: int = 3, min_anchors: int = 2):
    n = pb.groupby("donor").compartment.nunique()
    na = pb[pb.compartment.isin(anchors)].groupby("donor").compartment.nunique().reindex(n.index).fillna(0)
    bad = list(n.index[(n < min_tissues) | (na < min_anchors)])
    return pb[~pb.donor.isin(bad)].reset_index(drop=True), bad
