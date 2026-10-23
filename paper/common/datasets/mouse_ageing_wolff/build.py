"""Mouse ageing with known circadian time: 6 tissues x 3 ages x 12 circadian times (Wolff et al. 2023, GSE201207).

Source: `$SCCC_DATA/mouse_ageing/GSE201207/GSE201207_count_216.csv.gz` (raw counts, mouse symbols, 216 libraries =
pools of 2 mice; tissues adrenal, heart, hypothalamus, kidney, lung, muscle; ages young / aged / old; CT 18-62 h).
Legacy builder replaced: `desync/validation/mouse_ageing/scripts/a48_wolff.py:load` with `common.clock_table`.
Used by: A48 (H48.1-H48.3), the mouse ageing panels.

Choices (as legacy): unit (donor) = age x circadian time pool, `{age}_ct{CT}` (all tissues of a pool share the time),
compartment = tissue, depth = sum of all counts of the library, clock genes by mouse symbol mapped to the panel.
Tissue QC (`tissue_qc`: median count of every clock gene >= 5) and the unit rule
(`drop_units`: at least 3 tissues and 2 anchors per unit) are analysis steps and are provided as functions, not
applied. Donors table (index = unit): `age`, `ct` (h).
"""

from __future__ import annotations

import pandas as pd

import clocksync as cs
from common.paths import DATA

SRC = DATA / "mouse_ageing" / "GSE201207" / "GSE201207_count_216.csv.gz"
PERIPHERAL = ["adrenal", "heart", "kidney", "lung", "muscle"]  # anchors of A48 (hypothalamus is the tested clock)
AGES = ["young", "aged", "old"]
EXPR_MIN = 5.0


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    c = pd.read_csv(SRC, index_col=0)
    rec = []
    for s in c.columns:
        ct, age, tis = s.split("-")
        rec.append(dict(lib=s, donor=f"{age}_ct{int(ct[2:])}", compartment=tis, age=age, ct=int(ct[2:])))
    meta = pd.DataFrame(rec).set_index("lib")
    pb = cs.from_bulk(c.T, meta, donor="donor", compartment="compartment")
    donors = meta.drop_duplicates("donor").set_index("donor")[["age", "ct"]].sort_values(["age", "ct"])
    return pb, donors


def tissue_qc(pb: pd.DataFrame) -> pd.DataFrame:
    """Median count of the weakest clock gene per tissue; `kept` when >= 5 (legacy `common.tissue_qc`)."""
    med = pb.groupby("compartment")[list(cs.gtex().genes)].median()
    q = pd.DataFrame(dict(n_units=pb.groupby("compartment").donor.nunique(), min_gene_median=med.min(1),
                          weakest_gene=med.idxmin(1)))
    q["kept"] = q.min_gene_median >= EXPR_MIN
    return q


def drop_units(pb: pd.DataFrame, anchors, min_tissues: int = 3, min_anchors: int = 2):
    """Units with fewer than `min_tissues` tissues or `min_anchors` anchor tissues removed; returns (pb, dropped)."""
    n = pb.groupby("donor").compartment.nunique()
    na = pb[pb.compartment.isin(anchors)].groupby("donor").compartment.nunique().reindex(n.index).fillna(0)
    bad = list(n.index[(n < min_tissues) | (na < min_anchors)])
    return pb[~pb.donor.isin(bad)].reset_index(drop=True), bad
