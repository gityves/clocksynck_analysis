"""Lung ILD of Natri et al. 2024 (control and fibrotic donors), all genes: donor x lung class pseudobulks of every gene
(A286, lung fibrosis pairwise contrast). Raw input $SCCC_DATA/atlas_wg/natri2024_ild; the rules of `lung_natri` (classes
from annotation_level_1 and fibroblast lineages, >= 20 cells, stratum = majority processing site; donor column
`ild_natri` control / case, sarcoidosis excluded). Twin of `lung_natri` (clock-gene columns checked against it).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import lung as M  # noqa: E402

DATASET = "lung_natri"


def build(donors=None, genes="all"):
    """(pb, donors) of every gene from the public genome-wide group pseudobulk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)
