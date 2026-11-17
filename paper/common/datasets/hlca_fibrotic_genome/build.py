"""Lung, HLCA full, healthy and fibrotic donors, all genes: donor x compartment pseudobulks of every gene (A286, lung
fibrosis pairwise contrast). Raw input $SCCC_DATA/atlas_wg/lung_hlca. The rules of `lung_hlca_genome` (classes, cell and
donor thresholds, strata) with the fibrotic donors of `lung.FIBROTIC` added; donor columns `healthy` and `hlca_fibrotic`
(control / case). The clock-gene dataset `hlca_fibrotic` uses another raw build (hlca_tcell), so it is not an exact twin;
the healthy rows are checked against `lung_hlca_genome`. Choices: see `_atlas/genome.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import genome as M  # noqa: E402

DATASET = "hlca_fibrotic_genome"


def build(donors=None, genes="all"):
    """(pb, donors) of every gene from the public genome-wide group pseudobulk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)
