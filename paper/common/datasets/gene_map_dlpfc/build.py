"""Pooled normal adult DLPFC atlases of the gene-rhythm map (PsychAD cohorts, SEA-AD and other Census datasets), all genes: donor x cortical
class pseudobulks for `Dataset.gene_rhythms` (A50). Raw input $SCCC_DATA/gene_map/dlpfc.
Legacy builder replaced: desync/validation/gene_map/scripts/pool.py (`prep dlpfc`). Choices: see `_atlas/genemap.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import genemap as M  # noqa: E402

DATASET = "gene_map_dlpfc"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
