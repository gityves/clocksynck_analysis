"""Pooled normal adult lung atlases of the gene-rhythm map (HLCA healthy donors + 21 Census datasets), all genes: donor x lung class
pseudobulks for `Dataset.gene_rhythms` (A50, replication of A31b). Raw inputs $SCCC_DATA/gene_map/lung and atlas_wg/lung_hlca.
Legacy builder replaced: desync/validation/gene_map/scripts/pool.py (`prep lung`). Choices: see `_atlas/genemap.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import genemap as M  # noqa: E402

DATASET = "gene_map_lung"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
