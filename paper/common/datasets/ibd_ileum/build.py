"""IBDverse ileum: healthy donors and Crohn's disease inflamed / uninflamed ileal biopsies (scRNA-seq): donor x class
pseudobulks. Raw input $SCCC_DATA/geo/IBDverse/pseudobulk (group level; raw h5ad in geo/IBDverse/raw). Amendment A14.
A17 targets `ibd_ileum` (healthy), `ibd_crohn`, `ibd_crohn_inflamed`, `ibd_crohn_uninflamed` (donor columns of the same
names). Legacy builder replaced: desync/16_ibd/scripts/ibd.py. Choices: see `_atlas/ibd.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import ibd as M  # noqa: E402

DATASET = "ibd_ileum"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
