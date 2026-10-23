"""IBDverse healthy rectum biopsies (scRNA-seq): donor x class pseudobulks. Raw input $SCCC_DATA/geo/IBDverse/pseudobulk.
Amendment A14. A17 target `ibd_rectum` (healthy). Legacy builder replaced: desync/16_ibd/scripts/ibd.py.
Choices: see `_atlas/ibd.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import ibd as M  # noqa: E402

DATASET = "ibd_rectum"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
