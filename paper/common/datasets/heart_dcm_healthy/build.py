"""DCM cases (Reichart 2022) against the healthy donors of heart_global: a cross-study contrast, strata prefixed by
study (`reichart2022_dcm|...`, `heart_global|...`). Raw inputs $SCCC_DATA/atlas/{reichart2022_dcm,heart_global}.
Amendment A8 secondary. A17 target: `dcm_healthy`. Legacy builder replaced: desync/10_heart/scripts/heart.py
(`dcm_healthy`). Choices: see `_atlas/heart.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import heart as M  # noqa: E402

DATASET = "heart_dcm_healthy"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
