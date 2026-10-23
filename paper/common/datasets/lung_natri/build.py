"""Fibrotic lung (ILD) vs control lung (Natri et al. 2024, Census f14bc322-..., 10x 5'): donor x lung class pseudobulks.
Raw input $SCCC_DATA/atlas/natri2024_ild. Amendment A5. A17 target `ild_natri` (donor column `ild_natri`).
Legacy builder replaced: desync/07_lung_ild/scripts/ild.py. Choices: see `_atlas/lung.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import lung as M  # noqa: E402

DATASET = "lung_natri"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
