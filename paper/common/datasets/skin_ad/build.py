"""Atopic dermatitis skin, lesional vs control biopsies (GSE328048, 10x 5'): donor x class pseudobulks. Raw input
$SCCC_DATA/geo/GSE328048/pseudobulk. Amendment A11. A17 targets `skin_ad_control` (healthy) and `skin_ad`.
Legacy builder replaced: desync/14_skin/scripts/skin.py. Choices: see `_atlas/skin.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import skin as M  # noqa: E402

DATASET = "skin_ad"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
