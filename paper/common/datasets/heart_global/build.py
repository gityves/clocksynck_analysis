"""Healthy adult heart, left ventricle nuclei (Census heart_global: Litvinukova 2020, Kanemaru 2023, ...): donor x class
pseudobulks. Raw input $SCCC_DATA/atlas/heart_global. Amendment A8. A17 target: `heart` (healthy).
Legacy builder replaced: desync/10_heart/scripts/heart.py (`healthy_groups`). Choices: see `_atlas/heart.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import heart as M  # noqa: E402

DATASET = "heart_global"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
