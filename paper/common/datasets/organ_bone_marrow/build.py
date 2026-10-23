"""Healthy bone marrow (Zeng et al. 2025, Census): donor x class pseudobulks (A9 part 2; A17 healthy target `organ_bone_marrow`).
Raw input $SCCC_DATA/atlas/bone_marrow (Census group pseudobulk). Legacy builder replaced: desync/13_healthy_organs/scripts/
organs.py. Choices: see `_atlas/organs.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import organs as M  # noqa: E402

DATASET = "organ_bone_marrow"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
