"""AIDA phase 1 PBMC (Asian Immune Diversity Atlas, Census c838aec3-..., 625 healthy donors, 10x 5' v2): donor x immune
class pseudobulks. Raw input $SCCC_DATA/atlas/aida_a7. Amendment A9. A17 healthy target `blood_aida`.
Legacy builder replaced: desync/09_blood_onek1k/scripts/blood_a9.py. Choices: see `_atlas/blood.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import blood as M  # noqa: E402

DATASET = "blood_aida"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
