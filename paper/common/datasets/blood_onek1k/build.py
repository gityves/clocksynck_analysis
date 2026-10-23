"""OneK1K PBMC (Yazar et al. 2022, 981 donors, 10x 3' v2; Census 3faad104-...): donor x immune class pseudobulks. Raw input
$SCCC_DATA/atlas/onek1k. Amendment A6 (d). A17 healthy target `blood_onek1k`.
Legacy builder replaced: desync/09_blood_onek1k/scripts/blood.py. Choices: see `_atlas/blood.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import blood as M  # noqa: E402

DATASET = "blood_onek1k"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
