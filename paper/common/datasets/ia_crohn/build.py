"""Gut Cell Atlas ileum, Crohn's disease vs normal (Census gut_atlas_a7): donor x class pseudobulks. Raw input
$SCCC_DATA/atlas/gut_atlas_a7. Amendment A15 part 4. A17 target `ia_crohn` (donor column `ia_crohn`).
Legacy builder replaced: desync/17_immune_aging/scripts/immune_aging.py (`crohn`). Choices: see `_atlas/immune_aging.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import immune_aging as M  # noqa: E402

DATASET = "ia_crohn"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
