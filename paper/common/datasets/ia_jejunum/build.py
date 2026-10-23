"""Immune Aging atlas, normal jejunum (lamina propria and epithelium) (Census): donor x class pseudobulks. Raw input
$SCCC_DATA/atlas/immune_aging_a7. Amendment A15. A17 healthy target `ia_jejunum`.
Legacy builder replaced: desync/17_immune_aging/scripts/immune_aging.py. Choices: see `_atlas/immune_aging.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import immune_aging as M  # noqa: E402

DATASET = "ia_jejunum"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
