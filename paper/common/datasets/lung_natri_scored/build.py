"""Natri 2024 fibrotic vs control lung (A20 cohort `ild_natri`) with the A51 score genes: clock genes followed by 359 Hallmark TNFA/IFNG genes (immune-damping test, A51). Same rows as the
A17 table. Raw input $SCCC_DATA/immune_damping/{atlas,geo} (rebuilt group pseudobulks with the extended panel).
Legacy builder replaced: desync/validation/immune_damping/scripts/ext_tables.py. Choices: see `_atlas/scored.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import scored as M  # noqa: E402

DATASET = "lung_natri_scored"


def build(donors=None, genes="scored"):
    """(pb, donors) with the score genes from the rebuilt group pseudobulks; `donors`: slice (donor ids of the base dataset)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy table of the A51 run ($SCCC_DATA/immune_damping/tables)."""
    return M.from_legacy(DATASET)
