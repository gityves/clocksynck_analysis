"""Control lung of Natri et al. 2024 with the T class replaced by non-cycling T cells (natri2024_tcell build): donor x lung
class pseudobulks. Raw input $SCCC_DATA/atlas/natri2024_tcell. Amendment A6. A17 healthy target `lung_natri_control`.
Legacy builder replaced: desync/08_lung_t_residency/scripts/residency.py (`groups`, `load`).
Choices: see `_atlas/lung.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import lung as M  # noqa: E402

DATASET = "lung_natri_control"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
