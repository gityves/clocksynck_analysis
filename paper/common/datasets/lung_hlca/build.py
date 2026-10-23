"""Human Lung Cell Atlas (HLCA full, Sikkema et al. 2023, Census 9f222629-...; healthy donors and IPF): donor x lung class
pseudobulks. Raw input $SCCC_DATA/atlas/hlca_full. Amendments A4. A17 targets `lung_hlca` (healthy) and `lung_ipf_hlca`.
Legacy builder replaced: desync/06_lung/scripts/lung_atlas.py. Choices: see `_atlas/lung.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import lung as M  # noqa: E402

DATASET = "lung_hlca"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
