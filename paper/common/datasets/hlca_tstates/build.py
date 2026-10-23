"""HLCA lung with the T class split by per-cell T states (resident-like / circulating-like, core score, CD4 / CD8, non-cycling), healthy and
fibrotic donors: donor x lung class pseudobulks (A6, A7; figure 7 candidate breaks). Raw input $SCCC_DATA/atlas/hlca_tcell.
Legacy builder replaced: desync/08_lung_t_residency/scripts/residency.py. `build(mode=...)` selects the T compartments; no from_legacy
(no legacy table was stored): checked against `hlca_fibrotic`. Choices: see `_atlas/tstates.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import tstates as M  # noqa: E402

DATASET = "hlca_tstates"


def build(donors=None, genes="clock", mode="state", t_min=20):
    """(pb, donors) from the group pseudobulk; `mode`: T compartments (see `_atlas/tstates.py`); `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes, mode=mode, t_min=t_min)


def from_legacy():
    """Not available (no legacy table was stored for the T-state compartments)."""
    return M.from_legacy(DATASET)
