"""Periodontitis, gingivitis and normal gingiva / oral mucosa atlas (Census 9fcb0b73): donor x class pseudobulks (clock genes followed by the 359 A51 score genes); held-out cohort of A51.
Raw inputs: see `_atlas/heldout.py`. Legacy builders replaced: desync/validation/immune_damping/scripts/ext_heldout.py and the
loaders of desync/validation/ibd_taurus/scripts/a45_taurus.py, desync/validation/masld/scripts/a44_pooled.py.
Choices: see `_atlas/heldout.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import heldout as M  # noqa: E402

DATASET = "periodontitis_scored"


def build(donors=None, genes=None):
    """(pb, donors) from the group pseudobulks on disk; `donors`: slice (see `_atlas/heldout.py`)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the A51 table of the legacy run ($SCCC_DATA/immune_damping/tables)."""
    return M.from_legacy(DATASET)
