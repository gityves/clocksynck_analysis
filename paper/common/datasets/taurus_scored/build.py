"""TAURUS (Thomas et al. 2024; Zenodo 13768607) inflamed vs non-inflamed gut biopsies, unit = biopsy sample: donor x class pseudobulks (clock genes followed by the 359 A51 score genes); A45 and held-out cohort of A51.
Raw inputs: see `_atlas/heldout.py`. Legacy builders replaced: desync/validation/immune_damping/scripts/ext_heldout.py and the
loaders of desync/validation/ibd_taurus/scripts/a45_taurus.py, desync/validation/masld/scripts/a44_pooled.py.
Choices: see `_atlas/heldout.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import heldout as M  # noqa: E402

DATASET = "taurus_scored"


def build(donors=None, genes=None):
    """(pb, donors) from the group pseudobulks on disk; `donors`: slice (see `_atlas/heldout.py`)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the A51 table of the legacy run ($SCCC_DATA/immune_damping/tables)."""
    return M.from_legacy(DATASET)
