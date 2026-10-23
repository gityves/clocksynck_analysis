"""Pooled MASLD liver cohorts with the A96 gene panel (14 core clock genes + cue genes + locked genes; panel of the h5ad), unit = patient,
arms D1 (`arm`). Raw inputs $SCCC_DATA/masld/a96/<GSE>/pseudobulk (the A38/A43/A44 builders rerun with the A96 panel). Amendment A96 part 2.
Legacy builder replaced: desync/validation/masld/scripts/a96_counts.py (liver). Choices: see `_atlas/heldout.py` (as masld_pooled).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import heldout as M  # noqa: E402

DATASET = "masld_pooled_a96"


def build(donors=None, genes=None):
    """(pb, donors) from the group pseudobulks with the A96 panel; `donors`: slice (see `_atlas/heldout.py`)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from desync/validation/masld/output/a96_counts.csv (liver rows)."""
    return M.from_legacy(DATASET)
