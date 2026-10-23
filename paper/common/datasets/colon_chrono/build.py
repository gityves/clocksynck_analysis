"""Mouse colon chrono-atlas, 14 clock genes (CellBender counts): mouse x arm x cell class pseudobulks with known ZT. Raw input
$SCCC_DATA/colon_chrono_atlas/export. Amendment A21 (known-time validation of disease contrasts; arms C-DC, R-DC, R-PC).
Legacy builder replaced: desync/21_colon_injury/scripts/colon_injury.py::build. Choices: see `_atlas/colon_chrono.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import colon_chrono as M  # noqa: E402

DATASET = "colon_chrono"


def build(donors=None, genes=None):
    """(pb, donors) from the streamed export / raw matrices; `donors`: slice (mouse ids, e.g. "F1")."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy pseudobulk table of the known-time validation."""
    return M.from_legacy(DATASET)
