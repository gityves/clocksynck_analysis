"""Mouse colon chrono-atlas, every gene (CellBender counts): mouse x arm x cell class pseudobulks with known ZT; input of the
cell-type gene-rhythm truth test (A31). Raw input $SCCC_DATA/colon_chrono_atlas/export. Legacy builder replaced:
desync/validation/celltype_truth/scripts/build.py. Choices: see `_atlas/colon_chrono.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import colon_chrono as M  # noqa: E402

DATASET = "colon_chrono_genome"


def build(donors=None, genes=None):
    """(pb, donors) from the streamed export / raw matrices; `donors`: slice (mouse ids, e.g. "F1")."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy pseudobulk table of the known-time validation."""
    return M.from_legacy(DATASET)
