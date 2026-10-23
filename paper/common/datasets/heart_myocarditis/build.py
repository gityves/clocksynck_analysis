"""Myocarditis (incl. MIS-C, post-COVID, post-vaccination) vs normal hearts, left ventricle (Census heart_myocarditis):
donor x class pseudobulks. Raw input $SCCC_DATA/atlas/heart_myocarditis. Amendment A8. A17 target: `myocarditis`.
Legacy builder replaced: desync/10_heart/scripts/heart.py (`myocarditis`). Choices: see `_atlas/heart.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import heart as M  # noqa: E402

DATASET = "heart_myocarditis"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
