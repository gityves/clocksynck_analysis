"""Normal terminal ileum, epithelium and lamina propria fractions (Kong et al. 2023, Census): donor x class pseudobulks.
Raw inputs $SCCC_DATA/atlas/{kong_ti_epi,ileum_kong,kong_ti_str}. Amendment A16 (primary). A17 target `kong_ileum`.
Legacy builder replaced: desync/18_gut_kong/scripts/kong.py. Choices: see `_atlas/kong.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import kong as M  # noqa: E402

DATASET = "kong_ileum"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
