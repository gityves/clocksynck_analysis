"""PsychAD, RADC (ROSMAP) cohort (dorsolateral prefrontal cortex snRNA-seq): donor x cortical class pseudobulks.
Source: CELLxGENE Census PsychAD dataset (RADC subset), raw input $SCCC_DATA/atlas/psychad_radc.
Legacy builder replaced: desync/common/run_atlas.py::load + desync/03_cortex_conditions/scripts/brain_conditions.py.
Targets (A17): AD_radc, DM_radc (no healthy fit in the legacy run). Choices: see `_atlas/cortex.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import cortex as M  # noqa: E402

DATASET = "psychad_radc"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
