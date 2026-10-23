"""SEA-AD, dorsolateral prefrontal cortex (snRNA-seq, Gabitto 2024): donor x cortical class pseudobulks.
Source: CELLxGENE Census SEA-AD DLPFC dataset, raw input $SCCC_DATA/atlas/seaad_dlpfc.
Legacy builder replaced: desync/common/run_atlas.py::load + desync/03_cortex_conditions/scripts/brain_conditions.py.
Targets (A17): seaad_dlpfc (healthy), AD_seaad_dlpfc (dementia vs no dementia), ADNC_seaad_dlpfc (high vs not
AD/low/intermediate neuropathology). Choices: see `_atlas/cortex.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import cortex as M  # noqa: E402

DATASET = "seaad_dlpfc"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
