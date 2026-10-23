"""PsychAD, Mount Sinai cohort (dorsolateral prefrontal cortex snRNA-seq): donor x cortical class pseudobulks.
Source: CELLxGENE Census dataset 37a17b78-4864-4a42-b67b-31c00962795a (PsychAD, Lee/Roussos; HBCC, MSSM and RADC
cohorts are separate datasets here). Raw input: $SCCC_DATA/atlas/psychad_mssm (group pseudobulk of the Census build).
Legacy builder replaced: desync/common/run_atlas.py::load + desync/03_cortex_conditions/scripts/brain_conditions.py.
Targets (A17): psychad_mssm (healthy), AD_mssm, SCZ_mssm, DLBD_mssm, DLBpure_mssm, PD_mssm, VaD_mssm, age_mssm.
Donor columns: age, sex, disease, healthy, one case/control column per contrast. Choices: see `_atlas/cortex.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import cortex as M  # noqa: E402

DATASET = "psychad_mssm"


def build(donors=None, genes="clock"):
    """(pb, donors) from the public group-level pseudobulk on disk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """(pb, donors) from the legacy A17 table(s) ($SCCC_DATA/a17/tables)."""
    return M.from_legacy(DATASET)
