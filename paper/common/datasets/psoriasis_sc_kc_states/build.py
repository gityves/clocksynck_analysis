"""psoriasis_sc_kc_states: keratinocyte states kc_basal, kc_cycling, kc_differentiated (dermal classes come from psoriasis_sc).

Source: $SCCC_DATA/psoriasis_sc/pseudobulk_kc_states.csv (processed output of the legacy raw builders, see `_skin_sc/skin_sc.py`: the raw-data annotation
step is not re-implemented). Used by A61. Unit key `donor` = (study,patient,visit).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_skin_sc"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import skin_sc as S  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "psoriasis_sc/pseudobulk_kc_states.csv"


def build():
    """(pb, donors): see `_skin_sc/skin_sc.py` for the columns."""
    return S.read(SRC, key=("study","patient","visit"))
