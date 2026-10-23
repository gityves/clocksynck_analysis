"""healthy_skin_bins: healthy donors (GSE162183, GSE173706 healthy, Census Reynolds, palm-hip, sole-Boldo) and the uninvolved psoriasis units, dermal classes, kc_cycling and bins kc_d1..kc_d4.

Source: $SCCC_DATA/healthy_skin/pseudobulk_bins.csv (processed output of the legacy raw builders, see `_skin_sc/skin_sc.py`: the raw-data annotation
step is not re-implemented). Used by A63, A64, A69. Unit key `donor` = (study,patient,visit,status).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_skin_sc"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import skin_sc as S  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "healthy_skin/pseudobulk_bins.csv"


def build():
    """(pb, donors): see `_skin_sc/skin_sc.py` for the columns."""
    return S.read(SRC, key=("study","patient","visit","status"))
