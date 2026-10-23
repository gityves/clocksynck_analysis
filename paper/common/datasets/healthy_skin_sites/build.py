"""healthy_skin_sites: healthy skin by body site (GSE202352 hip, palm, sole; E-MTAB-13085 face), dermal classes and bins.

Source: $SCCC_DATA/healthy_skin/pseudobulk_bins_sites.csv (processed output of the legacy raw builders, see `_skin_sc/skin_sc.py`: the raw-data annotation
step is not re-implemented). Used by A66. Unit key `donor` = (study,patient,site).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_skin_sc"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import skin_sc as S  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "healthy_skin/pseudobulk_bins_sites.csv"


def build():
    """(pb, donors): see `_skin_sc/skin_sc.py` for the columns."""
    return S.read(SRC, key=("study","patient","site"))
