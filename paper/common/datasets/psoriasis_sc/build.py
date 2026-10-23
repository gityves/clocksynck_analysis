"""psoriasis_sc: class (keratinocyte, fibroblast, endothelial, mural, T, myeloid, melanocyte) of whole-skin psoriasis single cells (GSE278330 snRNA-seq, GSE228421 and GSE173706 whole skin).

Source: $SCCC_DATA/psoriasis_sc/pseudobulk.csv (processed output of the legacy raw builders, see `_skin_sc/skin_sc.py`: the raw-data annotation
step is not re-implemented). Used by A58, A60, A61 (status-split host). Unit key `donor` = (study,patient,visit).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_skin_sc"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import skin_sc as S  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "psoriasis_sc/pseudobulk.csv"


def build():
    """(pb, donors): see `_skin_sc/skin_sc.py` for the columns."""
    return S.read(SRC, key=("study","patient","visit"))
