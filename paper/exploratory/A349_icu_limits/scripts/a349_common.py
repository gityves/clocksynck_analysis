"""A349 shared setup: A109 Tissues (clock states, 14 genes and core loop), GTEx groups (aggregates only downstream)."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
for d in ("A109_coupling_nulls", "A132_brain_death", "A317_braindeath_confirm", "A336_two_readouts/scripts"):
    sys.path.insert(0, str(EXP / d if "/" in d else EXP / d / "scripts"))
sys.path.insert(0, str(EXP / "A336_two_readouts" / "scripts"))
import numpy as np  # noqa: E402

import _helpers as H  # noqa: E402
import compute_brain_death as B  # noqa: E402
import score_a317 as S  # noqa: E402
from a336_common import core as core_wf  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
SEEDS = {"14 genes": 0, "core": 1}


def tissues():
    """Tissues with clock ('14 genes') and core states; groups fast / ventilated post-mortem (donor indices into T.donors)."""
    T = H.Tissues(1, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    CORE = core_wf()
    T.Z["core"], T.V["core"] = {}, {}
    for t in H.BODY:
        dep = T.rows[t].depth.to_numpy(float)
        T.Z["core"][t], T.V["core"][t] = H._state_set(T.rows[t], [CORE], dep)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    h0 = T.group_donors("vent")
    return T, T.group_donors("fast"), h0[coh[h0] == B.PM]
