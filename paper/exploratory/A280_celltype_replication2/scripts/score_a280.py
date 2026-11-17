"""A280: second replication of the lymphocyte ranking (IBDverse ileum and rectum) and the third and fourth cortex
cohorts (PsychAD MSSM, RADC), with the offset-sign tie-break. The A274 script with these atlases.

Registered in paper/docs/registered_experiments.md (A280). `--stage design`: eligibility (twin rule) and pair counts.
`--stage score`: the A274 descriptive tables (family matrix, per atlas, pairs); its q95 columns decide nothing (the
verdicts are read by `calibrated_a280.py`, Q36 (b)).

    python score_a280.py --stage design
Outputs: paper/exploratory/A280_celltype_replication2/output/{twin,design}.csv
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("a280_a274", PAPER / "exploratory/A274_celltype_map/scripts/score_a274.py")
A274 = importlib.util.module_from_spec(spec)
sys.modules["a280_a274"] = A274
spec.loader.exec_module(A274)

A274.AMENDMENT = "A280_celltype_replication2"
A274.SEED = 280
A274.BODY = {"ibd_rectum": "ibd_rectum_genome", "ibd_ileum": "ibd_ileum_genome"}
A274.CORTEX = {"psychad_mssm": "psychad_mssm_genome"}
A274.ALL = {**A274.BODY, **A274.CORTEX}
_anchors = A274.A150.learned_anchors


def _anchors_or_none(name, comps):
    """The pairwise readout uses no anchors; atlases without learned anchors (the IBD targets, A153) stay eligible."""
    try:
        return _anchors(name, comps)
    except ValueError:
        return []


A274.A150.learned_anchors = _anchors_or_none

if __name__ == "__main__":  # verdicts come from calibrated_a280.py (Q36 (b)); `--stage score` gives descriptive tables only
    A274.main()
