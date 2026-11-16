"""PsychAD RADC cortex, all genes: genome-wide twin of `psychad_radc` (rules of `cortex.build`; cell-type program A280).
Raw input $SCCC_DATA/atlas_wg/psychad_radc (`exploratory/celltype_map_builds/run_psychad.sh`). Choices: see `_atlas/genome.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import genome as M  # noqa: E402

DATASET = "psychad_radc_genome"


def build(donors=None, genes="all"):
    """(pb, donors) of every gene from the public genome-wide group pseudobulk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    return M.from_legacy(DATASET)
