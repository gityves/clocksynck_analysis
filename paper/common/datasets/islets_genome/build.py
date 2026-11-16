"""Pancreatic islets (T2D atlas), all genes: donor x compartment pseudobulks of every gene (input of `Dataset.gene_rhythms`; cell-type gene rhythms; twin `islets`).
Raw input $SCCC_DATA/atlas_wg/islets_t2d. Legacy builder replaced: desync/validation/gonogo/scripts/sc_predict.py::load and
gonogo_new/scripts/atlases.py::table. Only build() exists (the legacy genome-wide tables were never stored); the
clock-gene columns are checked against the clock-gene dataset of the same rows. Choices: see `_atlas/genome.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import genome as M  # noqa: E402

DATASET = "islets_genome"


def build(donors=None, genes="all"):
    """(pb, donors) of every gene from the public genome-wide group pseudobulk; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, genes)


def from_legacy():
    """Not available: the legacy genome-wide tables were never stored (see `M.check_clock_twin` for the check)."""
    return M.from_legacy(DATASET)
