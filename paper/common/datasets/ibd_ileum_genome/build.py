"""IBDverse ileum (healthy and Crohn's), all genes: genome-wide twin of `ibd_ileum` (same groups, classes, conditions, minimum cells and
healthy rule; cell-type program A3). Raw input $SCCC_DATA/atlas_wg/ibdverse (built by
`exploratory/celltype_map_builds/build_ibd_genome.py` with the grouping of the panel build). Choices: see `_atlas/ibd.py`.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import census as C  # noqa: E402
import ibd as M  # noqa: E402

DATASET = "ibd_ileum"
PATH = C.ATLAS_WG / "ibdverse"


def build(donors=None, genes="all"):
    """(pb, donors) of every gene; `donors`: slice (Census donor ids)."""
    return M.build(DATASET, donors, "all", path=PATH)
