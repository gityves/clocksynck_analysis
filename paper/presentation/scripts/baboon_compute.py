"""Slide data (results, R1.4): host phase of each baboon (Mure 2018 atlas, one animal per ZT, 12 animals, 63 tissues).

Same Dataset as A36 (`results/figure_1/scripts/species/compute_baboon.py`: every tissue an anchor, GTEx waveform). A36 writes the
host-vs-ZT agreement (`a36_readouts.csv`) but not the per-animal phases: this script writes them (public data).

    python scripts/baboon_compute.py
Output: paper/presentation/output/baboon/host_phase.csv (animal, zt_h, host_h = host phase in the waveform frame)
"""
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(PAPER), str(PAPER / "results/figure_1/scripts/species")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import data  # noqa: E402

OUT = PAPER / "presentation" / "output" / "baboon"


def main():
    cs.config.limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    pb, donors = data.load("baboon")
    ds = cs.Dataset(pb, anchors=sorted(pb.compartment.unique()), donors=donors, name="baboon")
    h = ds.host()
    phase = np.mod(np.arctan2(h.y, h.x) * 24 / (2 * np.pi), 24)
    out = pd.DataFrame(dict(zt_h=donors.zt.astype(float).reindex(h.index), host_h=phase)).rename_axis("animal")
    out.to_csv(OUT / "host_phase.csv")
    print(out.round(2).to_string())


if __name__ == "__main__":
    main()
