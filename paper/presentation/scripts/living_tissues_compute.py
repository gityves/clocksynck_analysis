"""Per-person clock phases of MuTHER fat, skin and lymphoblastoid cell line, for the example panel of slide 1.6
(descriptive; illustrates A255 H255.1 and H255.2, changes no statistic).

States exactly as A255 (`exploratory/A255_living_tissues/scripts/compute_a255.py`): clock-plane state of the 14 genes
(GTEx waveform, unweighted projector), centred within tissue. Phase = angle of the centred state in hours; modulus = its
length. The script reproduces the scored r of each pair before writing.

Data: open (ArrayExpress E-TABM-1140), prepared by exploratory/A254_paired_lesions/scripts/prep_lesions_tissues.py.
Output: presentation/output/living_tissues/muther_phases.csv (one row per woman: phase_h and modulus per tissue),
presentation/output/living_tissues/check.csv (r reproduced against score_tests.csv).

    python scripts/living_tissues_compute.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
A255 = PAPER / "exploratory/A255_living_tissues"
sys.path.insert(0, str(A255 / "scripts"))
import compute_a255 as a  # noqa: E402

OUT = PAPER / "presentation" / "output" / "living_tissues"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m, _ = a.load_muther(placebo=False)
    z = m.Z()
    scored = pd.read_csv(A255 / "output/score_tests.csv").set_index("test")
    check = []
    for x, y in (("fat", "skin"), ("fat", "lcl"), ("skin", "lcl")):
        r, _, n = a.r_pair(z[x], z[y])
        check.append(dict(pair=f"{x}-{y}", r=r, persons=n, scored_r=scored.loc[f"muther {x}-{y}", "r"]))
    check = pd.DataFrame(check)
    assert np.allclose(check.r, check.scored_r, atol=1e-9), check
    d = pd.DataFrame({"person": m.ids})
    for k, v in z.items():
        d[f"{k}_phase_h"] = (np.angle(v) * a.HOURS) % 24
        d[f"{k}_modulus"] = np.abs(v)
    d.to_csv(OUT / "muther_phases.csv", index=False)
    check.to_csv(OUT / "check.csv", index=False)
    print(check.to_string(index=False))


if __name__ == "__main__":
    main()
