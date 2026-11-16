"""A243, registered descriptive (as `position_a242.py`): where the fetal signed coherence C of each test (`l1.csv`) lies among the adult
stand-ins of the feasibility (`feasibility_subsamples.csv`: 200 GTEx subsamples per test and setting, thinned to the
fetal design). Per test and setting: median adult C and the share (%) of adult subsamples with C at or below the
fetal value. Aggregates only.

    python position_a243.py

Output (`../output/`): position_in_adult.csv.
"""
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "output"

if __name__ == "__main__":
    fetal = pd.read_csv(OUT / "l1.csv").set_index("test").C
    sub = pd.read_csv(OUT / "feasibility_subsamples.csv")
    rows = [dict(test=t, setting=s, fetal_C=fetal[t], adult_median_C=g.C.median(),
                 pct_adult_at_or_below=100 * (g.C <= fetal[t]).mean())
            for (t, s), g in sub.groupby(["test", "setting"], sort=False)]
    pos = pd.DataFrame(rows)
    pos.to_csv(OUT / "position_in_adult.csv", index=False)
    print(pos.round(3).to_string(index=False))
