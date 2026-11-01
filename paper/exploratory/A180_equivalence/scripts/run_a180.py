"""A180: TOST equivalence of the A19 age null and of the best-powered A165 disease rows (registered 2026-11-01, secondary).

90 % interval (estimate +/- 1.645 SE) inside the bounds. Age: +/- 0.0375 coupling per decade (A164's +/- 0.15 over 40
years). A165: phase +/- 1.5 h; amplitude 1.5-fold, i.e. log z2 ratio within +/- 2 ln 1.5; only primary rows whose MDE is
at or below the bound.
    python run_a180.py   ->   output/tost.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
from common.paths import exploratory_dir  # noqa: E402

Z90 = 1.6448536269514722
BOUND = {"age": 0.15 / 4, "phase": 1.5, "amplitude": 2 * np.log(1.5)}


def tost(est, se, bound):
    lo, hi = est - Z90 * se, est + Z90 * se
    return dict(ci90_lo=lo, ci90_hi=hi, bound=bound, equivalent=bool(-bound < lo and hi < bound))


def main():
    age = pd.read_csv(PAPER / "results/figure_4/output/resilience/a19_age_pooled.csv").iloc[0]
    rows = [dict(test="A19 H19.1 gamma per decade", target="all compartments", readout="age", estimate=age.estimate,
                 se=age.se, mde=np.nan, **tost(age.estimate, age.se, BOUND["age"]))]
    a165 = pd.read_csv(PAPER / "exploratory/A165_disease_readouts/output/h1651.csv")
    a165 = a165[a165.primary & a165.readout.isin(["phase", "amplitude"])]
    for r in a165.itertuples():
        b = BOUND[r.readout]
        row = dict(test="A165 H165.1", target=r.target, compartment=r.compartment, readout=r.readout,
                   estimate=r.estimate, se=r.se, mde=r.mde)
        rows.append({**row, **tost(r.estimate, r.se, b)} if r.mde <= b else
                    {**row, "bound": b, "equivalent": None, "note": "MDE above the bound: not tested"})
    out = pd.DataFrame(rows)
    out.to_csv(exploratory_dir("A180_equivalence") / "tost.csv", index=False)
    tested = out[out.equivalent.notna()]
    print(tested.round(3).to_string())
    print(out.groupby("readout").equivalent.apply(lambda s: f"{int((s == True).sum())} of {int(s.notna().sum())} tested; "
                                                              f"{int(s.isna().sum())} not testable").to_string())


if __name__ == "__main__":
    main()
