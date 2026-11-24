"""Aggregates for slide 3.4b (why louder organs are more coupled): per GTEx peripheral tissue, the shared clock
variance S, the departure D, the phase scatter D implies, its split into amplitude and phase directions, and the
independent hour-scale departure of A259.

- S, D (= 2 psi), coupling: A294b table (`A294_floor/output/tissues.csv`).
- Implied phase scatter [h] = sqrt(D / (2 S)) x 24 / (2 pi): the tangential SD of the departure (D / 2 per axis)
  over the RMS radius of the shared state (sqrt S); small-angle approximation, descriptive (post hoc display).
- Phase fraction = tan / (rad + tan) of the departure (A311, registered): clock, and the median over its 100
  expression-matched sets per tissue.
- A259 (registered): per-target departure SD in hours with its 90 % interval (11 tissues).
GTEx aggregates only (one value per tissue).

    python scripts/organs_departure_compute.py
Output: presentation/output/organs_departure/{tissues,summary}.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "output" / "organs_departure"
H = 24 / (2 * np.pi)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(EXP / "A294_floor/output/tissues.csv", index_col="tissue")[["S", "own", "coupling"]]
    t = t.rename(columns={"own": "D"})
    t["amp"] = np.sqrt(t.S)
    t["phase_h"] = np.sqrt(t.D / (2 * t.S)) * H
    s = pd.read_csv(EXP / "A311_hopf/output/sets.csv")
    s["ftan"] = s.tan / (s.rad + s.tan)
    t["ftan_clock"] = s[s.kind == "clock"].set_index("tissue").ftan
    t["ftan_matched"] = s[s.kind == "matched"].groupby("tissue").ftan.median()
    a = pd.read_csv(EXP / "A259_tissue_deviation/output/score/targets.csv").query("mode == 'target'")
    t = t.join(a.set_index("target")[["hours", "hours_lo", "hours_hi"]].add_prefix("a259_"))
    t.to_csv(OUT / "tissues.csv")
    tests = pd.read_csv(EXP / "A294_floor/output/tests.csv").set_index("test")
    body = t.drop(index=["Whole Blood", "Testis"])
    summ = dict(tissues=len(t), median_D=t.D.median(), geomean_S=float(np.exp(np.log(t.S).mean())),
                slope=tests.loc["descriptive", "slope_log_own_clock_on_log_S"], S_fold=t.S.max() / t.S.min(),
                D_fold=t.D.max() / t.D.min(), phase_h_min=body.phase_h.min(), phase_h_max=body.phase_h.max(),
                ftan_clock_median=t.ftan_clock.median(), ftan_matched_median=t.ftan_matched.median(),
                a259_targets=int(t.a259_hours.notna().sum()),
                rho_phase_a259=float(t.dropna(subset=["a259_hours"])[["phase_h", "a259_hours"]].corr("spearman").iloc[0, 1]))
    pd.DataFrame([summ]).to_csv(OUT / "summary.csv", index=False)
    print(pd.Series(summ).round(3).to_string())
    print(t.dropna(subset=["a259_hours"])[["phase_h", "a259_hours", "a259_hours_lo", "a259_hours_hi"]].round(2).to_string())


if __name__ == "__main__":
    main()
