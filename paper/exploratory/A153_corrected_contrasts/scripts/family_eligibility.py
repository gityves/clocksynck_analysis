"""A153 note (2026-10-30, PI Q15 a): family-level eligibility of disease primaries from the A159 healthy-atlas fits.

C_F = mean anchored coupling of the clock over every compartment of family F in the 5 atlases; C_F,j the same mean under
pseudo set j (set j pooled across atlases). Eligible: C_F > q95 under the matched AND the responsive sets, with >= 3
compartments in >= 2 atlases. Sensitivity (descriptive): the 3 atlases that passed the A159 exact row check.

    python family_eligibility.py
Writes ../output/family_eligibility.csv (scope x family).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[3]
A159 = PAPER / "exploratory" / "A159_specific_coupling" / "output"
OUT = Path(__file__).resolve().parent.parent / "output"
SCOPES = {"all5": ("lung_hlca", "heart", "organ_breast", "organ_muscle", "organ_colon"),
          "rowcheck3": ("heart", "organ_breast", "organ_colon")}

spec = importlib.util.spec_from_file_location("anchor_rule", PAPER / "results/extended_data/scripts/host_identity/anchor_rule.py")
FAM = importlib.util.module_from_spec(spec)
spec.loader.exec_module(FAM)


def main():
    comp = pd.read_csv(A159 / "d159_compartments.csv")
    comp["family"] = comp.compartment.map(FAM.family)
    pseudo = {k: pd.read_csv(A159 / f"pseudo_{k}.csv") for k in ("matched", "responsive")}
    rows = []
    for scope, atlases in SCOPES.items():
        c = comp[comp.atlas.isin(atlases)]
        for fam, g in c.groupby("family"):
            row = dict(scope=scope, family=fam, n_compartments=len(g), n_atlases=g.atlas.nunique(), C=g.c_clock.mean())
            for kind, ps in pseudo.items():
                cols = [(a, k) for a, k in zip(g.atlas, g.compartment)]
                vals = np.vstack([ps.loc[ps.atlas == a].set_index("set")[k].sort_index().to_numpy() for a, k in cols])
                cj = np.nanmean(vals, axis=0)  # set j pooled across atlases; failed fits (NaN) dropped
                row.update({f"{kind}_median": np.median(cj), f"{kind}_q95": np.quantile(cj, 0.95),
                            f"{kind}_p": (1 + np.sum(cj >= row["C"])) / (1 + cj.size), f"{kind}_sets": cj.size})
            row["eligible"] = bool(row["C"] > row["matched_q95"] and row["C"] > row["responsive_q95"]
                                   and row["n_compartments"] >= 3 and row["n_atlases"] >= 2)
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "family_eligibility.csv", index=False)
    print(t.round(3).to_string())


if __name__ == "__main__":
    main()
