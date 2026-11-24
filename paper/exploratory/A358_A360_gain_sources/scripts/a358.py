"""A358 (registered, descriptive): is a human organ's gain lower where, in mouse, the organ weighs the drives evenly?

Mouse drive weights b and the amplitude kept under conflict (A322 Manella, A324 Xin) against the GTEx log gain
log sqrt(S) (A294b table), 5 organs mapped as registered. Spearman rho with the exact one-sided p over all orderings.

    python a358.py
Output: ../output/a358_organs.csv, a358_tests.csv.
"""
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
OUT = HERE.parents[1] / "output"
GTEX = {"liver": ["Liver"], "adipose": ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)"], "lung": ["Lung"],
        "heart": ["Heart - Left Ventricle", "Heart - Atrial Appendage"], "kidney": ["Kidney - Cortex"]}


def mouse():
    m = pd.read_csv(EXP / "A322_conflicting_drives/output/fit.csv").iloc[0]
    x = pd.read_csv(EXP / "A322_conflicting_drives/output/xin/fit.csv").iloc[0]
    rm = pd.read_csv(EXP / "A322_conflicting_drives/output/ratios.csv")
    rx = pd.read_csv(EXP / "A322_conflicting_drives/output/xin/ratios.csv")
    return m, x, rm, rx


def exact_p(x, y, less=True):
    r0 = spearmanr(x, y)[0]
    rs = np.array([spearmanr(x, np.asarray(y)[list(p)])[0] for p in permutations(range(len(y)))])
    return r0, float(np.mean(rs <= r0 + 1e-12) if less else np.mean(rs >= r0 - 1e-12))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m, x, rm, rx = mouse()
    b = {"liver": np.mean([m.b_Liver, x.b_liver]), "adipose": np.mean([m.b_WAT, x.b_VAT]), "lung": m.b_Lung,
         "heart": x.b_heart, "kidney": x.b_kidney}
    t = pd.read_csv(EXP / "A294_floor/output/tissues.csv", index_col="tissue")
    g = {o: float(np.mean(np.log(np.sqrt(t.S[ts])))) for o, ts in GTEX.items()}
    am, ax = rm.set_index("tissue").amp_ratio, rx.set_index("tissue").amp_ratio
    kept = {"liver": np.mean([am.Liver, ax.liver]), "adipose": np.mean([am.WAT, ax.VAT]), "lung": am.Lung,
            "heart": ax.heart, "kidney": ax.kidney}
    org = pd.DataFrame({"b": b, "kept": kept, "log_gain": g})
    org["balance"] = org.b * (1 - org.b)
    org.to_csv(OUT / "a358_organs.csv")
    r1, p1 = exact_p(org.balance.to_numpy(), org.log_gain.to_numpy(), less=True)
    r2, p2 = exact_p(org.kept.to_numpy(), org.log_gain.to_numpy(), less=False)
    tests = [dict(test="balance_vs_gain", rho=r1, p_one_sided=p1, consistent=bool(r1 <= -0.7)),
             dict(test="kept_vs_gain (reported)", rho=r2, p_one_sided=p2)]
    pd.DataFrame(tests).to_csv(OUT / "a358_tests.csv", index=False)
    print(org.round(3).to_string())
    print(pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
