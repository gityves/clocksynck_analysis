"""Roadmap item A1 (a): the organ-offset contrasts already scored, in hours, in one table. No new statistic.

Reads the output tables of A208 (sex), A231 (age), A233 (five chronic conditions) and A235 (T2D, obesity; registered host
A164 only) and writes one row per contrast: estimate (h), 90 % interval, detection limit (MDE at 80 % power), permutation
p, equivalence within 1 h as scored, and the pseudo-clock companion (|offset| q95 over 50 matched pseudo-clock sets).
A233 tables carry no interval or MDE for the offset (module offset and p only): left empty.

    python table_hours.py
Output: output/offsets_hours.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

EXP = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1] / "output"


def main():
    rows = []
    t = pd.read_csv(EXP / "A208_sex_coupling/output/tests.csv").set_index("hypothesis").loc["H208.2"]
    rows.append(dict(amendment="A208", contrast="women minus men", cohort="GTEx", readout="mean organ offset",
                     offset_h=t.estimate, lo90=t.lo90, hi90=t.hi90, mde80_h=t.mde_80, p=t.p_perm,
                     equivalent_1h=t.equivalent_within_1h, pseudo_abs_q95_h=t.pseudo_abs_q95,
                     source="A208_sex_coupling/output/tests.csv"))
    t = pd.read_csv(EXP / "A231_age_coupling/output/score_tests.csv").set_index("hypothesis").loc["H231.2"]
    rows.append(dict(amendment="A231", contrast="old minus young", cohort="GTEx", readout="mean organ offset",
                     offset_h=t.estimate, lo90=t.lo90, hi90=t.hi90, mde80_h=t.mde_80, p=t.p_perm,
                     equivalent_1h=t.equivalent_within_1h, pseudo_abs_q95_h=t.pseudo_abs_q95,
                     source="A231_age_coupling/output/score_tests.csv"))
    for _, t in pd.read_csv(EXP / "A233_chronic_conditions/output/score_tests.csv").iterrows():
        rows.append(dict(amendment="A233", contrast=f"{t.condition} minus control", cohort="GTEx",
                         readout="organ module offset", offset_h=t.module_offset_h, p=t.module_offset_p,
                         pseudo_abs_q95_h=t.offset_pseudo_abs_q95,
                         source="A233_chronic_conditions/output/score_tests.csv"))
    d = pd.read_csv(EXP / "A235_metabolic_offset/output/score_tests.csv")
    for _, t in d[d.host == "A164"].iterrows():
        rows.append(dict(amendment="A235", contrast=f"{t.condition} minus control", cohort="GTEx fast deaths",
                         readout="metabolic organ offset", offset_h=t.fast_offset_h, lo90=t.fast_lo90,
                         hi90=t.fast_hi90, mde80_h=t.fast_mde_80, p=t.fast_p_two_sided,
                         equivalent_1h=t.fast_equivalent_1h, pseudo_abs_q95_h=t.offset_pseudo_abs_q95,
                         source="A235_metabolic_offset/output/score_tests.csv"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "offsets_hours.csv", index=False)
    print(out.drop(columns="source").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
