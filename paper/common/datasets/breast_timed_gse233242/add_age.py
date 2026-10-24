"""One-off: add per-patient `age` (y) to resection_times.csv from PNAS SI Appendix Table S1 (Li et al. 2024, 43 pairs).

Source: PDF `pnas.2311854121.sapp-2.pdf` (public SI Appendix; local copy in the archive, `literature/`), text extracted with
`pdftotext -layout`; rows `W00xxxx Female <age> <hh:mm> ...`. Joined on patient id; the resection time of Table S1 is checked
against `resection_time` (must agree for every row). Parsed as text only (python -I), nothing executed.

    python -I add_age.py <sapp_text.txt>
"""
import re
import sys
from pathlib import Path

import pandas as pd

f = Path(__file__).resolve().parent / "resection_times.csv"
rows = [m.groups() for ln in Path(sys.argv[1]).read_text().splitlines()
        if (m := re.match(r"\s*(W00\d{4})\s+(?:Fe)?[Mm]ale\s+(\d+)\s+(\d{1,2}:\d{2})\s", ln))]
s1 = pd.DataFrame(rows, columns=["patient_id", "age", "t"]).drop_duplicates("patient_id").set_index("patient_id")
t = pd.read_csv(f).drop(columns=["age"], errors="ignore")
assert t.patient_id.isin(s1.index).all() and (s1.loc[t.patient_id, "t"].str.zfill(5).to_numpy() == t.resection_time.str.zfill(5).to_numpy()).all()
t.insert(t.columns.get_loc("resection_h") + 1, "age", s1.loc[t.patient_id, "age"].astype(int).to_numpy())
t.to_csv(f, index=False)
print(t.age.describe())
