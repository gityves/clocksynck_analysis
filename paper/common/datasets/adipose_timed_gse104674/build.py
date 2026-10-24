"""Timed human subcutaneous adipose, lean healthy vs obese type 2 diabetes: 12 men x 4 biopsies (GSE104674).

Source: `$SCCC_DATA/adipose_timed/GSE104674/` (Stenvers et al. 2019, Diabetologia 62:704; SOLiD RNA-seq, hg19, Lifescope raw gene
counts, 24,364 rows `gene:chr:RefSeq`; GEO series matrix for subject, status, sex, age, BMI, HbA1c). Biopsies on day 2 at ZT 15:30
and on day 3 at ZT 0:15, 5:45, 11:15 (ZT 0 = wake-up; its clock time is not in GEO, so the time is a zeitgeber time and a common
offset to wall clock is unknown). Standard feeding (three meals at ZT 0:30, 6:00, 11:30).
Used by: A113 (WP4), A130 (age and sex available).

Table (raw counts): `unit` = library, `block` = subject (unit of within-subject centring), `zt_h` (ZT of the biopsy mod 24, h),
`hours_since_first` (0, 8.75, 14.25, 19.75), `status`, `value_type` = counts, `depth` = total of all gene counts, the 14 genes as
counts (CIART is C1orf51 in the annotation). Donors table (index = subject): `status`, `sex`, `age`, `bmi`, `hba1c`.
"""

from __future__ import annotations

import gzip

import pandas as pd

import clocksync as cs
from common.paths import DATA

D = DATA / "adipose_timed" / "GSE104674"
ZT = {"D2_ZT_15:30": 15.5, "D3_ZT_0:15": 0.25, "D3_ZT_5:45": 5.75, "D3_ZT_11:15": 11.25}
SINCE = {"D2_ZT_15:30": 0.0, "D3_ZT_0:15": 8.75, "D3_ZT_5:45": 14.25, "D3_ZT_11:15": 19.75}
ALIAS = {"CIART": "C1orf51"}  # annotation of the SOLiD table (CHRONO)


def _annotation() -> pd.DataFrame:
    meta, desc = {}, None
    with gzip.open(D / "GSE104674_series_matrix.txt.gz", "rt") as fh:
        for line in fh:
            v = [s.strip('"') for s in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_characteristics_ch1"):
                meta[v[0].split(": ")[0]] = [s.split(": ", 1)[1] for s in v]
            elif line.startswith("!Sample_description"):
                desc = v
            elif line.startswith("!series_matrix_table_begin"):
                break
    ob = pd.DataFrame(meta)
    ob.index = desc  # library column names of the count table
    return ob


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(libraries, subjects): one row per biopsy with the 14 genes as raw counts, and the donor table."""
    ob = _annotation()
    c = pd.read_csv(D / "GSE104674_rawCounts.txt.gz", sep="\t", index_col=0)
    c.index = c.index.str.split(":").str[0]
    c = c[~c.index.duplicated()]
    genes = list(cs.gtex().genes)
    ob = ob.loc[c.columns]
    t = pd.DataFrame(dict(unit=c.columns, block=ob["subject id #"].to_numpy(), zt_h=ob.timepoint.map(ZT).to_numpy(),
                          hours_since_first=ob.timepoint.map(SINCE).to_numpy(), status=ob["subject status"].to_numpy(),
                          value_type="counts", depth=c.sum().to_numpy()))
    x = c.loc[[ALIAS.get(g, g) for g in genes]].T
    x.columns = genes
    t = pd.concat([t, x.reset_index(drop=True)], axis=1)
    donors = ob.drop_duplicates("subject id #").set_index("subject id #")[
        ["subject status", "gender", "age", "bmi (kg/m2)", "hba1c (mmol/mol)"]]
    donors.columns = ["status", "sex", "age", "bmi", "hba1c"]
    donors[["age", "bmi", "hba1c"]] = donors[["age", "bmi", "hba1c"]].astype(float)
    return t, donors
