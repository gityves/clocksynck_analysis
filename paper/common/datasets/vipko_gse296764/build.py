"""Central (VIP-Cre) Bmal1 knockout: SCN, liver and heart at CT0 and CT12, ad libitum or time-restricted feeding (GSE296764).

Source: `$SCCC_DATA/vipko_gse296764/GSE296764_SCN_Liver_Heart_count.xlsx` (raw counts, Ensembl gene id + symbol in the first
columns; 76 libraries) and the series matrix (genotype, feeding, CT of the liver and heart libraries, in the column order of the
sheet; the SCN block is out of order and is not used, see `build`). Library codes `S<n>gan` (liver), `S<n>xin` (heart), `S<n>SCN`: the same `S<n>` is the same
mouse across tissues (26 mice with liver and heart, 24 of them with SCN). Design: WT vs Bmal1-VIPKO (VIP-Cre; a central
clock perturbation, not a peripheral one), AL or TRF (night feeding), CT0 or CT12 only; 3-4 mice per cell. The study
publication is not identified in GEO; GEO gives no age and sex here (check before use).
Used by: scouting A129 (central clock loss; two time points only, so no timed truth beyond a CT0 vs CT12 contrast).

Choices: donor = `S<n>`, compartment = SCN / liver / heart, stratum = `{genotype}_{feeding}`, depth = library total
over all genes. Donors table: `ct` (0 / 12), `genotype`, `feeding`.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import mouse_ids  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "vipko_gse296764"
TISSUE = {"gan": "liver", "xin": "heart", "SCN": "SCN"}


def series_meta() -> pd.DataFrame:
    rows = {}
    with gzip.open(SRC / "GSE296764_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith(("!Sample_title", "!Sample_characteristics_ch1")):
                v = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = v[0].split(":")[0] if line.startswith("!Sample_char") else "title"
                rows[key] = v
    m = pd.DataFrame({k: ([s.split(": ", 1)[1] for s in v] if k != "title" else v) for k, v in rows.items()})
    return m


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    x = pd.read_excel(SRC / "GSE296764_SCN_Liver_Heart_count.xlsx", engine="openpyxl")
    x = x.loc[:, ~x.columns.duplicated()].set_index("name")
    cols = [c for c in x.columns if re.fullmatch(r"S\d+(gan|xin|SCN)", c)]
    cnt = x[cols]
    m = series_meta()
    assert len(m) == len(cols), (len(m), len(cols))
    # series-matrix order = sheet order within a tissue; check tissue blocks agree
    tis = [TISSUE[re.fullmatch(r"S\d+(gan|xin|SCN)", c)[1]] for c in cols]
    assert [t.lower() for t in tis] == [t.lower() for t in m.tissue], "tissue order differs from the series matrix"
    meta = pd.DataFrame(dict(lib=cols, donor=[re.match(r"S\d+", c)[0] for c in cols], compartment=tis,
                             genotype=m.genotype.map({"WT": "WT", "Bmal1-VIPKO": "VIPKO"}).to_numpy(),
                             feeding=m.treatment.map({"ad libitum": "AL", "time-restricted feeding": "TRF"}).to_numpy(),
                             ct=m.time.str[2:].astype(int).to_numpy())).set_index("lib")
    # The SCN block of the series matrix is NOT in the column order of the sheet (the same S<n> gets another genotype than
    # in liver and heart for 7 mice); liver and heart agree for every mouse. Mouse-level attributes therefore come from
    # the liver library of the mouse (heart where there is no liver), the SCN rows of the series matrix are not used.
    ref = meta[meta.compartment != "SCN"].groupby("donor")[["genotype", "feeding", "ct"]]
    chk = ref.nunique().max().max()
    assert chk == 1, "liver and heart disagree"
    ref = ref.first()
    for k in ("genotype", "feeding", "ct"):
        meta[k] = meta.donor.map(ref[k])
    meta = meta.dropna(subset=["genotype"])
    meta["stratum"] = meta.genotype + "_" + meta.feeding
    clock = mouse_ids.clock_counts(cnt)
    pb = cs.from_bulk(clock.T, meta, donor="donor", compartment="compartment", stratum="stratum",
                      depth=cnt.sum().reindex(meta.index))
    donors = meta.drop_duplicates("donor").set_index("donor")[["ct", "genotype", "feeding"]]
    return pb, donors


if __name__ == "__main__":
    from common import data
    pb, donors = build()
    data.save("vipko_gse296764", pb, donors)
    print(pb.groupby(["stratum", "compartment"]).size().unstack(), len(donors))
