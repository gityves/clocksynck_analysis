"""Jouffe et al. 2022 (PNAS 119:e2109183119): whole-body Bmal1 knockout x leptin deficiency, liver / muscle / eWAT (GSE150102).

Source: `$SCCC_DATA/jouffe2022_bmal1ko/GSE150102_counts_v2.txt.gz` (raw counts, 21,164 Ensembl genes x 68 libraries;
sample titles `Sample_{BmalWT|BmalKO}_{wt|ob}_{animal id}_{liver|MUSCLE|WAT}_ZT12`). Design: 20-week-old male C57BL/6,
all samples at ZT12 (one time point), Bmal1 WT vs Bmal1 KO (whole body, NOT tissue-specific), leptin wt / ob/ob; 4-6 mice per genotype x leptin arm, 1-3 tissues per mouse (45 mice).
Used by: scouting A129 (a whole-body clock-loss arm with no timed truth: single ZT; descriptive only).

Choices: donor = `{bmal}_{leptin}_m{animal id}`, compartment = liver / muscle / WAT, stratum =
`{bmal}_{leptin}`; depth = library total over all genes. Donors table: `zt` (12), `bmal`, `leptin`.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import mouse_ids  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "jouffe2022_bmal1ko" / "GSE150102_counts_v2.txt.gz"


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    x = pd.read_csv(SRC, sep="\t", index_col=0)
    x = x[[c for c in x.columns if c.startswith("Sample_")]]  # drop external_gene_name, description, gene_biotype
    depth = x.sum()
    clock = mouse_ids.clock_counts(x)
    rec = []
    for c in x.columns:
        m = re.fullmatch(r"Sample_Bmal(WT|KO)_(wt|ob)_(\d+)_?(liver|muscle|wat)_?zt(\d+)", c, flags=re.I)
        rec.append(dict(lib=c, bmal=m[1], leptin=m[2], donor=f"{m[1]}_{m[2]}_m{m[3]}",
                        compartment={"liver": "liver", "muscle": "muscle", "wat": "WAT"}[m[4].lower()], zt=int(m[5])))
    meta = pd.DataFrame(rec).set_index("lib")
    meta["stratum"] = meta.bmal + "_" + meta.leptin
    pb = cs.from_bulk(clock.T, meta, donor="donor", compartment="compartment", stratum="stratum",
                      depth=depth.reindex(meta.index))
    donors = meta.drop_duplicates("donor").set_index("donor")[["zt", "bmal", "leptin"]]
    return pb, donors


if __name__ == "__main__":
    from common import data
    pb, donors = build()
    data.save("jouffe2022_bmal1ko", pb, donors)
    print(pb.groupby(["stratum", "compartment"]).size().unstack(), len(donors))
