"""Manella et al. 2021 (Nat Metab 3:829): liver-clock mutants under ad libitum and daytime feeding, 3 tissues, 12 ZT (GSE159135).

Source: `$SCCC_DATA/manella2021/GSE159135_{Liver,WAT,Lung}_{AlbCre,BLKO}_reads.xlsx` (MARS-seq, 3' end; per file the
`*.raw` columns are read counts, `*.normalized` columns are ignored; genes with too few reads were already removed
by the authors: about 11.6 k genes, so `depth` = sum of the retained genes, not the library size).
Design (GEO + paper, Methods "Animals"): 3-4 month old male C57BL/6 Alb-Cre+ (control, "AlbCre") or Alb-Cre+ Bmal1 fl/fl
(BLKO = hepatocyte-specific Bmal1 knockout), ad libitum (AL) or light-phase feeding only (DF) for 30 days; ZT 0-22 every
2 h; 2 mice per ZT x genotype x feeding (labels `A`, `B`); 4 arms x 24 mice = 96 mice, 3 tissues each.
Assumption (not stated by the authors): the label `<ZT><A|B>` is the same mouse in the three tissue files; the sample
sheets are separate per tissue, only the labels match. A decision for the registration.
Used by: A111 (H111.1-H111.3), A129 (descriptive).

Choices: donor = `{genotype}_{feeding}_ZT{zt}_{A|B}`, compartment = Liver / WAT / Lung, stratum = `{genotype}_{feeding}`.
Donors table: `zt`, `genotype` (AlbCre / BLKO), `feeding` (AL / DF), `arm`, `mouse` (A / B).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "manella2021"
TISSUES = ["Liver", "WAT", "Lung"]
GENOTYPES = {"AlbCre": "albcre", "BLKO": "BLKO"}
COL = re.compile(r"(\d+)_(albcre|BLKO)_(AL|DF)_(\d+)([AB])\.raw$")


def read_tissue(tissue: str, genotype: str) -> pd.DataFrame:
    """Raw counts (genes x mice) of one tissue and genotype, columns `{genotype}_{feeding}_ZT{zt}_{A|B}`."""
    x = pd.read_excel(SRC / f"GSE159135_{tissue}_{genotype}_reads.xlsx", index_col=0, engine="openpyxl")
    x = x[[c for c in x.columns if str(c).endswith(".raw")]]
    new = []
    for c in x.columns:
        m = COL.search(c)
        if m is None:
            raise ValueError(f"unexpected column {c!r}")
        new.append(f"{genotype}_{m[3]}_ZT{int(m[4])}_{m[5]}")
    x.columns = new
    x.index = x.index.astype(str)
    return x


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    rec, cnt = [], []
    for t in TISSUES:
        for g in GENOTYPES:
            x = read_tissue(t, g)
            depth = x.sum()
            clock = x[x.index.map(canonical).isin([canonical(s) for s in cs.gtex().genes])]
            for d in x.columns:
                rec.append(dict(donor=d, compartment=t, depth=float(depth[d])))
            cnt.append(clock.T)
    cnt = pd.concat(cnt)
    meta = pd.DataFrame(rec)
    meta["genotype"] = meta.donor.str.split("_").str[0]
    meta["feeding"] = meta.donor.str.split("_").str[1]
    meta["zt"] = meta.donor.str.extract(r"_ZT(\d+)_")[0].astype(int)
    meta["mouse"] = meta.donor.str[-1]
    meta["stratum"] = meta.genotype + "_" + meta.feeding
    cnt.index = range(len(cnt))
    pb = cs.from_bulk(cnt, meta, donor="donor", compartment="compartment", stratum="stratum",
                      depth=meta.depth.set_axis(cnt.index))
    donors = meta.drop_duplicates("donor").set_index("donor")[["zt", "genotype", "feeding", "mouse"]]
    donors["arm"] = donors.genotype + "_" + donors.feeding
    return pb, donors


if __name__ == "__main__":
    from common import data
    pb, donors = build()
    data.save("manella2021", pb, donors)
    print(pb.groupby("compartment").donor.nunique(), len(donors))
