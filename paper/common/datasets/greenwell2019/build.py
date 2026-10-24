"""Greenwell et al. 2019 (eLife 8:e43031): mouse liver under ad libitum, night-restricted and arrhythmic feeding (GSE118967).

Source: `$SCCC_DATA/greenwell2019/GSE118967_RAW_Table.txt.gz` (raw 3' counts, 51,912 Ensembl genes x 72 livers; also
`GSE118967_TPM_Table.txt.gz`). Design (GEO): 12-week-old male C57BL/6J wild type, 4 feeding paradigms x 6 ZT (2, 6, 10, 14,
18, 22) x 3 mice, one liver per mouse (single tissue, so no host across tissues).
Library letter -> arm (GEO titles, in order): A = arrhythmic (AR), D = arrhythmic, "did not adapt fully" (DR),
L = ad libitum (LB), N = night-restricted (NR). Check: the order of the GEO titles/sources is asserted in `arms()`.
Used by: A112 (lock validation). Only a single tissue: the host is the liver's own clock-plane state per sample.

`build()` returns the standard table (compartment `liver`, donor = `{arm}{rep}_ZT{zt}`, depth = library total over all
genes, stratum = arm); `counts()` returns the whole genome (symbols x libraries) for the gene-set analyses of A112.
Donors table: `zt`, `arm`, `rep`.
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

SRC = DATA / "greenwell2019"
ARM = {"A": "arrhythmic", "D": "arrhythmic_unadapted", "L": "ad_libitum", "N": "night_restricted"}


def arms() -> dict:
    """Letter -> arm from the series matrix, asserting the mapping."""
    rows = {}
    with gzip.open(SRC / "GSE118967_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith(("!Sample_title", "!Sample_source_name_ch1")):
                rows[line.split("\t")[0]] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
    tit, src = rows["!Sample_title"], rows["!Sample_source_name_ch1"]
    kind = {"AR": "arrhythmic feeding paradigm;", "DR": "did not adapt", "LB": "ad-libitum", "NR": "night-restricted"}
    for t, s in zip(tit, src):
        k = t[:2]
        assert kind[k] in s, (t, s)
    return {"A": "AR", "D": "DR", "L": "LB", "N": "NR"}


def read_raw() -> pd.DataFrame:
    x = pd.read_csv(SRC / "GSE118967_RAW_Table.txt.gz", sep="\t", index_col=0)
    arms()
    return x


def meta_of(cols) -> pd.DataFrame:
    rec = []
    for c in cols:
        m = re.fullmatch(r"([ADLN])(\d)ZT(\d+)", c)
        rec.append(dict(lib=c, arm=ARM[m[1]], rep=int(m[2]), zt=int(m[3])))
    m = pd.DataFrame(rec).set_index("lib")
    m["donor"] = [f"{a}{r}_ZT{z}" for a, r, z in zip(m.arm, m.rep, m.zt)]
    return m


def counts() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Whole-genome raw counts (mouse symbols x libraries, duplicates summed) and library metadata."""
    x = read_raw()
    sym = mouse_ids.to_symbols(x.index)
    x = x.groupby(sym.fillna("").to_numpy()).sum()
    return x.drop(index="", errors="ignore"), meta_of(x.columns)


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    x = read_raw()
    depth = x.sum()
    clock = mouse_ids.clock_counts(x)
    meta = meta_of(x.columns)
    meta["compartment"] = "liver"
    meta["stratum"] = meta.arm
    pb = cs.from_bulk(clock.T, meta, donor="donor", compartment="compartment", stratum="stratum",
                      depth=depth.reindex(meta.index))
    donors = meta.set_index("donor")[["zt", "arm", "rep"]]
    return pb, donors


if __name__ == "__main__":
    from common import data
    pb, donors = build()
    data.save("greenwell2019", pb, donors)
    print(donors.groupby(["arm", "zt"]).size().unstack(), pb.shape)
