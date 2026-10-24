"""Xin et al. 2021 (iScience 24:102285): day- vs night-restricted feeding, 4 tissues, 24 h (GSE150380, 150381, 151221, 151228).

Source: `$SCCC_DATA/mouse_ageing/Xin2021/` (TPM only, no counts; Ensembl ids, kidney: Entrez ids). Design (GEO + paper):
9-week-old female C57BL/6J, 1 week of daytime (DRF) or nighttime (NRF) restricted feeding, sampled every 4 h at ZT 0-24
(ZT0 and ZT24 are different mice), n = 4 per arm x ZT. Tissues: liver, heart, kidney, visceral adipose tissue (VAT).
Assumption (not stated, not checkable on GEO): replicate `_k` is the same mouse in the four tissues; the four series
are separate. A decision for the registration. (kidney has a replicate `N24_5` and a missing one, see `read_tissue`.)
Used by: A111 replication (pseudo-counts as V25.5).

Choices: pseudo-counts = TPM x `N_READS` / 1e6 rounded (clock genes only; TPM lacks transcript lengths, so no length
correction, unlike V25.5), depth = `N_READS` (30 M, only scales the Poisson covariance); donor =
`{arm}_ZT{zt}_m{k}` (arm DRF/NRF), zt folded: ZT24 stays a separate donor with `zt` = 24 and `zt24` = 0.
Donors table: `zt` (0-24), `zt_fold` (0-20), `arm`, `mouse`.
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

SRC = DATA / "mouse_ageing" / "Xin2021"
FILES = {"liver": ("GSE150380_ncRNAseq-liver-gene-tpm.txt.gz", "ensembl"),
         "heart": ("GSE150381_RNAseq-heart-tpm.txt.gz", "ensembl"),
         "kidney": ("GSE151221_RNAseq-kidney-tpm.txt.gz", "entrez"),
         "VAT": ("GSE151228_RNAseq-VAT-tpm.txt.gz", "ensembl")}
N_READS = 30e6


def read_tissue(tissue: str) -> pd.DataFrame:
    """TPM (ids x libraries) of the 14 clock genes; columns = `{arm}_ZT{zt}_m{k}`."""
    f, how = FILES[tissue]
    x = pd.read_csv(SRC / f, sep="\t", dtype=str, skiprows=[1])
    ids = x["geneID"] if "geneID" in x else x.iloc[:, -1]
    cols = [c for c in x.columns if re.fullmatch(r"[DN]\d+_\d+", c)]
    t = x[cols].astype(float)
    t.index = ids.to_numpy()
    t.columns = [f"{'DRF' if c[0] == 'D' else 'NRF'}_ZT{int(re.search(r'[DN](\d+)_', c)[1])}_m{c.split('_')[1]}"
                 for c in cols]
    return mouse_ids.clock_counts(t, how)


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    parts, rec = [], []
    for t in FILES:
        x = read_tissue(t)
        parts.append((x.T * N_READS / 1e6).round())
        rec += [dict(donor=d, compartment=t) for d in x.columns]
    cnt = pd.concat(parts)
    meta = pd.DataFrame(rec)
    meta["arm"] = meta.donor.str[:3]
    meta["zt"] = meta.donor.str.extract(r"_ZT(\d+)_")[0].astype(int)
    meta["mouse"] = meta.donor.str.extract(r"_m(\d+)$")[0].astype(int)
    meta["stratum"] = meta.arm
    cnt.index = range(len(cnt))
    pb = cs.from_bulk(cnt, meta, donor="donor", compartment="compartment", stratum="stratum",
                      depth=pd.Series(float(N_READS), index=cnt.index))
    donors = meta.drop_duplicates("donor").set_index("donor")[["zt", "arm", "mouse"]]
    donors["zt_fold"] = donors.zt % 24
    return pb, donors


if __name__ == "__main__":
    from common import data
    pb, donors = build()
    data.save("xin2021", pb, donors)
    print(pb.groupby("compartment").donor.nunique(), len(donors))
