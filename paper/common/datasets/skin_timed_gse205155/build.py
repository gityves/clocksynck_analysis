"""Living human skin, epidermis and dermis sampled every 4 h (GSE205155, del Olmo et al. 2022): timed paired layers.

Source: `$SCCC_DATA/skin_paired/GSE205155/` (series matrix of an Agilent array, 11 subjects x 7 times 08-32 h x 2
layers; platform table GPL14550 for the symbols).
Legacy builder replaced: `skin_paired/scripts/skin_paired.py:load` and `table` (A37).
Used by: A37 (untimed dermis-vs-epidermis offset against the timed truth), A67, A76, A78 (timed cue sets), A82.

Choices (as legacy): log2 intensities are turned into pseudo-counts, n_reads x 2^x / sum(2^x) per array with
n_reads = 30e6 (the A37 sensitivities use 1e6 or no subject stratum), the 14 clock probes are unique in the platform
table. donor (`unit`) = subject_hhh (both layers are one sample of one subject at one time), compartment = layer
(epidermis, dermis), stratum = subject; donors table: subject, `label_h` (hours since the first day, 8-32), `clock_h`
(label modulo 24, the time of day). Only units with both layers exist in the series.
"""

from __future__ import annotations

import gzip
import io

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.genes import canonical
from common.paths import DATA

GEO = DATA / "skin_paired" / "GSE205155"
LAYERS = ("epidermis", "dermis")


def _read():
    meta = {}
    with gzip.open(GEO / "GSE205155_series_matrix.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith(("!Sample_title", "!Sample_geo_accession")):
                p = line.rstrip("\n").split("\t")
                meta[p[0]] = [v.strip('"') for v in p[1:]]
    s = pd.DataFrame(dict(gsm=meta["!Sample_geo_accession"], title=meta["!Sample_title"]))
    parts = s.title.str.extract(r"^(epidermis|dermis)_(\d+)h_subject(\d+)$")
    assert parts.notna().all().all(), s[parts.isna().any(axis=1)]
    s["layer"], s["label_h"], s["subject"] = parts[0], parts[1].astype(int), "P" + parts[2]
    s["clock_h"] = s.label_h % 24
    s["unit"] = s.subject + "_" + s.label_h.astype(str).str.zfill(2) + "h"
    m = pd.read_csv(GEO / "GSE205155_series_matrix.txt.gz", sep="\t", comment="!", index_col=0)
    lines = [ln for ln in open(GEO / "GPL14550_table.txt") if ln[0] not in "#^!"]
    g = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False, usecols=["ID", "GENE_SYMBOL"])
    return s, m, g.set_index("ID").GENE_SYMBOL.reindex(m.index)


def build(n_reads: float = 30e6, stratum: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pb, donors): pseudo-count pseudobulk (donor x layer) and the sampling-time table (index: donor)."""
    s, m, sym = _read()
    lin = np.power(2.0, m[s.gsm])
    cnt = lin / lin.sum(0) * n_reads
    names = sym.fillna("").map(canonical).to_numpy()
    cols = {}
    for g in cs.gtex().genes:
        hit = np.flatnonzero(names == g)
        assert len(hit) == 1, (g, len(hit))
        cols[g] = cnt.iloc[hit[0]].to_numpy()
    counts = pd.DataFrame(cols, index=s.gsm.to_numpy())
    samples = s.set_index("gsm")
    pb = cs.from_bulk(counts, samples, donor="unit", compartment="layer", stratum="subject" if stratum else None,
                      depth=pd.Series(float(n_reads), index=samples.index))
    if not stratum:
        pb["stratum"] = ""
    donors = s.drop_duplicates("unit").set_index("unit")[["subject", "label_h", "clock_h"]]
    return pb, donors
