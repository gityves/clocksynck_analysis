"""Living human whole blood sampled every ~3 h under sleep restriction or extension (GSE39445, Moller-Levet 2013).

Source: `$SCCC_DATA/blood_timed/GSE39445/` (Agilent array, quantile-normalised log2, 438 samples, 26 subjects, two
sleep protocols; series matrix and the GPL15331 table for symbols). Clock time, hours awake and the melatonin
circadian phase are given per sample.
Legacy builder replaced: `desync/validation/blood_composition/scripts/external.py:read_matrix` and
`blood_clock/scripts/blood_clock.py:load` (symbols averaged per gene, known aliases applied).
Used by: A32 (marker-score rhythms, external), A35 (blood clock reads clock time).

Table (arrays, no counts): `unit` = sample, `block` = subject|sleep protocol (the unit of within-block centring),
`condition` = sleep protocol, `clock_h` (time of the draw, h), `awake_h`, `cphase` (melatonin-referenced), the 14 clock
genes as log2 expression (probes of a symbol averaged; BMAL1 -> ARNTL, CHRONO/C1orf51 -> CIART), `value_type` =
log2_intensity (to be used as an array: natural log = log2 x ln 2). `gene_matrix()` returns the genes x samples log2
matrix over all symbols (legacy rule: probes averaged per symbol, genes with any missing value dropped) for the A32
marker scores. Donors table (index = block): `subject`, `protocol`.
"""

from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.genes import canonical
from common.paths import DATA

D = DATA / "blood_timed" / "GSE39445"


def read_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(probe x sample log2 matrix, sample annotation indexed by GSM)."""
    meta, rows, hdr, titles = {}, [], None, None
    with gzip.open(D / "GSE39445_series_matrix.txt.gz", "rt") as fh:
        intab = False
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                for j, v in enumerate(vals):
                    if ": " in v:
                        k, x = v.split(": ", 1)
                        meta.setdefault(k, {})[j] = x
            elif line.startswith("!series_matrix_table_begin"):
                intab = True
                hdr = next(fh).rstrip("\n").replace('"', "").split("\t")
            elif line.startswith("!series_matrix_table_end"):
                break
            elif intab:
                rows.append(line.rstrip("\n").split("\t"))
    X = pd.DataFrame([r[1:] for r in rows], index=[r[0].strip('"') for r in rows], columns=hdr[1:])
    X = X.replace("", np.nan).astype(float)
    ob = pd.DataFrame(meta).sort_index()
    ob.index = hdr[1:]
    ob["title"] = titles
    ob["subject"] = [re.search(r"(AF\d+)", t).group(1) for t in titles]
    hh = ob.timesampletaken.str.split(":", expand=True).astype(float)
    ob["clock_h"] = (hh[0] + hh[1] / 60) % 24
    ob["awake_h"] = ob.hoursawake.astype(float)
    ob["cphase"] = ob.circadianphase.astype(float)
    ob["block"] = ob.subject + "|" + ob.sleepprotocol
    return X, ob


def _symbols() -> pd.Series:
    lines = (D / "GPL15331_table.txt").read_text().splitlines()
    i0 = lines.index("!platform_table_begin") + 1
    i1 = lines.index("!platform_table_end") if "!platform_table_end" in lines else len(lines)
    hdr = lines[i0].split("\t")
    pl = pd.DataFrame([ln.split("\t")[:len(hdr)] for ln in lines[i0 + 1:i1]], columns=hdr).replace("", np.nan)
    return pl.set_index("ID").GENE_SYMBOL.dropna()


def gene_matrix(alias: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(genes x samples log2, sample annotation): probes averaged per symbol; alias=True maps BMAL1 etc. to the panel
    symbols (A35), alias=False keeps the platform symbols (A32)."""
    X, ob = read_matrix()
    sym = _symbols()
    X = X.loc[X.index.intersection(sym.index)]
    s = sym.loc[X.index]
    if alias:
        s = s.map(lambda g: canonical(g) if canonical(g) in cs.gtex().genes else g)
    G = X.groupby(s.to_numpy()).mean().dropna(axis=0, how="any")
    return G, ob


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(samples, blocks): one row per array with the 14 clock genes (log2), and the subject x protocol table."""
    G, ob = gene_matrix(alias=True)
    genes = list(cs.gtex().genes)
    t = pd.DataFrame({"unit": ob.index.to_numpy(), "block": ob.block.to_numpy(), "condition": ob.sleepprotocol.to_numpy(),
                      "subject": ob.subject.to_numpy(), "clock_h": ob.clock_h.to_numpy(), "awake_h": ob.awake_h.to_numpy(),
                      "cphase": ob.cphase.to_numpy(), "value_type": "log2_intensity", "depth": np.nan})
    t = pd.concat([t, G.loc[genes].T.reset_index(drop=True)], axis=1)
    blocks = t.drop_duplicates("block").set_index("block")[["subject", "condition"]].rename(
        columns={"condition": "protocol"})
    return t, blocks
