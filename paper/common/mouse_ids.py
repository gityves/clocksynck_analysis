"""Mouse gene identifiers -> symbols (NCBI `Mus_musculus.gene_info.gz`, `$SCCC_DATA/mouse_gene_map/`).

Orthologue convention of A48 (`clocksync.genes.canonical`): the 14 clock genes of the GTEx waveform are matched to
mouse genes by upper-cased symbol (Arntl/Bmal1 -> ARNTL, Ciart/Chrono/Gm129 -> CIART). Datasets that come with
Ensembl or Entrez ids (no symbols) are converted here first.
"""

from __future__ import annotations

import gzip
import re

import pandas as pd

from common.paths import DATA

GENE_INFO = DATA / "mouse_gene_map" / "Mus_musculus.gene_info.gz"


def gene_table() -> pd.DataFrame:
    """Columns entrez, symbol, ensembl (unversioned ENSMUSG) of every NCBI mouse gene."""
    rows = []
    with gzip.open(GENE_INFO, "rt") as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            m = re.search(r"Ensembl:(ENSMUSG\d+)", x[5])
            rows.append((x[1], x[2], m.group(1) if m else None))
    return pd.DataFrame(rows, columns=["entrez", "symbol", "ensembl"])


def to_symbols(ids, how: str = "ensembl") -> pd.Series:
    """Symbols of `ids` (Ensembl ids may carry a version suffix; `how` in {"ensembl", "entrez"}); NaN if unknown."""
    t = gene_table().dropna(subset=[how]).drop_duplicates(how).set_index(how).symbol
    key = pd.Series([str(i).split(".")[0] for i in ids])
    return pd.Series(key.map(t).to_numpy(), index=list(ids))


def clock_counts(counts: pd.DataFrame, how: str = "ensembl") -> pd.DataFrame:
    """Subset a genes x libraries matrix indexed by ids to the 14 clock genes (rows named by mouse symbol)."""
    import clocksync as cs
    from clocksync.genes import canonical
    sym = to_symbols(counts.index, how)
    want = {canonical(g) for g in cs.gtex().genes}
    keep = [i for i, s in zip(counts.index, sym) if isinstance(s, str) and canonical(s) in want]
    out = counts.loc[keep].copy()
    out.index = [sym[i] if not isinstance(sym[i], pd.Series) else sym[i].iloc[0] for i in keep]
    return out[~out.index.duplicated()]
