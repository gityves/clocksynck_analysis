"""Baboon, 12 animals sampled around the clock, 64 tissues (GSE98965, Mure et al. 2018, Science).

Source: `$SCCC_DATA/baboon/GSE98965/` (FPKM matrix, ENA read counts SRP107218, Ensembl 90 PapAnu2.0 GTF and Compara 90
homologies; PROVENANCE.md there). One animal per ZT0-ZT22 (every 2 h) in the series matrix; ZT label = time of day.
Legacy builder replaced: `desync/validation/baboon/scripts/prep.py` (A36; checks against its `clock_counts.csv`,
`tissue_qc.csv`, `gene_map.csv`).
Used by: A36, A40, A41 (baboon waveform, working-hour windows), A42 (cross-species brain vs periphery, couplings).

Choices (as legacy, no time used for any choice): approximate counts = FPKM x union-of-exons length of the gene
(Ensembl 90 PapAnu2.0 GTF), rescaled per library to its sequenced read count (ENA), clock-gene counts rounded;
depth = that read count. The 14 clock genes by Ensembl Compara 90 one-to-one orthologs of the human GRCh38 genes
(CIART = ENSPANG00000000498). donor = animal = `ZTnn`, compartment = tissue abbreviation. A tissue is kept when all 12
animals are present and the median count of every clock gene is at least 5 (63 of 64 tissues; `kept_only=False`
returns all). `approx` marks libraries whose FPKM column carries the series' ".approx" flag. No stratum.
Donors table (index ZTnn): `zt` (h).
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict

import numpy as np
import pandas as pd

import clocksync as cs
from common.paths import DATA

RAW = DATA / "baboon" / "GSE98965"
FPKM = RAW / "GSE98965_baboon_tissue_expression_FPKM.csv.gz"
GTF = RAW / "ensembl90" / "Papio_anubis.PapAnu2.0.90.gtf.gz"
HOMOL = RAW / "ensembl90" / "Compara.90.protein_default.homologies.tsv.gz"
ENA = RAW / "SRP107218_ena_read_counts.tsv"
HUMAN = {"ARNTL": "ENSG00000133794", "CLOCK": "ENSG00000134852", "NPAS2": "ENSG00000170485",
         "PER1": "ENSG00000179094", "PER2": "ENSG00000132326", "PER3": "ENSG00000049246",
         "CRY1": "ENSG00000008405", "CRY2": "ENSG00000121671", "NR1D1": "ENSG00000126368",
         "NR1D2": "ENSG00000174738", "DBP": "ENSG00000105516", "RORA": "ENSG00000069667",
         "TEF": "ENSG00000167074", "CIART": "ENSG00000159208"}
EXPR_MIN = 5.0


def gene_lengths() -> pd.Series:
    """Union-of-exons length per gene."""
    ex = defaultdict(list)
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if p[2] != "exon":
                continue
            g = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
            ex[g].append((int(p[3]), int(p[4])))
    out = {}
    for g, iv in ex.items():
        iv.sort()
        tot, (s, e) = 0, iv[0]
        for a, b in iv[1:]:
            if a > e:
                tot += e - s + 1
                s, e = a, b
            else:
                e = max(e, b)
        out[g] = tot + e - s + 1
    return pd.Series(out)


def gene_map(fpkm: pd.DataFrame) -> pd.DataFrame:
    """Baboon Ensembl IDs of the 14 clock genes (one-to-one orthologs in Compara 90) and their FPKM symbols."""
    want, rows = set(HUMAN.values()), []
    with gzip.open(HOMOL, "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            if p[7] == "homo_sapiens" and p[5] in want and p[4] == "ortholog_one2one":
                rows.append(dict(baboon=p[0], human=p[5]))
    hom = pd.DataFrame(rows)
    hom["symbol"] = hom.human.map({v: k for k, v in HUMAN.items()})
    bab = hom.drop_duplicates("symbol").set_index("symbol").baboon
    missing = [g for g in cs.gtex().genes if g not in bab.index]
    if missing:
        raise ValueError(f"no one-to-one baboon ortholog for {missing}")
    gm = pd.DataFrame(dict(symbol=list(cs.gtex().genes), baboon_id=[bab[g] for g in cs.gtex().genes]))
    gm["fpkm_symbol"] = gm.baboon_id.map(fpkm.set_index("EnsemblID").Symbol)
    return gm


def build(kept_only: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pb, donors): donor x tissue pseudobulk of approximate clock-gene counts, donors indexed by ZT label."""
    x = pd.read_csv(FPKM)
    gm = gene_map(x)
    x = x.set_index("EnsemblID").drop(columns="Symbol")
    ln = gene_lengths().reindex(x.index)
    x, ln = x[ln.notna()], ln[ln.notna()]
    rl = x.mul(ln, axis=0)
    ena = pd.read_csv(ENA, sep="\t").set_index("sample_title").read_count
    meta = []
    for c in x.columns:
        m = re.match(r"([A-Z]+)[._]ZT(\d+)(\.approx)?$", c)
        meta.append(dict(col=c, tissue=m.group(1), zt=int(m.group(2)), approx=bool(m.group(3)),
                         title=f"{m.group(1)}_ZT{m.group(2)}"))
    meta = pd.DataFrame(meta)
    meta["reads"] = meta.title.map(ena)
    if meta.reads.isna().any():
        raise ValueError(f"no read count for {meta[meta.reads.isna()].title.tolist()}")
    meta["animal"] = meta.zt.map(lambda z: f"ZT{z:02d}")  # one animal per ZT (series matrix)
    cnt = rl / rl.sum() * meta.reads.to_numpy()
    clk = cnt.loc[gm.baboon_id].T.round()
    clk.columns = gm.symbol.to_list()
    clk.index = meta.col.to_numpy()
    samples = meta.set_index("col")
    pb = cs.from_bulk(clk, samples, donor="animal", compartment="tissue", depth="reads")
    pb["approx"] = pb.set_index(["donor", "compartment"]).index.map(
        samples.set_index(["animal", "tissue"]).approx).to_numpy()
    if kept_only:
        g = pb.groupby("compartment")
        med = g[list(cs.gtex().genes)].median().min(1)
        keep = med.index[(g.donor.nunique() == 12) & (med >= EXPR_MIN)]
        pb = pb[pb.compartment.isin(keep)].reset_index(drop=True)
    donors = meta.drop_duplicates("animal").set_index("animal")[["zt"]].sort_index()
    return pb, donors
