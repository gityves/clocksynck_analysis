"""Genome-wide pseudobulks of the sheep and baboon atlases for the A263 companions (the dataset builders keep the 14 clock
genes only). Same animals, compartments, inclusion rules, depths and count construction as `common/datasets/sheep/build.py`
and `common/datasets/baboon/build.py`; every gene is kept, so that pseudo-clock sets can be drawn.

- sheep: SheepGTEx v0 counts (NCBI symbols as released), replicate libraries summed, the A34 inclusion rules applied to the
  animal x compartment table of the builder (rows identical to `data._builder("sheep").build()`).
- baboon: approximate counts (FPKM x union-of-exons length, rescaled to the ENA read count) of every baboon gene with a
  one-to-one human ortholog in Compara 90, named by the human gene symbol (GENCODE v29); the 63 tissues of the builder.
Checks: the clock-gene columns equal those of the builders.

    python -I build_genomewide.py
Writes `$SCCC_DATA/derived/A263/{sheep,baboon}_genomewide.parquet` (rows = builder rows; float32).
"""
import gzip
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

DEST = DATA / "derived" / "A263"
CLOCK = list(cs.gtex().genes)


def sheep():
    B = data._builder("sheep")
    pb, _ = B.build()
    meta = pd.read_excel(B.RAW / "SheepGTEx_v0.MetaData.xlsx", sheet_name="Table S1", header=1)
    meta = meta[meta.Filter == "pass"].set_index("Sample")
    comp = meta["Tissue/cell subtype"].replace(B.MERGE)
    key = (meta.Individual.astype(str) + "|" + comp.astype(str)).where(meta.Individual.notna())
    want = pb.donor + "|" + pb.compartment
    keep = key[key.isin(set(want))]
    parts = []
    for ch in pd.read_csv(B.RAW / "SheepGTEx_v0.Gene.count.txt.gz", sep="\t", index_col=0, chunksize=2000):
        ch = ch[[c for c in ch.columns if c in keep.index]].astype(np.float64)
        parts.append(ch.T.groupby(keep.reindex(ch.columns).to_numpy()).sum().T)
    x = pd.concat(parts).T.reindex(want.to_numpy())
    x = x.loc[:, ~x.columns.duplicated()]
    out = pd.concat([pb[["donor", "compartment", "stratum", "depth"]].reset_index(drop=True),
                     x.reset_index(drop=True).astype(np.float32)], axis=1)
    assert np.allclose(out[CLOCK].to_numpy(float), pb[CLOCK].to_numpy(float)), "sheep clock genes differ from the builder"
    return out


def baboon():
    B = data._builder("baboon")
    pb, _ = B.build()
    sym = {}
    with gzip.open(DATA / "gencode_v29" / "gencode.v29.basic.annotation.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] == "gene":
                sym[re.search(r'gene_id "(ENSG\d+)', p[8]).group(1)] = re.search(r'gene_name "([^"]+)"', p[8]).group(1)
    rows = []
    with gzip.open(B.HOMOL, "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            if p[7] == "homo_sapiens" and p[4] == "ortholog_one2one" and p[0].startswith("ENSPANG"):
                rows.append((p[0], sym.get(p[5])))
    hom = pd.DataFrame(rows, columns=["baboon", "symbol"]).dropna().drop_duplicates("baboon")
    hom = hom[~hom.symbol.duplicated(keep=False)].set_index("baboon").symbol
    x = pd.read_csv(B.FPKM).set_index("EnsemblID").drop(columns="Symbol")
    ln = B.gene_lengths().reindex(x.index)
    x, ln = x[ln.notna()], ln[ln.notna()]
    rl = x.mul(ln, axis=0)
    ena = pd.read_csv(B.ENA, sep="\t").set_index("sample_title").read_count
    title = [re.sub(r"[._]ZT(\d+)(\.approx)?$", lambda m: f"_ZT{m.group(1)}", c) for c in x.columns]
    cnt = rl / rl.sum() * ena.reindex(title).to_numpy()
    cnt = cnt.loc[cnt.index.intersection(hom.index)]
    cnt.index = hom.loc[cnt.index].to_numpy()
    cnt = cnt.T
    cnt["donor"] = [f"ZT{int(re.search(r'ZT(\d+)', c).group(1)):02d}" for c in x.columns]
    cnt["compartment"] = [re.match(r"([A-Z]+)", c).group(1) for c in x.columns]
    cnt = cnt.groupby(["donor", "compartment"]).sum()
    cnt = cnt.reindex(pd.MultiIndex.from_frame(pb[["donor", "compartment"]]))
    out = pd.concat([pb[["donor", "compartment", "depth"]].reset_index(drop=True),
                     cnt.reset_index(drop=True).round().astype(np.float32)], axis=1)
    assert np.allclose(out[CLOCK].to_numpy(float), pb[CLOCK].to_numpy(float)), "baboon clock genes differ from the builder"
    return out


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for name, f in (("baboon", baboon), ("sheep", sheep)):
        t = f()
        t.to_parquet(DEST / f"{name}_genomewide.parquet")
        print(name, t.shape, flush=True)


if __name__ == "__main__":
    main()
