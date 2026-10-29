"""A127 data preparation (blind to the outcome): read the GSE96058 gene matrix once, keep clock, G2M, purity and the
pseudo-clock pool as float32, parse the sample characteristics, drop replicates, write the feasibility table.
Run: python -I prep_scanb.py   (cwd-independent; does not execute anything downloaded)"""
import gzip, re, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "paper"))
sys.path.insert(0, str(ROOT / "paper/results/figure_5/scripts/cancer"))
import _shared as sh  # noqa

D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scanb_gse96058")
OUT = Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(exist_ok=True)


def series(f):
    rows = {}
    with gzip.open(D / f, "rt") as h:
        for line in h:
            if line.startswith("!Sample_title"):
                rows["title"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                rows["gsm"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                k, _, _ = v[0].partition(": ")
                rows[k] = [x.partition(": ")[2] for x in v]
    return pd.DataFrame(rows)


meta = pd.concat([series("GSE96058-GPL11154_series_matrix.txt.gz").assign(platform="HiSeq"),
                  series("GSE96058-GPL18573_series_matrix.txt.gz").assign(platform="NextSeq")], ignore_index=True)
meta = meta.replace({"NA": np.nan})
print(meta.shape, meta.columns.tolist())

g2m, GENES = sh.hallmark_g2m(), sh.GENES
want_clock = set(GENES)
need = set(GENES) | set(g2m) | set(sh.PURITY)
cols, chunks, names, mean_var = None, [], [], {}
allmean, allsd, allgenes = [], [], []
kept = {}
with gzip.open(D / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", "rt") as h:
    cols = [c.strip('"') for c in h.readline().rstrip("\n").split(",")][1:]
    for line in h:
        g, _, rest = line.partition(",")
        g = g.strip('"')
        x = np.array(rest.split(","), dtype=np.float32)
        allgenes.append(g); allmean.append(x.mean()); allsd.append(x.std())
        if g in need or x.min() > -3.3219 or True:
            kept[g] = x
print("genes", len(allgenes), "samples", len(cols))
E = pd.DataFrame(kept, index=cols).astype(np.float32)  # samples x genes (all genes; reduced below)
stats = pd.DataFrame(dict(mean=allmean, sd=allsd), index=allgenes)
stats = stats[~stats.index.duplicated()]
E = E.loc[:, ~E.columns.duplicated()]
floor = float(E.min().min())
print("floor", floor, "max", float(E.max().max()))

# first-listed sample per patient: replicates are titled "<id>repl"
meta["is_repl"] = meta.title.str.endswith("repl")
meta = meta.set_index("title")
meta = meta[~meta.index.duplicated()]
order = {c: i for i, c in enumerate(cols)}
meta["col_order"] = meta.index.map(order)
meta["patient"] = meta.index.str.replace("repl", "", regex=False)
# pool: genes with a non-floor value in >= 90 % of samples (variance/mean over all columns), clock/G2M/purity excluded later
detected = (E > floor + 1e-3).mean()
pool = detected.index[detected >= 0.9]
keep = sorted(set(pool) | (need & set(E.columns)))
E[keep].astype(np.float32).to_parquet(OUT / "scanb_genes_kept.parquet")
meta.to_csv(OUT / "scanb_samples_parsed.csv")
pd.Series({"floor_log2_fpkm_plus_0.1": floor, "n_genes_matrix": len(allgenes), "n_samples_columns": len(cols),
           "n_genes_kept": len(keep), "n_pool_candidates": len(pool)}).to_csv(OUT / "prep_scale.csv", header=False)
print("kept", len(keep))
