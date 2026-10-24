"""A134 scout: clock-gene loci, gene-set versions, D134.3 cell counts. Writes output/scout_{clock_loci,geneset_versions,d134_3}.csv."""
import gzip
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import pandas as pd  # noqa: E402

from common.paths import DATA  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
REF = DATA / "reference"
GENES = "ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART".split()
NCBI = {"ARNTL": "BMAL1"}
V2026 = Path("/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/validation/glucocorticoid/data/h.all.v2026.1.Hs.symbols.gmt")

# loci (NCBI gene_info map_location; GENCODE v46 position)
gi = pd.read_csv(REF / "Homo_sapiens.gene_info.gz", sep="\t", usecols=["#tax_id", "Symbol", "chromosome", "map_location"])
gi = gi[gi["#tax_id"] == 9606].drop_duplicates("Symbol").set_index("Symbol")
pos = pd.read_csv(REF / "gene_positions_gencode46.tsv", sep="\t", header=None, names=["id", "symbol", "chrom", "start", "end", "type"])
pos = pos.drop_duplicates("symbol").set_index("symbol")
rows = []
for g in GENES:
    s = NCBI.get(g, g)
    band = gi.loc[s, "map_location"]
    rows.append(dict(gene=g, ncbi_symbol=s, chrom=pos.loc[s, "chrom"], start_gencode46=int(pos.loc[s, "start"]), band=band,
                     arm=band[:band.index("p")+1] if "p" in band else band[:band.index("q")+1]))
pd.DataFrame(rows).to_csv(OUT / "scout_clock_loci.csv", index=False)


def gmt(p):
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in Path(p).read_text().splitlines()}


a, b = gmt(REF / "h.all.v2024.1.Hs.symbols.gmt"), gmt(V2026)
rows = []
for k in ("HALLMARK_MYC_TARGETS_V1", "HALLMARK_HYPOXIA"):
    x, y = set(a[k]), set(b[k])
    rows.append(dict(set=k, n_v2024_1=len(x), n_v2026_1=len(y), shared=len(x & y), only_2024_1=" ".join(sorted(x - y)),
                     only_2026_1=" ".join(sorted(y - x)), clock_genes_in_set=" ".join(sorted(x & set(GENES + ["BMAL1"])))))
pd.DataFrame(rows).to_csv(OUT / "scout_geneset_versions.csv", index=False)

# D134.3: cells per patient and group (this scout's pass) and the A18 per-cell groups (>= 100 cells)
import glob  # noqa: E402
c = pd.concat(pd.read_csv(f, dtype={"patient": str}, usecols=["patient", "group", "cells", "study"]) for f in glob.glob(str(OUT / "scout_cells_*.csv")))
w = c.pivot_table(index=["study", "patient"], columns="group", values="cells").reset_index()
g = pd.read_csv(PAPER / "results/figure_2/output/backbone/a18_cancer_groups.csv", dtype={"resection": str})
gm = g[g.compartment == "tumour:malignant"].set_index(["study", "resection"]).cells
gn = g[g.compartment == "adjacent_normal:epithelial"].set_index(["study", "resection"]).cells
w["a18_malignant"] = [gm.get((s, p)) for s, p in zip(w.study, w.patient)]
w["a18_normal"] = [gn.get((s, p)) for s, p in zip(w.study, w.patient)]
w["half_ok_scout"] = (w.malignant >= 100) & (w.normal >= 100)
w["half_ok_a18"] = w.a18_malignant.notna() & w.a18_normal.notna()
w.to_csv(OUT / "scout_d134_3.csv", index=False)
print(w.groupby("study")[["half_ok_scout", "half_ok_a18"]].sum().assign(n=w.groupby("study").size()))
print(w[["half_ok_scout", "half_ok_a18"]].sum(), "of", len(w))
print(pd.read_csv(OUT / "scout_clock_loci.csv").to_string())
print(pd.read_csv(OUT / "scout_geneset_versions.csv").iloc[:, :6].to_string())
