"""A161 extra overlap check (all probes, no clock readout): near-duplicate arrays between GSE51981 and GSE6364 (same lab, no shared patient ids).
Per sample pair, Spearman correlation over genes common to both gene-level matrices (rank-based, normalisation-independent); a re-deposited
sample would give a cross-series correlation as high as technical replicates. Reports the max cross-series correlation per GSE6364 sample
against the distribution of nearest-neighbour correlations within GSE6364 and within GSE51981. python -I s02b_overlap_expr.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import a161_lib as L
a = pd.read_parquet(L.EN / "derived/GSE51981_gene_log2.parquet"); b = pd.read_parquet(L.EN / "derived/GSE6364_gene_log2.parquet")
g = a.index.intersection(b.index); ra = a.loc[g].rank(); rb = b.loc[g].rank()
C = np.corrcoef(np.hstack([ra.values, rb.values]).T)
na = a.shape[1]; Cab = C[:na, na:]; Caa = C[:na, :na].copy(); Cbb = C[na:, na:].copy()
np.fill_diagonal(Caa, -1); np.fill_diagonal(Cbb, -1)
r = pd.DataFrame({"max_cross_corr_per_GSE6364_sample": Cab.max(0), "max_within_GSE6364": Cbb.max(0)}, index=b.columns)
r.to_csv(L.OUT / "expression_overlap_6364_vs_51981.csv")
s = pd.Series({"n_genes": len(g), "cross_max_overall": Cab.max(), "cross_median_of_per_sample_max": np.median(Cab.max(0)),
               "within_6364_median_nn": np.median(Cbb.max(0)), "within_6364_max": Cbb.max(), "within_51981_median_nn": np.median(Caa.max(0)),
               "within_51981_max": Caa.max()})
s.to_csv(L.OUT / "expression_overlap_summary.csv", header=["value"]); print(s)
