"""A222 descriptive: PER1 loss conditional on TP53 loss (17p13.1) and PER3 loss vs the 1p36 block (GISTIC, all samples)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "paper"))
from common.paths import DATA  # noqa: E402

x = pd.read_csv(DATA / "tcga" / "xena_gistic" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", sep="\t", index_col=0)
x = x[~x.index.duplicated()]
rows = []
for gene, ref in (("PER1", "TP53"), ("PER3", "CHD5"), ("PER3", "TP73")):
    a, b = x.loc[gene].dropna(), x.loc[ref]
    b = b[a.index]
    for lvl, name in ((-1, "loss (<= -1)"), (-2, "deep (-2)")):
        la, lb = a <= lvl, b <= lvl
        rows.append(dict(gene=gene, ref=ref, level=name, n=len(a), freq_gene=la.mean(), freq_ref=lb.mean(),
                         gene_loss_with_ref_loss=(la & lb).sum() / max(la.sum(), 1), gene_given_ref=(la & lb).sum() / max(lb.sum(), 1),
                         gene_given_ref_intact=(la & ~lb).sum() / max((~lb).sum(), 1)))
d = pd.DataFrame(rows).round(4)
d.to_csv(Path(__file__).resolve().parents[1] / "output" / "descriptive_17p_1p36.csv", index=False)
print(d.to_string(index=False))

# truncating share of non-silent calls: clock genes vs background genes with >= 30 non-silent calls (A222 gene cache)
g, t, _, _, _ = pd.read_pickle(Path(__file__).resolve().parents[1] / "output" / "gene_cache.pkl")
t = t.assign(gene=g.gene.to_numpy())
bg = t[(t.nonsilent >= 30) & ~t.gene.isin(["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"])]
share = bg.truncating / bg.nonsilent
ck = t[t.gene.isin(["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "RORA", "TEF", "CIART"])]
tr = pd.DataFrame(dict(gene=ck.gene, nonsilent=ck.nonsilent, truncating=ck.truncating, share=(ck.truncating / ck.nonsilent).round(3),
                       bg_percentile=[round(100 * (share < s).mean(), 1) for s in ck.truncating / ck.nonsilent]))
tr = pd.concat([tr, pd.DataFrame([dict(gene="background", nonsilent=len(bg), truncating=np.nan, share=round(share.median(), 3), bg_percentile=50.0)])])
tr.to_csv(Path(__file__).resolve().parents[1] / "output" / "descriptive_truncating.csv", index=False)
print(tr.to_string(index=False))
