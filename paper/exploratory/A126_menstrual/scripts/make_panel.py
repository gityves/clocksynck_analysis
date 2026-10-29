"""A126 gene panel for the Kumar Census extraction and Murrow parsing: 14 clock genes, luteal genes, Hallmark G2M, and a pool of candidate
pseudo-clock genes (abundance within the clock genes' range in Reed HR+ luminal pseudobulks; not clock-like, G2M, luteal; seed 126)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd, anndata as ad
import a126_lib as L
from a126_lib import CLOCK_LIKE
sys.path.insert(0, str(L.ROOT / "paper/results/figure_5/scripts/cancer"))
import _shared as sh
a = ad.read_h5ad(L.DATA / "atlas_wg/breast_hbca/pseudobulk.h5ad")
sel = (a.obs.cell_type == "luminal hormone-sensing cell of mammary gland").to_numpy()
X = a.X[sel]; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
dep = a.obs.depth.to_numpy()[sel]
sym = a.var.feature_name.astype(str).str.upper().to_numpy()
frac = np.log((X + .5) / dep[:, None])
mu = pd.Series(frac.mean(0), index=sym); mu = mu[~mu.index.duplicated()]
clock = L.GENES
alias = {"ARNTL": "BMAL1", "CIART": "CHRONO"}
ck = [g if g in mu.index else alias.get(g, g) for g in clock]
print("clock missing", [g for g in ck if g not in mu.index])
lo, hi = mu[[g for g in ck if g in mu.index]].min() - .5, mu[[g for g in ck if g in mu.index]].max() + .5
excl = set(L.canonical(g) for g in CLOCK_LIKE) | set(sh.hallmark_g2m()) | set(L.LUTEAL) | set(sh.PURITY) | set(clock) | set(ck)
excl = {g.upper() for g in excl} | {g for g in mu.index if L.canonical(g) in {L.canonical(x) for x in CLOCK_LIKE}}
pool = [g for g in mu.index[(mu >= lo) & (mu <= hi)] if g not in excl and g.isalnum() is not None and not g.startswith(("MT-", "ENSG"))]
rng = np.random.default_rng(126)
pool = list(rng.choice(pool, min(4000, len(pool)), replace=False))
panel = sorted(set(clock) | set(ck) | set(L.LUTEAL) | set(sh.hallmark_g2m()) | set(pool))
print(len(pool), "pool;", len(panel), "panel; abundance range", round(lo, 2), round(hi, 2))
(L.OUT / "panel_genes.txt").write_text("\n".join(panel) + "\n")
pd.Series(pool).to_csv(L.OUT / "pseudo_pool.csv", index=False, header=False)
