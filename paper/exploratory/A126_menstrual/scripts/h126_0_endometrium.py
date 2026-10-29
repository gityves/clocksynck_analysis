"""A126 H126.0: within-donor coherence of cycling vs non-cycling epithelial cells, endometrium GSE266265.
Cycling = A24 quartile rule (top quartile of max(S, G2M) among all epithelial cells, > 0; S/G2M = Tirosh lists with expression-
matched controls, `cells_cancer` builder). Non-cycling = every other epithelial cell. python -I h126_0_endometrium.py"""
import sys, importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd, scipy.sparse as sp
import a126_lib as L
E = L.DATA / "geo/GSE266265/export"
OUT = L.OUT
s = importlib.util.spec_from_file_location("cc", L.ROOT / "paper/common/datasets/cells_cancer/build.py")
sys.path.insert(0, str(L.ROOT / "paper/common/datasets"))
CC = importlib.util.module_from_spec(s); s.loader.exec_module(CC)

meta = pd.read_csv(E / "meta.csv", usecols=["cell", "sample", "AnnotationMain", "MenstrualCyclePhase_main", "EndometriosisStatus"], index_col=0)
epi = np.flatnonzero((meta.AnnotationMain == "epithelial").to_numpy())
m = meta.iloc[epi]
genes = open(E / "genes.txt").read().split("\n")[:-1]
p = np.load(E / "counts_p.npy"); I = np.memmap(E / "counts_i.int32", dtype="<i4", mode="r"); X = np.memmap(E / "counts_x.float32", dtype="<f4", mode="r")
cols = [sp.csc_matrix((np.asarray(X[p[c]:p[c + 1]], np.float64), np.asarray(I[p[c]:p[c + 1]]), np.array([0, p[c + 1] - p[c]])),
                      shape=(len(genes), 1)) for c in epi]
M = sp.hstack(cols).tocsr().T.tocsr()  # cells x genes
del cols
depth = np.asarray(M.sum(1)).ravel()
print("epithelial cells", M.shape, "min depth", depth.min(), flush=True)
sc = sp.diags(1e4 / np.maximum(depth, 1)) @ M; sc.data = np.log1p(sc.data)
log_mean = np.asarray(sc.mean(0)).ravel()
sets = CC.control_sets(genes, log_mean, np.random.default_rng(0), {"S": CC.S_GENES, "G2M": CC.G2M_GENES})
sg = CC.scores(M, depth, sets)
print("S/G2M genes found", len(sets["S"][0]), len(sets["G2M"][0]), flush=True)
mx = np.maximum(sg.S, sg.G2M).to_numpy()
q75, q50 = np.quantile(mx, .75), np.median(mx)
cyc = (mx >= q75) & (mx > 0)
g1 = (mx <= q50) & (sg.S.to_numpy() <= 0) & (sg.G2M.to_numpy() <= 0)  # Seurat G1 (both scores <= 0) and <= median
cc_ = L.canonical
gi = {cc_(g): i for i, g in reversed(list(enumerate(genes)))}
clock = [gi[cc_(g)] for g in L.GENES]
C = M[:, clock].toarray()
cells = pd.DataFrame(C, columns=L.GENES).assign(donor=m["sample"].to_numpy(), depth=depth, cyc=cyc, g1=g1, S=sg.S.to_numpy(), G2M=sg.G2M.to_numpy(),
                                                 phase=m.MenstrualCyclePhase_main.to_numpy(), status=m.EndometriosisStatus.to_numpy())
cells.to_parquet(OUT / "endometrium_epithelial_cells_clock.parquet")  # per-cell clock counts (no breast, public GEO data)
print("q75", q75, "cycling cells", cyc.sum(), flush=True)
# pseudo-clock pool statistics (log CPM over epithelial cells) for the companion, not registered for H126.0 (not used)
cell_n = cells.groupby("donor").agg(n_epi=("cyc", "size"), n_cyc=("cyc", "sum"), phase=("phase", "first"), status=("status", "first"))
cell_n["n_g1"] = cells.groupby("donor").g1.sum()
cell_n["n_non"] = cell_n.n_epi - cell_n.n_cyc
cell_n.to_csv(OUT / "h126_0_eligibility_counts.csv")
print(cell_n.describe().round(1)); print("donors with >=30 cycling:", (cell_n.n_cyc >= 30).sum(), "by phase:", cell_n[cell_n.n_cyc >= 30].phase.value_counts().to_dict())
