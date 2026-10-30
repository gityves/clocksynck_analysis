"""Per Pal patient: cells in tumour clusters (R code Tum/inferCNV clusters), micro-environment clusters (R code Sub), other.
Cluster ids from Chen 2022 R code (hb10x_annot/*.R), converted to 0-based seurat_clusters.
No CAF / endothelial label exists in these objects (meta.data has only integrated cluster ids)."""
import sys, pandas as pd
c = pd.read_csv(sys.argv[1]); c['cluster'] = c.cluster.astype(int)
TUM = {'ERTotal': {0, 4}, 'HER2': {0, 3}, 'TNBC': {0, 2}}          # Tum <- Cluster %in% c(1,5)/c(1,4)/c(1,3), 1-based
NONME = {'ERTotal': {0, 4, 5, 6}, 'HER2': {0, 3, 7}, 'TNBC': {0, 2}}  # complement of Sub <- !(Cluster %in% c(1,5,6,7)/c(1,4,8)/c(1,3))
lab = {'ERTotal': 'ER', 'HER2': 'HER2', 'TNBC': 'TN'}
rows = []
for (o, g), d in c.groupby(['object', 'group']):
    t = d[d.cluster.isin(TUM[o])].n.sum(); m = d[~d.cluster.isin(NONME[o])].n.sum(); tot = d.n.sum()
    rows.append(dict(object=o, sample=g, total=tot, tumour=t, microenv=m, other_nonME=tot - t - m,
                     n_tum_clusters_ge50=int((d[d.cluster.isin(TUM[o])].n >= 50).sum())))
r = pd.DataFrame(rows); r.to_csv(sys.argv[2], index=False)
pd.set_option('display.width', 200); print(r.to_string())
