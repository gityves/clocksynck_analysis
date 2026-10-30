"""Metadata-only: cells per (sample group x seurat cluster) from Pal/Chen Seurat objects (meta.data only; no expression).
python -I pal_cluster_counts.py out.csv obj1.rds obj2.rds ...  (cluster ids are seurat_clusters, 0-based)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rds_parse as R, numpy as np, pandas as pd

def fac(c):
    if 'levels' not in c:
        return np.array([x if isinstance(x, str) or x is None else str(x) for x in c.data], dtype=object)
    lv = list(c['levels'].data); codes = np.asarray(c.data)
    return np.array([lv[i - 1] if i > 0 else None for i in codes], dtype=object)

out = []
for path in sys.argv[2:]:
    o = R.read_rds(path)
    md = o['meta.data']; names = list(md['names'].data)
    cols = {n: c for n, c in zip(names, [md.get(i) for i in range(len(names))])} if False else None
    # list-columns are the Obj's payload (data list)
    vals = dict(zip(names, md.data))
    g, cl = fac(vals['group']), fac(vals['seurat_clusters'])
    obj = os.path.basename(path).replace('SeuratObject_', '').replace('.rds', '')
    df = pd.DataFrame({'group': g, 'cluster': cl}).value_counts().rename('n').reset_index()
    df.insert(0, 'object', obj); out.append(df)
    print(obj, 'cells', len(g), 'columns', names, flush=True)
    del o, md, vals
pd.concat(out).to_csv(sys.argv[1], index=False)
