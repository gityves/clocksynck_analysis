import sys,glob,h5py,numpy as np
from anndata.io import read_elem
from clocksync.genes import gene_columns,CLOCK_GENES
for p in sorted(glob.glob(sys.argv[1])):
    try:
        with h5py.File(p,'r') as f: v=read_elem(f['var'])
    except Exception as e: print(p,'ERR',e); continue
    fn='feature_name' in v
    cur=(gene_columns(v.index,CLOCK_GENES)>=0).sum()
    aft=(gene_columns(v.feature_name.astype(str),CLOCK_GENES)>=0).sum() if fn else None
    print(p.split('circadian/')[-1],len(v),list(v.columns)[:4],v.index[:2].tolist(),'cur',cur,'fn',aft)
