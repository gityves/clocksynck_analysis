import sys,glob,h5py,pandas as pd,numpy as np
sys.path.insert(0,'/home/maxine/Documents/paychere/scCircadianOverview/paper')
from anndata.io import read_elem
from clocksync.genes import gene_columns,CLOCK_GENES
ref=pd.read_csv('results/figure_3/output/gtex_cues/a81b_cue_reference.csv',index_col=0)
g=list(dict.fromkeys(list(CLOCK_GENES)+list(ref.index[ref.set=='cue'])))
print('requested',len(g))
for p in sorted(glob.glob('/mnt/max8tb/data/CoPhaser_data/cancer_circadian/atlas_wg/*/pseudobulk.h5ad')):
    with h5py.File(p) as f: v=read_elem(f['var'])
    print(p.split('/')[-2],'ensembl-index',(gene_columns(v.index,g)>=0).sum(),'feature_name',(gene_columns(v.feature_name.astype(str),g)>=0).sum())
