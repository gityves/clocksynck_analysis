import sys,h5py,numpy as np,pandas as pd,re
sys.path.insert(0,'/home/maxine/Documents/paychere/scCircadianOverview/paper')
from anndata.io import read_elem
from clocksync.genes import gene_columns,canonical
import clocksync as cs
from common import gtex
with h5py.File(gtex.LABELLED,'r') as f: v=read_elem(f['var']).index
def rep(name,g):
    g=list(g); c=gene_columns(v,g); print(name,len(g),'missing',[x for x,k in zip(g,c) if k<0])
rep('clock',cs.gtex().genes)
rep('GC',("FKBP5","TSC22D3","ZBTB16","KLF9","KLF15","DDIT4"))
rep('CONS',["BHLHE41","HLF","PPL","ADIRF-AS1","CFAP70","CYS1","MTHFD1L","MUC20","PLXDC1","POMZP3","RHOBTB3","TXNIP","VWF","WASH8P"])
rep('GEN',["ADIRF-AS1","BHLHE40","BHLHE41","CAMKK1","CHST7","FAM76A","HCP5","HLF","IL1R1","LAMA2","LFNG","LONRF1","MAMDC2","NFIL3","NIPSNAP1","PIK3IP1","PPL","RAB3IL1","SLC43A3"])
src=open('results/figure_2/scripts/gene_rhythms/compute_blood_clock.py').read()
m=re.search(r'SETS = \{(.*?)\n\}',src,re.S); d=eval('{'+m.group(1)+'}')
rep('blood',sorted({g for s in d.values() for g in s.split()}))
ref=pd.read_csv('results/figure_3/output/gtex_cues/a81b_cue_reference.csv',index_col=0)
print(ref.columns.tolist()[:6]); 
for k,x in ref.groupby('set'): rep('cue:'+k,x.index)
