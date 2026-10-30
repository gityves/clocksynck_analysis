"""Metadata-only: cells per patient x class for Wu 2021 (GSE176078). No expression is read."""
import sys, pandas as pd
md = pd.read_csv(sys.argv[1], index_col=0, usecols=['Unnamed: 0','orig.ident','subtype','celltype_major'])
t = md.groupby(['orig.ident','subtype','celltype_major']).size().unstack(fill_value=0)
t = t.reset_index().set_index('orig.ident')
t.to_csv(sys.argv[2])
cls=['Cancer Epithelial','CAFs','Endothelial','PVL']
print(md.celltype_major.value_counts().to_string())
for s,g in t.groupby('subtype'):
    print(s,len(g),'patients; median cells',{c:int(g[c].median()) for c in cls},'n>=50',{c:int((g[c]>=50).sum()) for c in cls},
          'all three (Cancer,CAF,Endo)>=50:',int(((g['Cancer Epithelial']>=50)&(g['CAFs']>=50)&(g['Endothelial']>=50)).sum()))
