"""Metadata-only: ER-positive vs ER-negative arms for Wu 2021 using GEO clinical_subtype; counts per class."""
import sys, gzip, pandas as pd
rows=[l.rstrip('\n').split('\t') for l in gzip.open(sys.argv[1],'rt') if l.startswith('!Sample_title') or l.startswith('!Sample_characteristics_ch1')]
title=[c.strip('"') for c in rows[0][1:]]
sub=[c.strip('"').split(': ',1)[1] for r in rows[1:] if r[1].strip('"').startswith('clinical_subtype') for c in r[1:]]
print(set(c.strip('"').split(':')[0] for r in rows[1:] for c in r[1:2]))
clin=pd.Series(sub,index=title,name='clinical_subtype')
t=pd.read_csv(sys.argv[2],index_col=0).join(clin)
t['ER']=t.clinical_subtype.map(lambda s:'ER+' if 'ER+' in s else 'ER-')
t.to_csv(sys.argv[3])
print(t.groupby(['clinical_subtype','subtype']).size())
cls=['Cancer Epithelial','CAFs','Endothelial']
for s,g in t.groupby('ER'):
    ok=(g[cls]>=50).all(axis=1)
    print(s,len(g),'pts; median',{c:int(g[c].median()) for c in cls},'n>=50',{c:int((g[c]>=50).sum()) for c in cls},'all3>=50:',int(ok.sum()), 'all3>=30:',int((g[cls]>=30).all(axis=1).sum()))
