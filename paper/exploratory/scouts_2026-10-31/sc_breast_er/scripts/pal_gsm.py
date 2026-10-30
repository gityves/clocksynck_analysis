"""Metadata-only: Pal 2021 (GSE161529) sample table from series matrix."""
import sys, gzip, collections
rows=[l.rstrip('\n').split('\t') for l in gzip.open(sys.argv[1],'rt') if l.startswith('!Sample_characteristics_ch1') or l.startswith('!Sample_title')]
get=lambda key:[r for r in rows if any(c.strip('"').startswith(key) for c in r[1:2])][0][1:]
ct=get('cancer type'); pop=get('cell population'); pt=get('patient'); nc=get('number of cells')
c=collections.defaultdict(set)
for a,b,p in zip(ct,pop,pt): c[(a.strip('"')[13:],b.strip('"')[17:])].add(p.strip('"')[9:])
for k,v in sorted(c.items()): print(k,len(v))
