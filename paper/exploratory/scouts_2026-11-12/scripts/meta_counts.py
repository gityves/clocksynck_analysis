"""Count samples/subjects/fields from GEO series-matrix metadata (no expression values read)."""
import sys, gzip, collections
p = sys.argv[1]
rows = collections.defaultdict(list)
with gzip.open(p, 'rt') as f:
    for l in f:
        if l.startswith('!series_matrix_table_begin'): break
        if l.startswith('!Sample_'):
            k, *v = l.rstrip('\n').split('\t'); rows[k].append([x.strip('"') for x in v])
n = len(rows['!Sample_geo_accession'][0]); print('n_samples', n)
ch = rows['!Sample_characteristics_ch1']
import re
for r in ch:
    key = r[0].split(':')[0]
    c = collections.Counter(x.split(': ',1)[-1] for x in r)
    print(key, len(c), dict(c.most_common(8)) if len(c) <= 30 else '(many values)')
