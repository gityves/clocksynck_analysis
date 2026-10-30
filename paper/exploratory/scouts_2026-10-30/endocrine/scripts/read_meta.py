"""Metadata-only reader for GEO series_matrix files. Stops at the expression table: no values are read.
Usage: python -I read_meta.py <series_matrix.txt.gz> [out.tsv]   (prints header summary, writes sample table)"""
import gzip, sys, csv, collections

def read_meta(path):
    meta = collections.OrderedDict()
    with gzip.open(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break  # never read the expression block
            if not line.startswith('!'):
                continue
            k, *v = line.rstrip('\n').split('\t')
            meta.setdefault(k, []).append([x.strip('"') for x in v])
    return meta

def sample_table(meta):
    n = len(meta['!Sample_geo_accession'][0])
    cols = {}
    for k in ('!Sample_geo_accession', '!Sample_title', '!Sample_platform_id', '!Sample_source_name_ch1'):
        if k in meta: cols[k[8:]] = meta[k][0]
    for i, v in enumerate(meta.get('!Sample_characteristics_ch1', [])):
        cols[f'char{i}'] = v
    return n, cols

if __name__ == '__main__':
    meta = read_meta(sys.argv[1])
    n, cols = sample_table(meta)
    print(sys.argv[1], 'n_samples', n, 'license/status:', meta.get('!Series_status'), meta.get('!Series_contact_institute'))
    names = list(cols)
    rows = list(zip(*[cols[c] for c in names]))
    if len(sys.argv) > 2:
        with open(sys.argv[2], 'w', newline='') as o:
            w = csv.writer(o, delimiter='\t'); w.writerow(names); w.writerows(rows)
    for r in rows[:6]: print(r)
