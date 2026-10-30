"""Which clock-gene probes exist in each platform annotation, and which of them have a row in each series matrix.
Reads only the annotation (ID, Gene symbol) and the FIRST FIELD (probe id) of series-matrix rows: no expression value is parsed.
Usage: python -I probe_check.py <scout_dir>"""
import gzip, sys, collections
D = sys.argv[1]
GENES = ["ARNTL","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","DBP","RORA","TEF","CIART","BHLHE40","BHLHE41"]
ALIAS = {"BMAL1":"ARNTL","CHRONO":"CIART","C1ORF51":"CIART","DEC1":"BHLHE40","DEC2":"BHLHE41","BHLHB2":"BHLHE40","BHLHB3":"BHLHE41"}
def annot(gpl):
    m = collections.defaultdict(set); hdr = None
    with gzip.open(f'{D}/platforms/{gpl}.annot.gz', 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('ID\t'): hdr = line.rstrip('\n').split('\t'); continue
            if hdr is None or line[0] in '!^#': continue
            r = line.rstrip('\n').split('\t')
            for s in r[hdr.index('Gene symbol')].split('///'):
                s = s.strip().upper(); s = ALIAS.get(s, s)
                if s in GENES: m[s].add(r[0])
    return m
def rowids(path):
    ids = set(); on = False
    with gzip.open(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'): on = True; next(f); continue
            if line.startswith('!series_matrix_table_end'): break
            if on: ids.add(line.split('\t', 1)[0].strip('"'))
    return ids
series = {'GSE20181':'GPL96','GSE5462':'GPL96','GSE59515':'GPL10558','GSE55374':'GPL10558','GSE111563':'GPL10558','GSE105777':'GPL10558',
          'GSE33658':'GPL570','GSE87411':'GPL6480','GSE93204':'GPL6480'}
A = {g: annot(g) for g in sorted(set(series.values()))}
for g, m in A.items():
    print(f'\n{g}: genes with >=1 probe {len([x for x in GENES if x in m])}/{len(GENES)}; missing: {[x for x in GENES if x not in m]}; probes/gene {dict((k, len(v)) for k, v in sorted(m.items()))}')
for s, g in series.items():
    ids = rowids(f'{D}/{s}/{s}_series_matrix.txt.gz')
    pres = {x: len(A[g][x] & ids) for x in GENES if x in A[g]}
    print(f'{s} ({g}): {len(ids)} probe rows; genes with a row in matrix {sum(v>0 for v in pres.values())}/{len(GENES)}; absent from matrix: {[x for x in GENES if pres.get(x,0)==0]}')
