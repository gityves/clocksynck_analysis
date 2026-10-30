"""Sample submission dates and GSM ranges by time-point class (batch-confound check, metadata only). Usage: python -I batch_dates.py <scout_dir>"""
import sys, re, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_meta import read_meta
D = sys.argv[1]
def cls(g, t, ch):
    if g == 'GSE20181': return re.match(r'\d+([A-C])', t).group(1)
    if g == 'GSE105777': return ('control' if 'reanal' in t else 'AI') + (':' + re.search(r'([BSD])( |$)', t).group(1))
    if g == 'GSE111563': return 'tp' + t[-1]
    if g == 'GSE93204': return re.search(r'_(\w+)$', t).group(1)
    if g == 'GSE87411': return re.search(r'_(\w[\w -]*)$', t).group(1)
    if g == 'GSE33658': return 'pre' if 'pre' in t else 'post'
    if g in ('GSE59515',): return t.split('.')[2] if t.count('.') > 2 else t
    if g == 'GSE55374': return [c for c in ch if c.startswith('treatment')][0]
for g in ['GSE20181','GSE59515','GSE55374','GSE105777','GSE111563','GSE93204','GSE87411','GSE33658']:
    m = read_meta(f'{D}/{g}/{g}_series_matrix.txt.gz')
    gsm = m['!Sample_geo_accession'][0]; t = m['!Sample_title'][0]; sub = m['!Sample_submission_date'][0]
    chs = list(zip(*m.get('!Sample_characteristics_ch1', [[]*len(t)]))) if m.get('!Sample_characteristics_ch1') else [()]*len(t)
    d = collections.defaultdict(lambda: [set(), []])
    for i in range(len(t)):
        k = cls(g, t[i], list(chs[i])); d[k][0].add(sub[i]); d[k][1].append(int(gsm[i][3:]))
    print(g, {k: (sorted(v[0]), (min(v[1]), max(v[1]), len(v[1]))) for k, v in sorted(d.items())})
