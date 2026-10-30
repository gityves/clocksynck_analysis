"""Second metadata pass: GSE87411 patient/time-point classes, patient-ID overlaps between Edinburgh series, GSE105777 control = GSE73237 reanalysis,
GSE73237 GPL6884 part. Metadata only. Usage: python -I design_checks.py <scout_dir>"""
import sys, re, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_meta import read_meta
D = sys.argv[1]
# GSE87411
m = read_meta(f'{D}/GSE87411/GSE87411_series_matrix.txt.gz'); t = m['!Sample_title'][0]
cls = {'BL': 'base', 'BL Redo': 'base', 'BL_Redo': 'base', 're-analyzed': 'reanalyzed(?)', 'w2': 'post', 'W2': 'post', 'w2 Redo': 'post', 'w2_Redo': 'post', 'W2_Redo': 'post', 'M': 'post', 'M_Redo': 'post'}
by = collections.defaultdict(set)
for x in t:
    p, tp = re.search(r'Patient ID_(\w+?)_(.*)$', x).groups(); by[p].add(cls[tp])
print('GSE87411 patients', len(by), 'with base+post', sum({'base', 'post'} <= s for s in by.values()), 'with reanalyzed+post', sum({'reanalyzed(?)', 'post'} <= s for s in by.values()),
      'with base+post or reanalyzed+post', sum(('post' in s) and (s & {'base', 'reanalyzed(?)'}) != set() for s in by.values()))
print('  other sample characteristics rows:', [r[0][:60] for r in m.get('!Sample_characteristics_ch1', [])], ' treatment/arm field present:', any('arm' in k.lower() or 'treat' in k.lower() for k in m))
print('  GSE87411 distinct descriptions', set(m['!Sample_description'][0]))
# Edinburgh patient ids
ids = {}
for g, f in [('GSE59515', 0), ('GSE55374', 0)]:
    mm = read_meta(f'{D}/{g}/{g}_series_matrix.txt.gz'); ids[g] = {re.findall(r'\d+', x)[0] for x in mm['!Sample_characteristics_ch1'][0] if re.findall(r'\d+', x)}
mm = read_meta(f'{D}/GSE20181/GSE20181_series_matrix.txt.gz'); ids['GSE20181'] = {re.match(r'(\d+)', x).group(1) for x in mm['!Sample_title'][0]}
mm = read_meta(f'{D}/GSE111563/GSE111563_series_matrix.txt.gz'); ids['GSE111563'] = {re.search(r'PatientX?(\d+)\w?_', x).group(1) for x in mm['!Sample_title'][0]}
for a in ids:
    for b in ids:
        if a < b: print('numeric patient-id overlap', a, b, len(ids[a] & ids[b]), '(ID spaces may differ: weak evidence)')
print('id ranges', {g: (min(map(int, v)), max(map(int, v))) for g, v in ids.items()})
# GSE105777 vs GSE73237
a = read_meta(f'{D}/GSE105777/GSE105777_series_matrix.txt.gz')['!Sample_title'][0]
b = read_meta(f'{D}/GSE73237/GSE73237-GPL10558_series_matrix.txt.gz')
bt = set(b['!Sample_title'][0]); ra = {x.replace(' [reanalysis]', '') for x in a if 'reanalysis' in x}
print('GSE105777 control titles', len(ra), 'matching GSE73237 GPL10558 titles', len(ra & bt))
print('GSE73237 GPL10558 extra titles (not in 105777):', sorted(bt - ra))
print('GSE73237 series summary/design:', b.get('!Series_summary'), b.get('!Series_overall_design'))
c = read_meta(f'{D}/GSE105777/GSE105777_series_matrix.txt.gz')
trt = [x for x in a if x.startswith('AI.Treated')]
print('GSE105777 treated pairs', len(trt)//2, 'HER2/subtype fields:', collections.Counter(v for row in c['!Sample_characteristics_ch1'] for v in row if v.startswith(('her2', 'subtype'))))
