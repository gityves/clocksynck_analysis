"""Counts of patients x time points per endocrine-therapy series, from sample metadata only (no expression values read).
Usage: python -I summarise_cohorts.py <scout_dir>   -> prints tables, writes cohort_counts.tsv and gsm_overlap.tsv next to this script's parent."""
import sys, re, collections, itertools, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_meta import read_meta
D = sys.argv[1]
def load(g, fn=None):
    return read_meta(f'{D}/{g}/{fn or g+"_series_matrix.txt.gz"}')
def rows(m):
    gsm = m['!Sample_geo_accession'][0]; title = m['!Sample_title'][0]
    chars = m.get('!Sample_characteristics_ch1', [])
    return gsm, title, chars
out = collections.OrderedDict(); gsms = {}
def report(name, pat_tp, extra=''):
    by = collections.defaultdict(set)
    for p, t in pat_tp: by[p].add(t)
    tps = sorted({t for s in by.values() for t in s})
    print(f'\n== {name}: {len(by)} patients, {len(pat_tp)} samples; {extra}')
    for k in range(1, len(tps)+1):
        for combo in itertools.combinations(tps, k):
            n = sum(1 for s in by.values() if set(combo) <= s)
            if k <= 2 or n: print('  patients with', combo, n) if k<=3 else None
    out[name] = (len(by), len(pat_tp), tps)
    return by

# GSE20181 / GSE5462 (title: "10A;pretreatment..." A=pre, B=10-14 d, C=90 d)
for g in ['GSE20181', 'GSE5462']:
    m = load(g); gsm, title, ch = rows(m); gsms[g] = set(gsm)
    pt = []
    for t in title:
        mm = re.match(r'(\d+)([A-C])[;:]', t)
        pt.append((mm.group(1), {'A': 'pre', 'B': '2wk', 'C': '3m'}[mm.group(2)]))
    report(g, pt, 'responder labels in title: ' + str(collections.Counter(t.split('; ')[-1].strip() for t in title)))
    # note: A/B/C letter mapping is from the title text; verified below for C
    print('  C-title example:', [t for t in title if re.match(r'\d+C;', t)][:1])
# GSE59515
m = load('GSE59515'); gsm, title, ch = rows(m); gsms['GSE59515'] = set(gsm)
pid = [x.split(': ')[1] for x in ch[0]]; tp = [x.split(': ')[1] for x in ch[3]]
by = report('GSE59515', list(zip(pid, tp)))
print('  response:', collections.Counter(zip(pid, [x.split(': ')[1] for x in ch[2]])).most_common(0), collections.Counter(x.split(': ')[1] for x, t in zip(ch[2], ch[3]) if 'pre' in t))
pat59 = set(pid)
# GSE55374
m = load('GSE55374'); gsm, title, ch = rows(m); gsms['GSE55374'] = set(gsm)
pid = [x.split('patient ')[1] for x in ch[0]]; hist = [x.split(': ')[1] for x in ch[1]]; tr = [x.split(': ')[1] for x in ch[2]]
report('GSE55374', list(zip(pid, tr)), str(collections.Counter(hist)))
pat55 = set(pid)
# GSE111563
m = load('GSE111563'); gsm, title, ch = rows(m); gsms['GSE111563'] = set(gsm)
pt = [re.search(r'Patient(\w+)_Timepoint(\d)', t).groups() for t in title]
report('GSE111563', pt)
pat111 = {p for p, _ in pt}
print('  patient-ID overlap (string) 111563 vs 59515:', len(pat111 & pat59), ' vs 55374:', len(pat111 & pat55), ' | 111563 example IDs', sorted(pat111)[:5], ' 59515 example', sorted(pat59)[:5])
# GSE33658
m = load('GSE33658'); gsm, title, ch = rows(m); gsms['GSE33658'] = set(gsm)
ch = {c[0].split(':')[0]: [x.split(': ',1)[1] for x in c] for c in ch}
keys = list(ch); print(keys)
pidk = [k for k in keys if k.startswith('patient id')][0]
by = report('GSE33658', list(zip(ch[pidk], ch['pre_or_post_treatment'])), str(collections.Counter(zip(ch['type_of_rx'], ch['pre_or_post_treatment']))))
pa = {p: ch['type_of_rx'][ch[pidk].index(p)] for p in by}; print('  arms per patient', collections.Counter(pa.values()), ' response', collections.Counter(ch['response'][i] for i in range(len(title)) if ch['pre_or_post_treatment'][i]=='pre'))
# GSE87411
m = load('GSE87411'); gsm, title, ch = rows(m); gsms['GSE87411'] = set(gsm)
pt = []
for t in title:
    mm = re.search(r'Patient ID_(\w+?)_(.*)$', t); pt.append(mm.groups())
print('  tp labels', collections.Counter(t for _, t in pt))
by = report('GSE87411', pt)
# GSE93204
m = load('GSE93204'); gsm, title, ch = rows(m); gsms['GSE93204'] = set(gsm)
pt = [re.search(r'Patient ID_(\w+?)_(\w+)$', t).groups() for t in title]
report('GSE93204', pt)
# GSE105777 + GSE73237
m = load('GSE105777'); gsm, title, ch = rows(m); gsms['GSE105777'] = set(gsm)
print('  105777 title patterns', collections.Counter(re.sub(r'\d+', 'N', t) for t in title))
pt = []; 
for t in title:
    mm = re.match(r'(AI\.Treated|[A-Za-z.]+?)\.?(\d+)([BS])$', t)
    pt.append((mm.group(1) + mm.group(2), mm.group(3)) if mm else ('?' + t, '?'))
by = report('GSE105777', pt)
grp = collections.defaultdict(lambda: collections.Counter())
for p, tpp in pt: grp[re.sub(r'\d+', '', p)][tpp] += 1
print('  per group sample counts', dict(grp))
m = read_meta(f'{D}/GSE73237/GSE73237-GPL10558_series_matrix.txt.gz'); gsm, title, ch = rows(m); gsms['GSE73237_GPL10558'] = set(gsm)
print('  73237 (GPL10558) groups', collections.Counter(c.split(': ')[1] for row in ch for c in row if c.startswith('group')), 'sampling', collections.Counter(c.split(': ')[1] for row in ch for c in row if c.startswith('sampling time')))
t73 = set(title); t105 = set(title_ for title_ in read_meta(f'{D}/GSE105777/GSE105777_series_matrix.txt.gz')['!Sample_title'][0])
print('  title overlap 73237 vs 105777:', len(t73 & t105), sorted(t73)[:6], sorted(t105)[:3])
# GSM overlaps
print('\nGSM overlaps')
for a, b in itertools.combinations(gsms, 2):
    k = len(gsms[a] & gsms[b])
    if k: print(' ', a, b, k)
