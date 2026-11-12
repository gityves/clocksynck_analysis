"""Sample sheets, design tables and clock-gene detection summaries for the scouted psych / transplant datasets.
No clock statistic is computed: only gene presence, fraction non-zero and median within-sample percentile among all genes.
Usage: python -I scout_build.py   (reads $SCCC_DATA, writes ../output)"""
import gzip, collections, io, re, os, sys
import numpy as np, pandas as pd

D = '/mnt/max8tb/data/CoPhaser_data/cancer_circadian/'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output')
SCR = sys.argv[1] if len(sys.argv) > 1 else '/tmp'          # scratch dir holding GPL570.txt / GPL10558.txt (GEO platform soft)
CLOCK = ['ARNTL', 'CLOCK', 'NPAS2', 'PER1', 'PER2', 'PER3', 'CRY1', 'CRY2', 'NR1D1', 'NR1D2', 'DBP', 'RORA', 'TEF', 'CIART']
ALIAS = {'CIART': ['C1orf51'], 'ARNTL': ['BMAL1']}

def sm_samples(f):
    rows = collections.defaultdict(list)
    for l in gzip.open(f, 'rt', errors='replace'):
        if l.startswith('!series_matrix_table_begin'): break
        p = l.rstrip('\n').split('\t')
        if p[0].startswith('!Sample_'): rows[p[0][8:]].append([x.strip('"') for x in p[1:]])
    d = {'gsm': rows['geo_accession'][0], 'title': rows['title'][0], 'source': rows['source_name_ch1'][0],
         'platform': rows['platform_id'][0]}
    if 'description' in rows: d['description'] = rows['description'][0]
    for c in rows['characteristics_ch1']:
        k = collections.Counter(x.split(':', 1)[0].strip() for x in c).most_common(1)[0][0]
        d[k] = [x.split(':', 1)[1].strip() if ':' in x else x for x in c]
    return pd.DataFrame(d)

def sm_table(f):
    with gzip.open(f, 'rt', errors='replace') as fh:
        for l in fh:
            if l.startswith('!series_matrix_table_begin'): break
        txt = []
        for l in fh:
            if l.startswith('!series_matrix_table_end'): break
            txt.append(l)
    t = pd.read_csv(io.StringIO(''.join(txt)), sep='\t', index_col=0)
    t.index = t.index.astype(str).str.strip('"'); t.columns = t.columns.str.strip('"')
    return t

def ens2sym():
    m = {}
    for l in gzip.open(D + 'gencode_v29/gencode.v29.basic.annotation.gtf.gz', 'rt'):
        if l.startswith('#'): continue
        p = l.split('\t')
        if p[2] != 'gene': continue
        gid = re.search('gene_id "([^"]+)"', p[8]).group(1).split('.')[0]; m[gid] = re.search('gene_name "([^"]+)"', p[8]).group(1)
    return m

def probe2sym(gpl, col):
    seen = False; rows = []
    for l in open(f'{SCR}/{gpl}.txt', errors='replace'):
        if l.startswith('!platform_table_begin'): seen = True; hdr = None; continue
        if not seen: continue
        if l.startswith('!platform_table_end'): break
        p = l.rstrip('\n').split('\t')
        if hdr is None: hdr = p; ci = hdr.index(col); continue
        if len(p) > ci and p[ci]: rows += [(p[0], q.strip()) for q in p[ci].split(' /// ')]
    return rows

def clock_summary(expr, label, count_like=True):
    """expr: genes (symbols, unique) x samples. Percentile of each gene within each sample among all genes (rank of values, ties averaged)."""
    pct = expr.rank(axis=0, pct=True) * 100
    res = []
    for g in CLOCK:
        k = g if g in expr.index else next((a for a in ALIAS.get(g, []) if a in expr.index), None)
        if k is None: res.append(dict(dataset=label, gene=g, present=False)); continue
        v = expr.loc[k]
        res.append(dict(dataset=label, gene=g, present=True, matrix_name=k, n_genes_total=expr.shape[0], n_samples=expr.shape[1],
                        frac_samples_nonzero=(float((v > 0).mean()) if count_like else np.nan), median_pctile=float(pct.loc[k].median())))
    return pd.DataFrame(res)

def collapse(t, pairs):
    pr = pd.DataFrame(pairs, columns=['probe', 'sym']); pr = pr[pr.probe.isin(t.index)]
    pr['m'] = t.mean(axis=1).reindex(pr.probe).values
    pr = pr.sort_values('m', ascending=False).drop_duplicates('sym')          # highest-mean probe per symbol
    return pd.DataFrame(t.loc[pr.probe].values, index=pr.sym.values, columns=t.columns)

def design(s, label, group, region, donor):
    n = s.groupby([group, donor])[region].nunique().rename('n_regions').reset_index()
    tab = pd.crosstab(n[group], n['n_regions']); tab.columns = [f'donors_with_{c}_regions' for c in tab.columns]
    tab['donors_total'] = n.groupby(group).size()
    allr = s[region].nunique(); tab['donors_with_all_%d_regions' % allr] = (n['n_regions'] == allr).groupby(n[group]).sum()
    tab['donors_with_ge2'] = (n['n_regions'] >= 2).groupby(n[group]).sum(); tab['donors_with_ge3'] = (n['n_regions'] >= 3).groupby(n[group]).sum()
    tab.insert(0, 'dataset', label); return tab.reset_index()

def save(df, name): df.to_csv(f'{OUT}/{name}', index=False)

def num(x): return pd.to_numeric(x, errors='coerce')
designs, clocks = [], []
E2S = ens2sym()
P = D + 'psych_multiregion/'

# ---- GSE80655 (Ramaker 2017), RNA-seq counts, Ensembl ids
s = sm_samples(P + 'GSE80655/GSE80655_series_matrix.txt.gz')
s = pd.DataFrame(dict(sample_id=s.gsm, library_id=s.description, donor_id=s.title.str.split('_').str[0], donor_id_derived='title prefix (X####); age/sex/diagnosis/PMI identical across a donor\'s regions',
    region=s['brain region'], group=s['clinical diagnosis'], age=num(s['age at death']), sex=s.gender, PMI=num(s['post-mortem interval']), pH=num(s['brain ph']),
    RIN=np.nan, time_of_death='not recorded', platform=s.platform, data_type='RNA-seq counts (gene)', batch_fields='library_id (SL number) only'))
save(s, 'psych_GSE80655_samples.csv'); designs.append(design(s, 'GSE80655', 'group', 'region', 'donor_id'))
x = pd.read_csv(P + 'GSE80655/GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz', sep='\t', index_col=0)
x.index = x.index.map(lambda g: E2S.get(g.split('.')[0], g)); x = x[~x.index.duplicated()]
assert set(s.library_id) <= set(x.columns); clocks.append(clock_summary(x, 'GSE80655'))

# ---- GSE102556 (Labonte 2017), RNA-seq FPKM, Ensembl ids
s = sm_samples(P + 'GSE102556/GSE102556-GPL11154_series_matrix.txt.gz')
s = pd.DataFrame(dict(sample_id=s.gsm, library_id=np.nan, donor_id=s.title.str.split(':').str[0], donor_id_derived='title prefix (donor number); age/sex/phenotype/PMI identical across a donor\'s regions',
    region=s.tissue, group=s.phenotype, age=num(s.age), sex=s.gender, PMI=num(s.pmi), pH=num(s.ph), RIN=num(s.rin), time_of_death='not recorded', platform=s.platform,
    data_type='RNA-seq FPKM (log-scale values, supplementary fpkmtab)', batch_fields='cause of death; alcohol; drugs; medication; smoking (see md)',
    cause_of_death=s['Cause of death'], medication=s.medication, medication_type=s['medication type'], smoking=s.smoking))
save(s, 'psych_GSE102556_samples.csv'); designs.append(design(s, 'GSE102556', 'group', 'region', 'donor_id'))
x = pd.read_csv(P + 'GSE102556/GSE102556_HumanMDD_fpkmtab.txt.gz', sep='\t', index_col=0)
x = x.drop(columns=['gene_name', 'biotype', 'length'])
x.index = x.index.map(lambda g: E2S.get(g.split('.')[0], g)); x = x[~x.index.duplicated()]
cm = pd.Series(x.columns).str.split('.', n=1, expand=True); cm.columns = ['donor', 'rcode']
print('GSE102556 matrix columns per region code', cm.rcode.value_counts().to_dict(), 'GEO sheet samples', len(s))
cm.to_csv(f'{OUT}/psych_GSE102556_matrix_columns.csv', index=False)
clocks.append(clock_summary(x, 'GSE102556'))

# ---- GSE202537 (Ketchesin), RNA-seq counts, Ensembl ids, striatal subregions, time of death
s = sm_samples(P + 'GSE202537/GSE202537_series_matrix.txt.gz')
s['colid'] = s.title.str.extract(r'\[(\w+)\]')[0]
s = pd.DataFrame(dict(sample_id=s.gsm, library_id=s.colid, donor_id=s['patients_id'], donor_id_derived='GEO field patients_id',
    region=s.tissue, group=s['disease state'], age=num(s.age), sex=s.gender, PMI=num(s.pmi), pH=num(s.ph), RIN=num(s.rin),
    time_of_death='corrected.tod (hours; public GEO field, as deposited by authors)' , tod_value=num(s['corrected.tod']),
    platform=s.platform, data_type='RNA-seq counts (gene)', batch_fields='sequence id; library size; tissue storage time',
    BMI=num(s.bmi), manner_of_death=s['manner of death'], sequence_id=s['sequence id'], library_size=num(s['library size']), tissue_storage_time=num(s.tissuestoragetime), race=s.race))
save(s, 'psych_GSE202537_samples.csv'); designs.append(design(s, 'GSE202537', 'group', 'region', 'donor_id'))
x = pd.read_csv(P + 'GSE202537/GSE202537_NAc_Cau_Put_counts_corrected.csv.gz', index_col=0)
x.index = x.index.map(lambda g: E2S.get(str(g).split('.')[0], g)); x = x[~x.index.duplicated()]
print('GSE202537 matrix cols matched', len(set(s.library_id) & set(x.columns)), 'of', len(s)); clocks.append(clock_summary(x, 'GSE202537'))

# ---- GSE53987 (Lanz), arrays GPL570
s = sm_samples(P + 'GSE53987/GSE53987_series_matrix.txt.gz')
g, r, i = s.title.str.split('_').str[0], s.title.str.split('_').str[1], s.title.str.split('_').str[2]
s = pd.DataFrame(dict(sample_id=s.gsm, library_id=np.nan, donor_id=g + '_' + i, donor_id_derived='title <group>_<index>; age/sex/PMI/pH identical across regions of the same index (inferred, no id in GEO)',
    region=s.tissue, group=s['disease state'], age=num(s.age), sex=s.gender, PMI=num(s.pmi), pH=num(s.ph), RIN=num(s.rin), time_of_death='not recorded',
    platform=s.platform, data_type='Affymetrix HG-U133 Plus 2 (series matrix, normalised)', batch_fields='race', race=s.race))
save(s, 'psych_GSE53987_samples.csv'); designs.append(design(s, 'GSE53987', 'group', 'region', 'donor_id'))
x = collapse(sm_table(P + 'GSE53987/GSE53987_series_matrix.txt.gz'), probe2sym('GPL570', 'Gene Symbol')); clocks.append(clock_summary(x, 'GSE53987', False))

# ---- GSE43974 kidney (transplant): no donor/graft id in the sample sheet
T = D + 'transplant/GSE43974/GSE43974_series_matrix.txt.gz'
s = sm_samples(T)
s = pd.DataFrame(dict(sample_id=s.gsm, library_id=s.description, donor_id=np.nan, donor_id_derived='NOT available: title suffix is a per-group/timepoint counter, not a graft id',
    tissue='kidney biopsy', group=s['donor type (0=living, 1=bd, 2=nhb)'].map({'0': 'living', '1': 'brain-dead', '2': 'DCD (nhb)'}),
    timepoint=s['timepoint biopsy'].map({'1': 'T1 donation (before retrieval)', '2': 'T2 after cold ischaemia', '3': 'T3 60 min after reperfusion'}),
    pair_code=s['deceased paired biopsies'], dgf=s['delayed graft function (0=no, 1=yes)'], title=s.title, age=np.nan, sex=np.nan, PMI=np.nan, pH=np.nan, RIN=np.nan,
    time_of_death='not recorded', cold_ischaemia_time='not recorded', platform=s.platform, data_type='Illumina HT-12 v4 (series matrix, normalised)', batch_fields='chip/position in description'))
save(s, 'transplant_GSE43974_samples.csv')
tab = pd.crosstab([s.group, s.timepoint], s.pair_code); tab.columns = ['pair_code_' + c for c in tab.columns]
tab.insert(0, 'dataset', 'GSE43974'); designs.append(tab.reset_index().rename(columns={'timepoint': 'region'}))
x = collapse(sm_table(T), probe2sym('GPL10558', 'Symbol')); clocks.append(clock_summary(x, 'GSE43974', False))

pd.concat(designs[:4]).to_csv(f'{OUT}/psych_design_tables.csv', index=False)
designs[4].to_csv(f'{OUT}/transplant_GSE43974_design.csv', index=False)
pd.concat(clocks).to_csv(f'{OUT}/clock_gene_detection.csv', index=False)
print(pd.concat(clocks).round(1).query("dataset in ['GSE53987','GSE43974']").to_string())
