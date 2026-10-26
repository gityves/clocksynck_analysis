"""Inventory of GSE266919 (Zhang 2025) processed SCE .rds files: parse only (pure python, no R).
Usage: python -I inventory_GSE266919.py <raw_dir> <out_dir>. Writes output/inventory_tables/*.tsv and summary.json."""
import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rds_parse import read_rds, Obj, BigVec
import numpy as np, pandas as pd

raw, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
named = lambda x: dict(zip(x['names'].data, x.data))
files = ['Bcell', 'CD4Tcell', 'CD8Tcell', 'Myeloid', 'NKcell']
summ = {}
for name in files:
    o = read_rds(f'{raw}/GSE266919_{name}.rds.gz')
    cd = o['colData']
    cols = named(cd['listData'])
    df = {}
    for k, v in cols.items():
        d = v.data
        if v.kind == 13 and 'levels' in v:
            lv = v['levels'].data
            df[k] = [lv[i - 1] if i > 0 else None for i in d]
        else:
            df[k] = list(d)
    df = pd.DataFrame(df)
    cnt = named(o['assays']['data']['listData'])
    s = {'n_cells': int(cd['nrows'].data[0]), 'columns': list(df.columns), 'assays': list(cnt)}
    for an, a in cnt.items():
        dim = a['dim'].data
        dn = a['dimnames'].data
        s[f'assay_{an}'] = {'dim': [int(x) for x in dim], 'gene_example': list(dn[0].data[:5]) if dn[0] is not None else None,
                            'colnames_match_barcode': bool(dn[1] is not None and list(dn[1].data[:1000]) == list(df['Barcode'][:1000])),
                            'type': a.kind}
        d = a.data
        if isinstance(d, BigVec):
            s[f'assay_{an}'].update(integer_valued=d.integer_valued, min=d.vmin, max=d.vmax, head=[float(x) for x in d.head[:10]])
        else:
            dd = np.asarray(d); s[f'assay_{an}'].update(integer_valued=bool(np.all(dd == np.round(dd))), min=float(dd.min()), max=float(dd.max()))
        if dn[0] is not None:
            pd.Series(list(dn[0].data)).to_csv(f'{out}/{name}_genes.txt', index=False, header=False)
    s['value_counts'] = {c: df[c].value_counts(dropna=False).head(60).to_dict() for c in df.columns
                         if df[c].dtype == object and df[c].nunique() < 400 and c != 'Barcode'}
    s['n_samples'] = int(df['Sample'].nunique()); s['n_patients'] = int(df['Patient'].nunique())
    summ[name] = s
    for lev in ['MajorCluster', 'SubCluster', 'Subset']:
        if lev in df:
            ct = pd.crosstab(df['Sample'], df[lev]); ct.to_csv(f'{out}/{name}_sample_by_{lev}.tsv', sep='\t')
    meta = df[['Sample', 'Patient', 'Group', 'Tissue', 'Treatment', 'Efficacy', 'Response']].drop_duplicates()
    meta.to_csv(f'{out}/{name}_sample_meta.tsv', sep='\t', index=False)
    df.drop(columns=['Barcode']).groupby('Sample').size().rename('n_cells').to_csv(f'{out}/{name}_cells_per_sample.tsv', sep='\t')
    print(name, s['n_cells'], s['assay_counts']['dim'], s['assay_counts'].get('integer_valued'), flush=True)
    del o
json.dump(summ, open(f'{out}/summary.json', 'w'), indent=1, default=str)
