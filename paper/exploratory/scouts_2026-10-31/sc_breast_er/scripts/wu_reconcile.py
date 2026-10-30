"""Metadata-only: reconcile Wu 2021 ER/HER2/subtype per patient: paper Suppl. Table 1 (Nat Genet ESM4 xlsx),
GEO clinical_subtype and metadata.csv 'subtype'. Usage: python -I wu_reconcile.py MOESM4.xlsx wu_patient_x_class_ER.csv out.csv"""
import sys, openpyxl, pandas as pd
ws = openpyxl.load_workbook(sys.argv[1], read_only=True)['Supplementary Table 1']
rows = list(ws.iter_rows(values_only=True))
hdr = rows[3]
p = pd.DataFrame([r for r in rows[4:] if r[0] is not None], columns=hdr)
# cell-table ids drop the hyphen/suffix of the paper (CID4404-1 -> CID44041, CID4290 -> CID4290A, CID4530 -> CID4530N)
IDMAP = {'CID4404-1': 'CID44041', 'CID4497-1': 'CID44971', 'CID4499-1': 'CID44991', 'CID4517-1': 'CID45171', 'CID4290': 'CID4290A', 'CID4530': 'CID4530N'}
p['orig.ident'] = ('CID' + p['Case ID'].astype(str)).replace(IDMAP)
p = p.set_index('orig.ident')
t = pd.read_csv(sys.argv[2], index_col=0)
o = pd.DataFrame({'paper_ER_raw': p['ER'].astype(str), 'paper_PR': p['PR'].astype(str), 'paper_HER2_IHC': p['HER2 IHC'].astype(str),
                  'paper_HER2_ISH': p['HER2 ISH (ratio)'], 'paper_subtype_IHC': p['Subtype by IHC'], 'treatment': p['Treatment status'],
                  'cancer_type': p['Cancer Type']})
o['ER_paper'] = o.paper_subtype_IHC.map(lambda s: 'ER+' if 'ER+' in s else 'ER-')
o['GEO_clinical_subtype'] = t['clinical_subtype']; o['metadata_subtype'] = t['subtype']
o['ER_GEO'] = t['ER']
o['final_label'] = o['paper_subtype_IHC']; o['final_ER'] = o['ER_paper']; o['source'] = 'Wu 2021 Suppl Table 1 (Subtype by IHC)'
cls = ['Cancer Epithelial', 'CAFs', 'Endothelial', 'PVL', 'Normal Epithelial']
o = o.join(t[cls])
o.to_csv(sys.argv[3])
print('n patients', len(o), 'in paper not in cells:', set(p.index) - set(t.index), 'in cells not in paper:', set(t.index) - set(p.index))
print('ER paper vs GEO disagreements:', list(o.index[o.ER_paper != o.ER_GEO]))
print('HER2/ER label disagreements paper vs metadata subtype:')
print(o[o.paper_subtype_IHC.str.replace('HER2+/ER+', 'ER+').ne(o.metadata_subtype.replace({'ER+': 'ER+'}))][['paper_ER_raw', 'paper_subtype_IHC', 'GEO_clinical_subtype', 'metadata_subtype']])
print(o.groupby('final_ER').size())
