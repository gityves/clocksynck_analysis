"""Pooled per-patient table (Wu 2021 + Pal 2021) and per-ER-arm summary. Metadata only (counts of cells per class).
Pal: only 'tumour' (R-code Tum clusters, inferCNV-defined) and 'microenv' (non-tumour, non-epithelial clusters, unlabelled) exist."""
import pandas as pd
w = pd.read_csv('wu_reconciled.csv', index_col=0)
p = pd.read_csv('pal_sample_tumour_microenv.csv')
w2 = pd.DataFrame({'cohort': 'Wu', 'patient': w.index, 'label': w.final_label.values, 'ER_arm': w.final_ER.values, 'malignant': w['Cancer Epithelial'].values,
                   'CAF': w.CAFs.values, 'endothelial': w.Endothelial.values, 'PVL': w.PVL.values, 'microenv_unlabelled': pd.NA,
                   'ER_source': 'Wu Suppl Table 1 (IHC)'})
def lab(s):
    if s.startswith('ER_'): return 'ER+' if s != 'ER_0319' else 'ER-low/PR+'
    if s.startswith('HER2'): return 'HER2+ (ER not reported)'
    return 'TNBC (BRCA1)' if s.startswith('TN_B1') else 'TNBC'
arm = lambda l: {'ER+': 'ER+', 'ER-low/PR+': 'ER+', 'TNBC': 'ER-', 'TNBC (BRCA1)': 'ER-'}.get(l, 'unknown')
p2 = pd.DataFrame({'cohort': 'Pal', 'patient': p['sample'], 'label': p['sample'].map(lab), 'malignant': p.tumour, 'CAF': pd.NA, 'endothelial': pd.NA,
                   'PVL': pd.NA, 'microenv_unlabelled': p.microenv,
                   'ER_source': 'GEO title + Table EV3/EV2 (HER2: ER not reported)'})
p2['ER_arm'] = p2.label.map(arm)
pool = pd.concat([w2, p2], ignore_index=True); pool.to_csv('pooled_patient_table.csv', index=False)
print(pool.to_string())
print('\nWU per arm'); 
for a, g in w2.groupby('ER_arm'):
    print(a, len(g), 'malig>=50', int((g.malignant >= 50).sum()), 'CAF>=50', int((g.CAF >= 50).sum()), 'endo>=50', int((g.endothelial >= 50).sum()),
          'all3>=50', int(((g.malignant >= 50) & (g.CAF >= 50) & (g.endothelial >= 50)).sum()), 'all3>=30', int(((g.malignant >= 30) & (g.CAF >= 30) & (g.endothelial >= 30)).sum()))
print('\nPAL per arm (samples; unique patients)')
for a, g in p2.groupby('ER_arm'):
    print(a, len(g), 'tumour>=50', int((g.malignant >= 50).sum()), 'tumour>=50 & microenv>=50', int(((g.malignant >= 50) & (g.microenv_unlabelled >= 50)).sum()),
          'tumour>=30 & microenv>=30', int(((g.malignant >= 30) & (g.microenv_unlabelled >= 30)).sum()), 'median tumour', int(g.malignant.median()), 'median microenv', int(g.microenv_unlabelled.median()))
