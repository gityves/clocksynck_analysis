"""METABRIC feasibility scout: metadata and gene coverage only (no clock statistic). Run with python -I."""
import sys, gzip
import pandas as pd
D = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/metabric/brca_metabric/"
GMT = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/reference/h.all.v2024.1.Hs.symbols.gmt"
CLOCK = ['ARNTL','CLOCK','NPAS2','PER1','PER2','PER3','CRY1','CRY2','NR1D1','NR1D2','DBP','RORA','TEF','CIART']
PUR = ["COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRB","PTPRC","CD3E","CD14","LYZ","CD74"]
g2m = [g for l in open(GMT) for f in [l.rstrip("\n").split("\t")] if f[0]=="HALLMARK_G2M_CHECKPOINT" for g in f[2:] if g not in CLOCK]
ex = pd.read_csv(D+"data_mrna_illumina_microarray.txt", sep="\t", index_col=0)
print("expr shape genes x samples", ex.shape, "dup symbols", ex.index.duplicated().sum())
ent = ex.pop("Entrez_Gene_Id")
v = ex.to_numpy(float)
print("value range", float(pd.DataFrame(v).min().min()), float(pd.DataFrame(v).max().max()), "NaN frac", float(pd.isna(v).mean()))
print("sample mean-of-means / sd-of-means (non-z check):", float(ex.mean().mean()), float(ex.mean().std()), "gene-mean sd", float(ex.mean(axis=1).std()))
for name, gs in [("clock", CLOCK), ("G2M(no clock)", g2m), ("purity", PUR)]:
    miss = [g for g in gs if g not in ex.index]
    print(name, f"{len(gs)-len(miss)}/{len(gs)} present; missing:", miss)
alias = {"ARNTL":"BMAL1","CIART":"CHRONO","NPAS2":"","NR1D1":"REV-ERBa"}
print("alias check:", {a: (a in ex.index) for a in ["BMAL1","CHRONO","C1orf51"]})
for name, gs in [("clock", CLOCK), ("G2M", g2m), ("purity", PUR)]:
    d = [g for g in gs if g in ex.index and ex.index.get_loc(g) is not None and (ex.index==g).sum()>1]
    print(name, "symbols with >1 row (probe duplicates):", d)
print("clock genes: n rows, NaN count:", {g: (int((ex.index==g).sum()), int(ex.loc[g].isna().sum().sum())) for g in CLOCK})
print("cells NaN total", int(ex.isna().sum().sum()), "samples with NaN", int(ex.isna().any().sum()))
pt = pd.read_csv(D+"data_clinical_patient.txt", sep="\t", comment="#")
sm = pd.read_csv(D+"data_clinical_sample.txt", sep="\t", comment="#")
print("patients", len(pt), "samples", len(sm), "sample types", sm.SAMPLE_TYPE.value_counts().to_dict(), "patients with >1 sample", int(sm.PATIENT_ID.duplicated().sum()))
cl = sm.merge(pt, on="PATIENT_ID")
print("sample id == patient id for all:", bool((sm.SAMPLE_ID==sm.PATIENT_ID).all()))
mrna = [l.strip() for l in open(D+"case_lists/cases_mRNA.txt") if l.startswith("case_list_ids")][0].split(":")[1].split("\t")
print("cases_mRNA list", len(mrna), "expr cols", ex.shape[1], "cols in clinical", int(ex.columns.isin(cl.SAMPLE_ID).sum()))
c = cl[cl.SAMPLE_ID.isin(ex.columns)].copy()
print("N tumours with expression and clinical:", len(c))
for col in ["ER_IHC","ER_STATUS","HER2_SNP6","HER2_STATUS","PR_STATUS","CLAUDIN_SUBTYPE","COHORT","THREEGENE","INTCLUST","CHEMOTHERAPY","HORMONE_THERAPY","RADIO_THERAPY","OS_STATUS","VITAL_STATUS","RFS_STATUS","GRADE","CELLULARITY","INFERRED_MENOPAUSAL_STATE","SEX"]:
    print("\n##", col); print(c[col].value_counts(dropna=False).to_string())
for col in ["OS_MONTHS","RFS_MONTHS","AGE_AT_DIAGNOSIS"]:
    print(col, "non-missing", int(c[col].notna().sum()), "median", float(c[col].median()))
print("\n## ER_IHC vs ER_STATUS (sample-level)"); print(pd.crosstab(c.ER_IHC.fillna("NA"), c.ER_STATUS.fillna("NA")))
print("\n## PAM50 x ER_IHC"); print(pd.crosstab(c.CLAUDIN_SUBTYPE.fillna("NA"), c.ER_IHC.fillna("NA"), margins=True))
print("\n## PAM50 x COHORT"); print(pd.crosstab(c.CLAUDIN_SUBTYPE.fillna("NA"), c.COHORT.fillna("NA"), margins=True))
print("\n## ER_IHC x COHORT"); print(pd.crosstab(c.ER_IHC.fillna("NA"), c.COHORT.fillna("NA"), margins=True))
print("\n## ER x HER2_SNP6 (ER-)"); print(pd.crosstab(c.ER_IHC.fillna("NA"), c.HER2_SNP6.fillna("NA")))
print("\n## OS events (any death) / disease-specific:"); print(pd.crosstab(c.OS_STATUS.fillna("NA"), c.VITAL_STATUS.fillna("NA")))
print("\n## ER+ x hormone therapy x chemo"); print(pd.crosstab([c.ER_IHC.fillna("NA"), c.HORMONE_THERAPY.fillna("NA")], c.CHEMOTHERAPY.fillna("NA")))
