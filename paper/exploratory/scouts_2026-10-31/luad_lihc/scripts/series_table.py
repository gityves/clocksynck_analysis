"""Per-series scout table: counts from the header-derived sample tables (metadata only), platform coverage from platform_coverage.csv."""
from pathlib import Path
import sys
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import parse_hdr as P
R = Path(__file__).resolve().parents[1]; D = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_luad_lihc"
cov = pd.read_csv(R / "platform_coverage.csv").set_index("platform")
def T(name):
    return P.table(P.parse(R / "hdr" / f"{name}_series_matrix.hdr"))
def eq(t, col, vals): return int(t[col].isin(vals).sum())
rows = []
def add(series, cancer, nT, nN, plat, scale, ann, overlap, dl, hdrname=None, tag=""):
    c = cov.loc[plat]
    rows.append(dict(series=series, cancer=cancer, n_tumours=nT, n_normal=nN, platform=plat, scale=scale, n_clock_genes=f"{c.n_clock}/{c.n_clock_total}",
                     clock_missing=c.clock_missing if isinstance(c.clock_missing, str) else "", n_g2m_genes=f"{c.n_g2m}/{c.n_g2m_total}",
                     n_purity_genes=f"{c.n_purity}/{c.n_purity_total}", annotations=ann, overlap_notes=overlap, downloaded="yes" if dl else "no",
                     path=f"{D}/{series}/" if dl else ""))
t = T("GSE72094"); add("GSE72094", "LUAD", len(t), 0, "GPL15048", "log2-like (IRON-normalised, RNA-quality batch-corrected; min -0.31 max 15.0)",
    f"KRAS mut {eq(t,'kras_status',['Mut'])}, EGFR mut {eq(t,'egfr_status',['Mut'])}, STK11 mut {eq(t,'stk11_status',['Mut'])}, TP53 mut {eq(t,'tp53_status',['Mut'])}; smoking Ever {eq(t,'smoking_status',['Ever'])}/Never {eq(t,'smoking_status',['Never'])}/missing {eq(t,'smoking_status',['Missing'])}; stage; age; sex; survival; no TRU/PP/PI label (derivable from expression, not deposited); no adjacent normal",
    "Moffitt Cancer Center (Schabath/Eschrich); no TCGA samples expected (independent institutional cohort)", True)
t = T("GSE68465"); add("GSE68465", "LUAD", eq(t,'disease_state',['Lung Adenocarcinoma']), eq(t,'disease_state',['Normal']), "GPL96", "linear MAS5 (max 130623), log2 needed",
    "4 sites (MI 178, MSKCC 104, DFCI 82, HLM 79); 19 Stratagene 'Normal' arrays (reference RNA, not matched lung); stage, smoking history, grade, adjuvant chemo/RT, survival; no mutation calls",
    "Director's Challenge: multi-site batch structure (site as batch); no TCGA overlap expected; NOT in GSE72094/GSE31210", True)
t = T("GSE31210"); add("GSE31210", "LUAD", eq(t,'tissue',['primary lung tumor']), eq(t,'tissue',['normal lung']), "GPL570", "linear MAS5 (max 53814), log2 needed",
    "stage I-II only (IA/IB/II); EGFR mut 127, KRAS mut 20, ALK fusion 11, triple-negative 68; smoking ever/never 123/123; Brinkman index; relapse/death",
    "NCC Japan; single institution; early stage only; 20 normal lung (15 labelled as patient-linked, 'for KHxxx')", True)
t = T("GSE41271"); add("GSE41271", "LUAD", eq(t,'histology',['Adenocarcinoma']), 0, "GPL6884", "log2-like (Illumina WG6-v3; min 1.5 max 15.2), 275 NSCLC of which 80 squamous",
    "histology (183 ADC of 275), stage IA-IV, smoking Y/N, race, sex, recurrence; no mutation calls; no normals",
    "MD Anderson (Wistuba/Girard/Heymach); subset of the MDACC cohort; no TCGA overlap expected", True)
t = T("GSE30219"); add("GSE30219", "LUAD", eq(t,'histology',['ADC']), eq(t,'histology',['NTL']), "GPL570", "not downloaded (RMA expected, to be verified)",
    "NSCLC+NE tumours (293 total; ADC only 85); stage; relapse; 14 non-tumoural lung", "CIT programme (France); <150 ADC", False)
t = T("GSE50081"); add("GSE50081", "LUAD", eq(t,'histology',['adenocarcinoma','NSCLC-favor adenocarcinoma']), 0, "GPL570", "not downloaded",
    "stage I-II (JBR.10-related UHN181); smoking Never 24/Ex 79/Current 57; 127 ADC + 1 'favor ADC' of 181 NSCLC", "UHN Toronto; <150 ADC", False)
t = T("GSE11969"); add("GSE11969", "LUAD", eq(t,'Histology',['AD']), eq(t,'Histology',['Normal lung tissue']), "GPL7015", "not downloaded (Agilent 21.6K custom, 12/14 clock genes)",
    "AD 90 of 163; EGFR/KRAS/TP53 status; TRU-a/TRU-b/non-TRU labels; BI smoking; stage", "Japan; <150 ADC; platform lacks NR1D2 and CIART", False)
# LIHC
t = T("GSE14520-GPL3921"); t2 = T("GSE14520-GPL571")
nT = eq(t,'Tissue',['Liver Tumor Tissue']) + eq(t2,'Tissue',['Liver Tumor Tissue']); nN = eq(t,'Tissue',['Liver Non-Tumor Tissue']) + eq(t2,'Tissue',['Liver Non-Tumor Tissue']) + 2
add("GSE14520", "LIHC", nT, nN, "GPL3921", "log2-like (MAS5->log2, min 2.5 max 14.2), two platforms: GPL3921 (225 T/220 NT) + GPL571 (22 T/19 NT/2 healthy); gene coverage shown for GPL3921 (GPL571: 13/14, 192/200 G2M)",
    "HBV viral status (HBV-related throughout: CC 160, AVR-CC 58), cirrhosis Y 223/N 19, BCLC, TNM, CLIP, AFP, tumour size, multinodular, survival/recurrence; paired non-tumour (239); 2 healthy donor livers; CGH survival group",
    "Chinese Liver Cancer Group (Fudan Liver Cancer Institute, Roessler 2010); two array batches; some overlap with other HBV Chinese cohorts not checkable from metadata; no TCGA overlap expected", True)
t = T("GSE25097"); add("GSE25097", "LIHC", eq(t,'tissue',['tumor liver']), eq(t,'tissue',['non_tumor liver','healthy liver']), "GPL10687", "NOT log2-like (min 0.005, median 0.11, max 66): scale to be established before use (ratio/linear?)",
    "tumour 268, adjacent non-tumour 243, cirrhotic 40 (not HCC-paired), healthy 6; no clinical fields in GEO", "Origin not stated in GEO record; check publications; no TCGA overlap expected", True)
t = T("GSE36376"); add("GSE36376", "LIHC", eq(t,'tissue',['liver tumor']), eq(t,'tissue',['adjacent non-tumor liver']), "GPL10558", "log2-like (Illumina HT-12v4, min 5.15 max 14.8; 0.002% NaN)",
    "AJCC stage (1 184, 2 179, 3a 44, 3b 17, 3c 9 over 433 arrays); adjacent non-tumour; aetiology not in GEO fields (Samsung Medical Center, Korea; HBV-endemic, per publication)",
    "Samsung Medical Center; training + validation sets in one series; possible overlap with GSE39791 (Korean, Baylor/Lee) not verifiable", True)
t = T("GSE63898"); add("GSE63898", "LIHC", eq(t,'tissue subtype',['primary hepatocellular carcinoma']), eq(t,'tissue subtype',['non-tumor cirrhosis']), "GPL13667", "log2-like (RMA, min 3.0 max 14.6)",
    "sex (tumours 180 M/47 F), BCLC (A 178, 0 19, B 17, C 11; clinical fields only on the 228 tumours); comparison tissue is cirrhosis (168), not paired normal; aetiology mostly HCV per publication (not in GEO fields)",
    "Mount Sinai/Barcelona (Villanueva, Llovet); publication PMID 31344396 shared with Hoshida GSE10141/10143 meta-analysis; possible patient overlap with GSE10143 not verifiable", True)
t = T("GSE76427"); add("GSE76427", "LIHC", eq(t,'tissue',['primary hepatocellular carcinoma tumor']), eq(t,'tissue',['adjacent non-tumor liver tissue']), "GPL10558", "not downloaded (lumi RSN, log2 expected)",
    "BCLC, TNM, OS/RFS, age/sex; 52 adjacent non-tumour; Singapore (NCCS)", "<150 tumours", False)
t = T("GSE54236"); add("GSE54236", "LIHC", eq(t,'tissue type',['Biopsy of tumor tissue']), eq(t,'tissue type',['Biopsy of cirrhotic non-malignant  tissue']), "GPL6480", "not downloaded (Agilent, one-colour; two labelling kits)",
    "tumour doubling time, survival, sex; HCV-cirrhosis (Italy, biopsies); labelling-kit batches", "<150 tumours; multiple kits", False)
t = T("GSE10143"); add("GSE10143", "LIHC", 80, 307, "GPL5474", "not downloaded (DASL, 6k panel)",
    "80 HCC tumours + 307 non-tumour; HBV/HCV/alcohol in GSE10141 header; DASL panel lacks PER3 CRY1 NR1D1 TEF CIART (9/14 clock genes)", "Mount Sinai/Japan (Hoshida); GSE10141 is the 80-tumour training subset of the same series; unsuitable platform", False)
t = T("GSE64041"); add("GSE64041", "LIHC", eq(t,'tissue',['tumor']), eq(t,'tissue',['liver parenchyma','liver']), "GPL6244", "not downloaded (RMA)", "paired tumour/parenchyma biopsies (Basel), 5 normal livers", "<150 tumours", False)
t = T("GSE39791"); add("GSE39791", "LIHC", eq(t,'tissue',['tumor']), eq(t,'tissue',['non-tumor']), "GPL10558", "not downloaded", "paired non-tumour, age, sex", "<150 tumours", False)
df = pd.DataFrame(rows); df.to_csv(R / "series_table.csv", index=False)
pd.set_option("display.width", 250); print(df.drop(columns=["annotations", "overlap_notes", "path", "scale", "clock_missing"]).to_string())
